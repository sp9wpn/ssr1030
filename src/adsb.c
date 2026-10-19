/* 1090 MHz Mode S downlink decoder used during short ADS-B windows to build
 * an address -> callsign table. Sources:
 *   DF17/18 TC 1-4 (identification), CRC must be zero;
 *   DF11 all-call replies (address presence), PI syndrome within II/SI range;
 *   DF20/21 Comm-B with BDS 2,0 (identification), AP syndrome = address,
 *   accepted only for addresses already known. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ssr1030.h"

#define CS_HASH 16384

struct cs_entry {
    uint32_t addr;
    char cs[9];
    double last;
    bool used;
};
static struct cs_entry tab[CS_HASH];
static int ncs, nknown;
static unsigned cs_ver;

static struct cs_entry *cs_slot(uint32_t a, bool create)
{
    uint32_t h = (a * 2654435761u) >> 18;
    for (int k = 0; k < CS_HASH; k++) {
        struct cs_entry *e = &tab[(h + k) % CS_HASH];
        if (!e->used) {
            if (!create || nknown > CS_HASH * 3 / 4) return NULL;
            e->used = true;
            e->addr = a;
            nknown++;
            return e;
        }
        if (e->addr == a) return e;
    }
    return NULL;
}

const char *cs_lookup(uint32_t a)
{
    struct cs_entry *e = cs_slot(a, false);
    return e && e->cs[0] ? e->cs : NULL;
}

/* display name: callsign when known, else the hex address */
const char *ac_name(uint32_t a)
{
    static char bufs[8][9];
    static int n;
    const char *c = cs_lookup(a);
    if (c) return c;
    char *b = bufs[n++ % 8];
    snprintf(b, 9, "%06X", a & 0xFFFFFF);
    return b;
}

static void trim(char *s)
{
    int n = strlen(s);
    while (n > 0 && s[n - 1] == ' ') s[--n] = 0;
}

int64_t cs_find(const char *callsign)
{
    char want[9];
    snprintf(want, sizeof(want), "%s", callsign);
    for (char *p = want; *p; p++) if (*p >= 'a' && *p <= 'z') *p -= 32;
    trim(want);
    for (int k = 0; k < CS_HASH; k++)
        if (tab[k].used && tab[k].cs[0] && !strcmp(tab[k].cs, want)) return tab[k].addr;
    return -1;
}

int cs_count(void) { return ncs; }
unsigned cs_version(void) { return cs_ver; }

static void cs_set(uint32_t a, const char *cs, double when)
{
    struct cs_entry *e = cs_slot(a, true);
    if (!e) return;
    if (strcmp(e->cs, cs)) {
        if (!e->cs[0]) ncs++;
        snprintf(e->cs, sizeof(e->cs), "%s", cs);
        cs_ver++;
    }
    e->last = when;
}

static bool known(uint32_t a)
{
    return cs_slot(a, false) != NULL || analysis_addr_known(a);
}

int cs_load(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char line[128];
    int n = 0;
    while (fgets(line, sizeof(line), f)) {
        unsigned a; char cs[16]; double t = 0;
        if (sscanf(line, "%x %15s %lf", &a, cs, &t) >= 2) {
            cs[8] = 0;
            cs_set(a & 0xFFFFFF, cs, t);
            n++;
        }
    }
    fclose(f);
    return n;
}

int cs_save(const char *path)
{
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) return -1;
    for (int k = 0; k < CS_HASH; k++)
        if (tab[k].used && tab[k].cs[0]) fprintf(f, "%06X %s %.0f\n", tab[k].addr, tab[k].cs, tab[k].last);
    fclose(f);
    return rename(tmp, path);
}

/* ---------------- demodulator ---------------- */

static double fs, spus;
static float *mag;                 /* tail + current block */
static size_t mcap, mlen, resume;
static int o_hi[4], o_lo[6];
static struct adsb_stats st;

void adsb_init(double rate)
{
    fs = rate;
    spus = fs / 1e6;
    const double hi[4] = { 0.25, 1.25, 3.75, 4.75 };
    const double lo[6] = { 0.75, 2.25, 2.75, 3.25, 5.75, 6.75 };
    for (int k = 0; k < 4; k++) o_hi[k] = (int)lround(hi[k] * spus);
    for (int k = 0; k < 6; k++) o_lo[k] = (int)lround(lo[k] * spus);
}

void adsb_window_start(void)
{
    mlen = 0;
    resume = 0;
    st.windows++;
    st.win_frames = st.win_clip = st.win_samples = 0;
    st.win_peak_db = -99;
}

void adsb_get_stats(struct adsb_stats *s) { *s = st; }

static const char charset[65] = "#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######";

static bool callsign(const uint8_t *msg, int first, char out[9])
{
    int letters = 0;
    for (int k = 0; k < 8; k++) {
        char c = charset[getbits(msg, first + 6 * k, first + 6 * k + 5)];
        if (c == '#') return false;
        if (c != ' ') letters++;
        out[k] = c;
    }
    out[8] = 0;
    trim(out);
    return letters > 0 && out[0] != ' ';
}

static double half_bit(const float *m, double start, double len)
{
    int a = (int)lround(start + 1), b = (int)lround(start + len - 1);
    double s = 0;
    for (int k = a; k < b; k++) s += m[k];
    return s / (b - a > 0 ? b - a : 1);
}

/* returns message length in samples when decoded, 0 otherwise */
static int try_decode(const float *m, float thr)
{
    uint8_t msg[14] = { 0 };
    double base = 8.0 * spus, bl = spus, hb = 0.5 * spus;
    int nbits = 112;
    for (int i = 0; i < nbits; i++) {
        double a = half_bit(m, base + i * bl, hb), b = half_bit(m, base + i * bl + hb, hb);
        if (a < thr * 0.5 && b < thr * 0.5) return 0;     /* signal dropped */
        if (a > b) msg[i / 8] |= 0x80 >> (i % 8);
        if (i == 4) {
            int df = msg[0] >> 3;
            if (df != 11 && df != 17 && df != 18 && df != 20 && df != 21) return 0;
            nbits = df >= 16 ? 112 : 56;
        }
    }
    int df = msg[0] >> 3;
    uint32_t syn = modes_crc(msg, nbits) ^ getbits(msg, nbits - 23, nbits);
    double now = (double)time(NULL);
    char cs[9];

    if (df == 17 || df == 18) {
        if (syn) return 0;
        st.frames++; st.df17++;
        uint32_t aa = getbits(msg, 9, 32);
        cs_slot(aa, true);
        int tc = getbits(msg, 33, 37);
        if (tc >= 1 && tc <= 4 && callsign(msg, 41, cs)) { cs_set(aa, cs, now); st.idents++; }
    } else if (df == 11) {
        if (syn & ~0x7Fu) return 0;
        st.frames++; st.df11++;
        cs_slot(getbits(msg, 9, 32), true);
    } else {                                          /* DF20/21 */
        uint32_t aa = syn;
        if (!known(aa)) return 0;
        st.frames++;
        if (getbits(msg, 33, 40) == 0x20 && callsign(msg, 41, cs)) { cs_set(aa, cs, now); st.commb_idents++; }
    }
    return (int)lround((8 + nbits) * spus);
}

void adsb_process(const int8_t *iq, size_t ns)
{
    size_t need = mlen + ns;
    if (need > mcap) {
        mcap = need + 65536;
        mag = realloc(mag, mcap * sizeof(float));
    }
    double si = 0, sq = 0;
    for (size_t k = 0; k < ns; k++) { si += iq[2 * k]; sq += iq[2 * k + 1]; }
    float di = (float)(si / ns), dq = (float)(sq / ns);
    double acc = 0;
    for (size_t k = 0; k < ns; k++) {
        if (iq[2 * k] >= 127 || iq[2 * k] <= -128 || iq[2 * k + 1] >= 127 || iq[2 * k + 1] <= -128) st.win_clip++;
        float i = iq[2 * k] - di, q = iq[2 * k + 1] - dq;
        float a = sqrtf(i * i + q * q);
        mag[mlen + k] = a;
        acc += a;
    }
    mlen += ns;
    st.seconds += ns / fs;
    st.win_samples += ns;
    float thr = (float)(acc / ns) * 4.0f;           /* ~12 dB over mean */

    size_t msg_span = (size_t)((8 + 112 + 2) * spus);
    if (mlen < msg_span) return;
    size_t end = mlen - msg_span, p = resume;
    while (p < end) {
        const float *m = mag + p;
        if (m[o_hi[0]] < thr || m[o_hi[1]] < thr || m[o_hi[2]] < thr || m[o_hi[3]] < thr) { p++; continue; }
        float h = 0, l = 0;
        for (int k = 0; k < 4; k++) h += m[o_hi[k]];
        for (int k = 0; k < 6; k++) l += m[o_lo[k]];
        h /= 4; l /= 6;
        if (h < 3 * l) { p++; continue; }
        /* refine: best of the next quarter microsecond */
        size_t best = p; float bs = h - l;
        for (size_t q = p + 1; q < p + (size_t)(0.25 * spus) && q < end; q++) {
            const float *mq = mag + q;
            float hq = (mq[o_hi[0]] + mq[o_hi[1]] + mq[o_hi[2]] + mq[o_hi[3]]) / 4;
            float lq = 0;
            for (int k = 0; k < 6; k++) lq += mq[o_lo[k]];
            lq /= 6;
            if (hq - lq > bs) { bs = hq - lq; best = q; }
        }
        int len = try_decode(mag + best, h);
        if (len) {
            st.win_frames++;
            float db = 20.0f * log10f(h / 127.0f);
            if (db > st.win_peak_db) st.win_peak_db = db;
        }
        p = len ? best + len : p + 1;
    }
    resume = p > end ? p - end : 0;
    memmove(mag, mag + end, (mlen - end) * sizeof(float));
    mlen -= end;
}

/* ---------------- adsb.lol bootstrap ---------------- */

/* Fetch aircraft within radius_nm of lat/lon from api.adsb.lol and add every
 * hex/callsign pair to the table. Uses the curl binary. Returns the number of
 * callsigns added, or -1 on failure. */
int cs_fetch_adsblol(double lat, double lon, int radius_nm)
{
    char cmd[256];
    snprintf(cmd, sizeof(cmd),
             "curl -s -f --max-time 20 'https://api.adsb.lol/v2/lat/%.5f/lon/%.5f/dist/%d'",
             lat, lon, radius_nm);
    FILE *p = popen(cmd, "r");
    if (!p) return -1;
    size_t cap = 1 << 16, len = 0;
    char *buf = malloc(cap);
    size_t r;
    while ((r = fread(buf + len, 1, cap - len - 1, p)) > 0) {
        len += r;
        if (cap - len < 4096) { cap *= 2; buf = realloc(buf, cap); }
    }
    int rc = pclose(p);
    buf[len] = 0;
    if (rc != 0 || !strstr(buf, "\"ac\"")) { free(buf); return -1; }

    /* each aircraft object starts at "hex"; its "flight" lies before the next "hex" */
    int added = 0;
    double now = (double)time(NULL);
    char *h = strstr(buf, "\"hex\":\"");
    while (h) {
        char *v = h + 7;
        char *next = strstr(v, "\"hex\":\"");
        unsigned a;
        char *f = strstr(v, "\"flight\":\"");
        if (*v != '~' && sscanf(v, "%6x", &a) == 1 && f && (!next || f < next)) {
            char cs[9] = { 0 };
            f += 10;
            int k = 0;
            while (k < 8 && f[k] && f[k] != '"') { cs[k] = f[k]; k++; }
            trim(cs);
            if (cs[0]) { cs_set(a & 0xFFFFFF, cs, now); added++; }
        }
        h = next;
    }
    free(buf);
    return added;
}
