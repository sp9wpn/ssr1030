/* bds_logger: reads a Beast stream (dump1090, readsb, ...) and prints the
 * Comm-B replies (DF20/21) of ADS-B aircraft that are BDS 5,0 (track and
 * turn) or BDS 6,0 (heading and speed), one line each.
 *
 * A Comm-B reply does not say which register it holds; that is in the
 * interrogation, which a 1090 receiver does not hear. So the register is
 * inferred from the content, as in pyModeS:
 *  - every field with its status bit clear must be all zeros, and the
 *    values must be in range (roll, speeds, rates, Mach);
 *  - the aircraft must be known from its own ADS-B (DF17) in the last minute,
 *    and the content must agree with its ADS-B velocity: BDS 5,0 ground
 *    speed and track, BDS 6,0 heading (wind correction) and Mach / IAS
 *    against its altitude and ground speed;
 *  - a reply that fits both registers is printed with both decodes, each
 *    marked '?'.
 * The address of a DF20/21 reply is its parity syndrome (AP field).
 *
 * Altitude: the last pressure altitude (DF20 reply, DF0/4/16, ADS-B), made
 * AMSL with --qnh, or else with a QNH read from BDS 4,0 replies (barometric
 * setting) of aircraft below 4000 ft that are not set to standard.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <math.h>
#include <netdb.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define AC_N        4096            /* aircraft table */
#define KNOWN_AGE   60.0            /* s: an aircraft is known this long after its last DF17 */
#define VEL_AGE     30.0            /* s: an ADS-B velocity is used this long */
#define ALT_AGE     60.0            /* s: a pressure altitude is shown this long */
#define QNH_N       16              /* barometric settings kept for the derived QNH */
#define QNH_AGE     1800.0          /* s */
#define QNH_MAX_FT  4000            /* settings only from aircraft below this (above the transition: standard) */

struct ac {
    uint32_t addr;
    bool used;
    double seen;                    /* last DF17 */
    char cs[9];
    double vel_t;                   /* ADS-B velocity: ground speed kt, track deg, vertical rate ft/min */
    double gs, trk, vr;
    bool vr_ok;
    double alt_t;
    int alt;                        /* pressure altitude, ft */
};

static struct ac act[AC_N];
static double qnh_cmd = NAN;        /* --qnh */
static int max_alt = INT32_MAX;     /* --max-alt: pressure altitude, ft */
#define SEL_N 32
static char sel[SEL_N][9];          /* -A: hex addresses or callsigns */
static int nsel;
static double qset[QNH_N], qset_t[QNH_N];
static int qset_n;
static bool verbose;
static long n_commb, n_unknown, n_50, n_60, n_both, n_none, n_high, n_other;

static double now_s(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static unsigned getbits(const uint8_t *m, int first, int last)
{
    unsigned v = 0;
    for (int b = first; b <= last; b++) v = v << 1 | ((m[(b - 1) >> 3] >> (7 - ((b - 1) & 7))) & 1);
    return v;
}

/* Mode S CRC-24 remainder over the whole message (0 for a good DF17, the address for DF20/21) */
static uint32_t crc(const uint8_t *m, int nbits)
{
    uint32_t r = 0;
    for (int b = 0; b < nbits - 24; b++) {
        uint32_t bit = (m[b >> 3] >> (7 - (b & 7))) & 1;
        uint32_t top = (r >> 23) & 1;
        r = (r << 1 | bit) & 0xFFFFFF;
        if (top) r ^= 0xFFF409;
    }
    for (int b = 0; b < 24; b++) {
        uint32_t top = (r >> 23) & 1;
        r = (r << 1) & 0xFFFFFF;
        if (top) r ^= 0xFFF409;
    }
    return r ^ getbits(m, nbits - 23, nbits);
}

static struct ac *ac_get(uint32_t a, bool create)
{
    uint32_t h = (a * 2654435761u) >> 20;
    for (int k = 0; k < AC_N; k++) {
        struct ac *e = &act[(h + k) % AC_N];
        if (!e->used) {
            if (!create) return NULL;
            memset(e, 0, sizeof(*e));
            e->used = true;
            e->addr = a;
            e->seen = e->vel_t = e->alt_t = -1e300;
            return e;
        }
        if (e->addr == a) return e;
    }
    return NULL;
}

/* ---------------- altitude ---------------- */

/* 13-bit field (C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4) -> 0xABCD, one octal digit per nibble */
static unsigned abcd13(unsigned v)
{
#define B13(n) ((v >> (13 - (n))) & 1)
    return B13(6) << 14 | B13(4) << 13 | B13(2) << 12 | B13(12) << 10 | B13(10) << 9 | B13(8) << 8 |
           B13(5) << 6 | B13(3) << 5 | B13(1) << 4 | B13(13) << 2 | B13(11) << 1 | B13(9);
#undef B13
}

/* Gillham code -> 100 ft units, INT32_MIN if invalid */
static int gillham(unsigned a)
{
    if ((a & 0xFFFF8889) || !(a & 0xF0)) return INT32_MIN;
    unsigned h = 0, f = 0;
    if (a & 0x10) h ^= 7;
    if (a & 0x20) h ^= 3;
    if (a & 0x40) h ^= 1;
    if ((h & 5) == 5) h ^= 2;
    if (h > 5) return INT32_MIN;
    if (a & 0x0002) f ^= 0xFF;
    if (a & 0x0004) f ^= 0x7F;
    if (a & 0x1000) f ^= 0x3F;
    if (a & 0x2000) f ^= 0x1F;
    if (a & 0x4000) f ^= 0x0F;
    if (a & 0x0100) f ^= 0x07;
    if (a & 0x0200) f ^= 0x03;
    if (a & 0x0400) f ^= 0x01;
    if (f & 1) h = 6 - h;
    return (int)(f * 5 + h) - 13;
}

/* 13-bit altitude code (DF0/4/16/20) -> feet, INT32_MIN if absent, metric or invalid */
static int ac13_ft(unsigned v)
{
    if (!v || (v & 0x40)) return INT32_MIN;
    if (v & 0x10) return (int)(((v & 0x1F80) >> 2 | (v & 0x20) >> 1 | (v & 0x0F)) * 25) - 1000;
    int h = gillham(abcd13(v));
    return h == INT32_MIN ? INT32_MIN : h * 100;
}

/* ADS-B 12-bit altitude: the 13-bit code without its M bit */
static int ac12_ft(unsigned v) { return ac13_ft((v & 0xFC0) << 1 | (v & 0x3F)); }

static void set_alt(struct ac *a, int ft, double t)
{
    if (ft == INT32_MIN) return;
    a->alt = ft;
    a->alt_t = t;
}

/* QNH: --qnh, else the median of recent BDS 4,0 settings; NAN when none */
static double qnh_now(double t)
{
    if (!isnan(qnh_cmd)) return qnh_cmd;
    double v[QNH_N];
    int n = 0;
    for (int k = 0; k < qset_n && k < QNH_N; k++)
        if (t - qset_t[k] < QNH_AGE) v[n++] = qset[k];
    if (!n) return NAN;
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && v[j] < v[j - 1]; j--) { double x = v[j]; v[j] = v[j - 1]; v[j - 1] = x; }
    return n & 1 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* pressure altitude -> AMSL with the ISA pressure-height relation */
static double amsl_ft(int pa, double qnh)
{
    double h_qnh = 145366.45 * (1 - pow(qnh / 1013.25, 0.190284));   /* pressure altitude of the QNH level */
    return pa - h_qnh;
}

/* ---------------- ADS-B (DF17) ---------------- */

static const char ais[] = "#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######";

static void df17(const uint8_t *m, double t)
{
    struct ac *a = ac_get(getbits(m, 9, 32), true);
    if (!a) return;
    a->seen = t;
    int tc = getbits(m, 33, 37);
    if (tc >= 1 && tc <= 4) {
        for (int k = 0; k < 8; k++) a->cs[k] = ais[getbits(m, 41 + 6 * k, 46 + 6 * k)];
        a->cs[8] = 0;
        for (int k = 7; k >= 0 && (a->cs[k] == ' ' || a->cs[k] == '#'); k--) a->cs[k] = 0;
    } else if (tc >= 9 && tc <= 18) {
        set_alt(a, ac12_ft(getbits(m, 41, 52)), t);
    } else if (tc == 19) {
        int st = getbits(m, 38, 40);
        if (st == 1 || st == 2) {
            int ew = getbits(m, 47, 56), ns = getbits(m, 58, 67);
            if (!ew || !ns) return;
            double k = st == 2 ? 4 : 1;
            double vx = (getbits(m, 46, 46) ? -1 : 1) * (ew - 1) * k;
            double vy = (getbits(m, 57, 57) ? -1 : 1) * (ns - 1) * k;
            a->gs = hypot(vx, vy);
            a->trk = fmod(atan2(vx, vy) * 180 / M_PI + 360, 360);
            int vr = getbits(m, 70, 78);
            a->vr_ok = vr != 0;
            a->vr = (getbits(m, 69, 69) ? -1 : 1) * (vr - 1) * 64.0;
            a->vel_t = t;
        }
    }
}

/* ---------------- Comm-B inference ---------------- */

/* signed field with its sign bit just before it */
static double sgn(const uint8_t *m, int sbit, int first, int last)
{
    int v = (int)getbits(m, first, last);
    return getbits(m, sbit, sbit) ? v - (1 << (last - first + 1)) : v;
}

/* a field whose status bit is clear must be all zeros */
static bool field_ok(const uint8_t *m, int status, int last)
{
    return getbits(m, status, status) || !getbits(m, status + 1, last);
}

static double angdiff(double a, double b)
{
    double d = fmod(fabs(a - b), 360);
    return d > 180 ? 360 - d : d;
}

struct b50 { bool roll_ok, trk_ok, gs_ok, rate_ok, tas_ok; double roll, trk, gs, rate, tas; };
struct b60 { bool hdg_ok, ias_ok, mach_ok, baro_ok, ins_ok; double hdg, ias, mach, baro, ins; };

/* -A lists this address as hex: its replies are decoded even without ADS-B */
static bool hex_given(uint32_t addr)
{
    for (int k = 0; k < nsel; k++) {
        char *e;
        unsigned long v = strtoul(sel[k], &e, 16);
        if (strlen(sel[k]) == 6 && !*e && v == addr) return true;
    }
    return false;
}

/* ref: an ADS-B velocity to check against. Without one (an aircraft given as
 * hex with -A and no recent ADS-B), every field must be present instead */
static bool is50(const uint8_t *m, const struct ac *a, bool ref, struct b50 *o)
{
    if (!field_ok(m, 33, 43) || !field_ok(m, 44, 55) || !field_ok(m, 56, 66) || !field_ok(m, 67, 77) || !field_ok(m, 78, 88))
        return false;
    o->roll_ok = getbits(m, 33, 33); o->roll = sgn(m, 34, 35, 43) * 45.0 / 256;
    o->trk_ok = getbits(m, 44, 44);  o->trk = fmod(sgn(m, 45, 46, 55) * 90.0 / 512 + 360, 360);
    o->gs_ok = getbits(m, 56, 56);   o->gs = getbits(m, 57, 66) * 2.0;
    o->rate_ok = getbits(m, 67, 67); o->rate = sgn(m, 68, 69, 77) * 8.0 / 256;
    o->tas_ok = getbits(m, 78, 78);  o->tas = getbits(m, 79, 88) * 2.0;
    if (!o->trk_ok || !o->gs_ok) return false;                  /* the fields checked against ADS-B */
    if (o->roll_ok && fabs(o->roll) > 50) return false;
    if (o->gs < 40 || o->gs > 700) return false;
    if (o->tas_ok && (o->tas < 40 || o->tas > 600 || fabs(o->gs - o->tas) > 200)) return false;
    if (o->rate_ok && fabs(o->rate) > 8) return false;
    if (!ref) return o->roll_ok && o->tas_ok;
    return fabs(o->gs - a->gs) <= 25 && angdiff(o->trk, a->trk) <= 12;
}

/* Mach from calibrated airspeed (kt) at pressure altitude (ft), ISA */
static double cas_to_mach(double cas, double alt_ft)
{
    double p = 101325 * pow(1 - 6.8756e-6 * alt_ft, 5.2559);
    if (alt_ft > 36089) p = 22632 * exp(-4.8063e-5 * (alt_ft - 36089));
    double qc = 101325 * (pow(1 + 0.2 * pow(cas / 661.47, 2), 3.5) - 1);
    return sqrt(5 * (pow(qc / p + 1, 2.0 / 7) - 1));
}

/* ISA speed of sound at pressure altitude, kt */
static double sound_kt(double alt_ft)
{
    double tk = alt_ft < 36089 ? 288.15 - 1.98e-3 * alt_ft : 216.65;
    return 38.967854 * sqrt(tk);
}

static bool is60(const uint8_t *m, const struct ac *a, double t, bool ref, struct b60 *o)
{
    if (!field_ok(m, 33, 44) || !field_ok(m, 45, 55) || !field_ok(m, 56, 66) || !field_ok(m, 67, 77) || !field_ok(m, 78, 88))
        return false;
    o->hdg_ok = getbits(m, 33, 33);  o->hdg = fmod(sgn(m, 34, 35, 44) * 90.0 / 512 + 360, 360);
    o->ias_ok = getbits(m, 45, 45);  o->ias = getbits(m, 46, 55);
    o->mach_ok = getbits(m, 56, 56); o->mach = getbits(m, 57, 66) * 2.048 / 512;
    o->baro_ok = getbits(m, 67, 67); o->baro = sgn(m, 68, 69, 77) * 32;
    o->ins_ok = getbits(m, 78, 78);  o->ins = sgn(m, 79, 80, 88) * 32;
    if (!o->hdg_ok || (!o->ias_ok && !o->mach_ok)) return false;
    if (o->ias_ok && (o->ias < 40 || o->ias > 500)) return false;
    if (o->mach_ok && (o->mach < 0.1 || o->mach > 0.95)) return false;
    if (o->baro_ok && fabs(o->baro) > 6000) return false;
    if (o->ins_ok && fabs(o->ins) > 6000) return false;
    if (o->baro_ok && o->ins_ok && fabs(o->baro - o->ins) > 2000) return false;
    bool alt = t - a->alt_t <= ALT_AGE;
    if (!ref)
        return o->ias_ok && o->mach_ok && (o->baro_ok || o->ins_ok) &&
               (!alt || fabs(cas_to_mach(o->ias, a->alt) - o->mach) <= 0.04);
    if (angdiff(o->hdg, a->trk) > 40) return false;             /* wind correction angle */
    if (a->vr_ok && o->ins_ok && fabs(o->ins - a->vr) > 1500) return false;
    if (a->vr_ok && !o->ins_ok && o->baro_ok && fabs(o->baro - a->vr) > 2000) return false;
    if (o->ias_ok && o->mach_ok && alt && fabs(cas_to_mach(o->ias, a->alt) - o->mach) > 0.04) return false;
    double m_ = o->mach_ok ? o->mach : alt ? cas_to_mach(o->ias, a->alt) : NAN;
    if (alt && !isnan(m_) && fabs(m_ * sound_kt(a->alt) - a->gs) > 200) return false;   /* TAS vs ground speed: wind */
    return true;
}

/* BDS 4,0 barometric setting, for the QNH: reserved bits zero, setting not standard */
static void bds40(const uint8_t *m, const struct ac *a, double t)
{
    if (!isnan(qnh_cmd) || getbits(m, 72, 79) || getbits(m, 84, 85)) return;
    if (!field_ok(m, 33, 45) || !field_ok(m, 46, 58) || !field_ok(m, 59, 71) || !field_ok(m, 80, 83) || !field_ok(m, 86, 88))
        return;
    if (!getbits(m, 59, 59) || t - a->alt_t > ALT_AGE || a->alt > QNH_MAX_FT) return;
    double q = getbits(m, 60, 71) * 0.1 + 800;
    if (q < 950 || q > 1060 || fabs(q - 1013.2) < 0.4) return;
    if (getbits(m, 33, 33) && (getbits(m, 34, 45) * 16 > 45000)) return;
    double before = qnh_now(t);
    qset[qset_n % QNH_N] = q;
    qset_t[qset_n % QNH_N] = t;
    qset_n++;
    double after = qnh_now(t);
    if (isnan(before) || fabs(after - before) >= 0.5)
        fprintf(stderr, "# QNH %.1f hPa (from BDS 4,0 of %06x %s, %d ft)\n", after, a->addr, a->cs, a->alt);
}

static void print_head(const struct ac *a, double t, double rssi, char *buf, size_t n)
{
    time_t s = (time_t)t;
    struct tm tm;
    localtime_r(&s, &tm);
    int o = snprintf(buf, n, "%02d:%02d:%02d.%03d  %06x  %-8s  ", tm.tm_hour, tm.tm_min, tm.tm_sec,
                     (int)((t - floor(t)) * 1000), a->addr, a->cs[0] ? a->cs : "-");
    double q = qnh_now(t);
    if (t - a->alt_t > ALT_AGE) o += snprintf(buf + o, n - o, "%14s", "alt ?");
    else if (isnan(q)) o += snprintf(buf + o, n - o, "%8d ft STD", a->alt);
    else o += snprintf(buf + o, n - o, "%9.0f ft   ", amsl_ft(a->alt, q));
    if (isnan(rssi)) snprintf(buf + o, n - o, "  %9s", "-");
    else snprintf(buf + o, n - o, "  %5.1f dBFS", rssi);
}

/* -A: the aircraft's hex address or callsign is in the list */
static bool selected(const struct ac *a)
{
    if (!nsel) return true;
    char hex[8];
    snprintf(hex, sizeof(hex), "%06x", a->addr);
    for (int k = 0; k < nsel; k++)
        if (!strcasecmp(sel[k], hex) || (a->cs[0] && !strcasecmp(sel[k], a->cs))) return true;
    return false;
}

static void commb(const uint8_t *m, int df, double t, double rssi)
{
    n_commb++;
    uint32_t addr = crc(m, 112);
    struct ac *a = ac_get(addr, false);
    bool given = hex_given(addr);
    if ((!a || t - a->seen > KNOWN_AGE) && !given) { n_unknown++; return; }
    if (!a && !(a = ac_get(addr, true))) return;
    if (df == 20) set_alt(a, ac13_ft(getbits(m, 20, 32)), t);
    bds40(m, a, t);
    struct b50 f5;
    struct b60 f6;
    bool ref = t - a->vel_t <= VEL_AGE;
    if (!ref && !given) { n_none++; return; }
    bool r5 = is50(m, a, ref, &f5), r6 = is60(m, a, t, ref, &f6);
    if (!r5 && !r6) { n_none++; return; }
    /* altitude filter: above it, or not known, not shown */
    if (max_alt != INT32_MAX && (t - a->alt_t > ALT_AGE || a->alt > max_alt)) { n_high++; return; }
    if (!selected(a)) { n_other++; return; }
    /* fits both registers: both decodes, each marked '?' */
    const char *mark = r5 && r6 ? "?" : " ";
    if (r5 && r6) n_both++;
    else if (r5) n_50++;
    else n_60++;
    char head[96], body[256];
    print_head(a, t, rssi, head, sizeof(head));
    if (r5) {
        int o = snprintf(body, sizeof(body), "BDS5,0%s", mark);
        if (f5.roll_ok) o += snprintf(body + o, sizeof(body) - o, " roll %+5.1f°", f5.roll);
        o += snprintf(body + o, sizeof(body) - o, "  track %5.1f°  GS %3.0f kt", f5.trk, f5.gs);
        if (f5.tas_ok) o += snprintf(body + o, sizeof(body) - o, "  TAS %3.0f kt", f5.tas);
        if (f5.rate_ok) snprintf(body + o, sizeof(body) - o, "  turn %+.2f°/s", f5.rate);
        printf("%s  %s%s\n", head, body, ref ? "" : "  (no ADS-B check)");
    }
    if (r6) {
        int o = snprintf(body, sizeof(body), "BDS6,0%s hdg %5.1f°", mark, f6.hdg);
        if (f6.ias_ok) o += snprintf(body + o, sizeof(body) - o, "  IAS %3.0f kt", f6.ias);
        if (f6.mach_ok) o += snprintf(body + o, sizeof(body) - o, "  M%.3f", f6.mach);
        if (f6.baro_ok) o += snprintf(body + o, sizeof(body) - o, "  baro %+5.0f ft/min", f6.baro);
        if (f6.ins_ok) snprintf(body + o, sizeof(body) - o, "  inertial %+5.0f ft/min", f6.ins);
        printf("%s  %s%s\n", head, body, ref ? "" : "  (no ADS-B check)");
    }
    fflush(stdout);
}

/* rssi: received level from the Beast signal byte, dBFS (NAN when the receiver gives none) */
static void frame(const uint8_t *m, int len, double t, double rssi)
{
    int df = m[0] >> 3;
    if (len == 14 && (df == 17) && crc(m, 112) == 0) df17(m, t);
    else if (len == 14 && (df == 20 || df == 21)) commb(m, df, t, rssi);
    else if (len == 7 && (df == 0 || df == 4)) {
        uint32_t addr = crc(m, 56);
        struct ac *a = ac_get(addr, false);
        if (!a && hex_given(addr)) a = ac_get(addr, true);
        if (a && (t - a->seen <= KNOWN_AGE || hex_given(addr))) set_alt(a, ac13_ft(getbits(m, 20, 32)), t);
    }
}

/* ---------------- Beast stream ---------------- */

static int dial(const char *host, const char *port)
{
    struct addrinfo hints = { .ai_socktype = SOCK_STREAM }, *res;
    if (getaddrinfo(host, port, &hints, &res)) return -1;
    int fd = -1;
    for (struct addrinfo *r = res; r; r = r->ai_next) {
        fd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);
        if (fd < 0) continue;
        if (!connect(fd, r->ai_addr, r->ai_addrlen)) break;
        close(fd);
        fd = -1;
    }
    freeaddrinfo(res);
    return fd;
}

/* frame time: the receiver's 12 MHz timestamp, placed on the wall clock by the
 * smallest (wall - timestamp) seen, since frames reach us in batches */
static double frame_time(const uint8_t *ts)
{
    static double off = NAN, last;
    uint64_t v = 0;
    for (int k = 0; k < 6; k++) v = v << 8 | ts[k];
    double w = now_s(), r = v / 12e6;
    if (!v) return w;                                   /* no timestamp */
    if (isnan(off) || r < last - 1 || w - r < off) off = w - r;      /* first, receiver restarted, or less delay */
    last = r;
    return r + off;
}

/* Beast: 0x1a, type ('1' A/C 2 bytes, '2' short 7, '3' long 14), 6-byte timestamp,
 * 1 signal byte, the message; 0x1a inside is doubled */
static void beast(int fd)
{
    uint8_t buf[4096], f[32];
    int st = 0, need = 0, got = 0;
    bool esc = false;
    double last_stat = now_s();
    for (;;) {
        ssize_t n = read(fd, buf, sizeof(buf));
        if (n <= 0) return;
        for (ssize_t i = 0; i < n; i++) {
            uint8_t c = buf[i];
            if (st == 0) { if (c == 0x1a) st = 1; continue; }
            if (st == 1) {
                need = c == '1' ? 2 : c == '2' ? 7 : c == '3' ? 14 : 0;
                if (!need) { st = c == 0x1a ? 1 : 0; continue; }
                need += 7; got = 0; esc = false; st = 2;
                continue;
            }
            if (c == 0x1a && !esc) { esc = true; continue; }
            if (esc && c != 0x1a) { st = 1; i--; esc = false; continue; }   /* unescaped 0x1a: a new frame */
            esc = false;
            f[got++] = c;
            if (got == need) {
                /* signal byte: sqrt of the signal power (full scale 1) x 255, as dump1090 writes it */
                if (need > 9) frame(f + 7, need - 7, frame_time(f), f[6] ? 20 * log10(f[6] / 255.0) : NAN);
                st = 0;
            }
        }
        double t = now_s();
        if (verbose && t - last_stat > 60) {
            last_stat = t;
            fprintf(stderr, "# Comm-B %ld: unknown aircraft %ld, BDS5,0 %ld, BDS6,0 %ld, above --max-alt %ld, other aircraft (-A) %ld, ambiguous %ld, neither %ld\n",
                    n_commb, n_unknown, n_50, n_60, n_high, n_other, n_both, n_none);
        }
    }
}

static void usage(void)
{
    fprintf(stderr,
            "usage: bds_logger [HOST:PORT] [--qnh HPA] [--max-alt FT] [-A LIST] [-v]\n"
            "  HOST:PORT   Beast output (default localhost:30005)\n"
            "  --qnh HPA   QNH for AMSL altitudes; default: from BDS 4,0 settings of aircraft below %d ft,\n"
            "              or the pressure altitude (STD) until one is seen\n"
            "  --max-alt FT  hide aircraft above this pressure altitude (or with no altitude known)\n"
            "  -A LIST     only these aircraft: hex addresses or callsigns, separated by commas\n"
            "              (e.g. -A 48af10,LOT3TD; a callsign matches once its ADS-B identification is heard).\n"
            "              An aircraft given as hex is decoded also without ADS-B: then all fields of the\n"
            "              register must be present, and the line ends '(no ADS-B check)'\n"
            "  -v          statistics every minute on stderr\n"
            "Prints Comm-B replies of ADS-B aircraft that are BDS 5,0 or 6,0, one per line.\n", QNH_MAX_FT);
}

int main(int argc, char **argv)
{
    const char *hp = "localhost:30005";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--qnh") && i + 1 < argc) {
            qnh_cmd = atof(argv[++i]);
            if (qnh_cmd < 900 || qnh_cmd > 1100) { fprintf(stderr, "--qnh: %s is not a QNH in hPa\n", argv[i]); return 1; }
        } else if (!strcmp(argv[i], "--max-alt") && i + 1 < argc) {
            char *e;
            long v = strtol(argv[++i], &e, 10);
            if (*e || v < -1000 || v > 60000) { fprintf(stderr, "--max-alt: %s is not an altitude in ft\n", argv[i]); return 1; }
            max_alt = (int)v;
        } else if (!strcmp(argv[i], "-A") && i + 1 < argc) {
            char list[512];
            snprintf(list, sizeof(list), "%s", argv[++i]);
            for (char *sv, *tok = strtok_r(list, ", ", &sv); tok; tok = strtok_r(NULL, ", ", &sv)) {
                if (nsel == SEL_N || strlen(tok) > 8) { fprintf(stderr, "-A: too many entries or '%s' too long\n", tok); return 1; }
                snprintf(sel[nsel++], sizeof(sel[0]), "%s", tok);
            }
        } else if (!strcmp(argv[i], "-v")) verbose = true;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(); return 0; }
        else if (argv[i][0] != '-') hp = argv[i];
        else { usage(); return 1; }
    }
    char host[256], port[16] = "30005";
    snprintf(host, sizeof(host), "%s", hp);
    char *c = strrchr(host, ':');
    if (c) { snprintf(port, sizeof(port), "%s", c + 1); *c = 0; }
    for (;;) {
        int fd = dial(host, port);
        if (fd < 0) { fprintf(stderr, "# cannot connect to %s:%s, retrying in 5 s\n", host, port); sleep(5); continue; }
        fprintf(stderr, "# connected to %s:%s\n", host, port);
        beast(fd);
        close(fd);
        fprintf(stderr, "# connection to %s:%s lost, reconnecting\n", host, port);
        sleep(1);
    }
}
