/* Statistics, main-beam burst detection, radar attribution, filtering and
 * interleaving analysis.
 *
 * Stage 1 (all interrogations, delayed by LOOKAHEAD): main-beam burst
 * detection. The background level is the median, over the last 20 s, of the
 * strongest radar interrogation in each 10 ms bin (empty bins count as the
 * detection threshold) - i.e. the sidelobe level, or the floor when only the
 * main beam is audible. A burst is a run of interrogations >= member_db above
 * background, spaced < burst_gap, containing one >= core_db above background.
 * While a burst is open everything is held in an output queue so the burst's
 * radar identity (majority of II/SI codes it carries) is known before stage 2.
 *
 * Stage 2 (in time order): radar attribution, filters, statistics, logging.
 * Attribution: code carried by the message (UF11 II/SI, roll-call SIS/IIS),
 * else the code of the burst it belongs to, else for an all-call the code of
 * a UF11 sent < 300 us before it. Per-burst statistics are built from the
 * interrogations that pass the filters. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "analysis.h"

#define LOOKAHEAD   0.3
#define BIN_S       0.01
#define NBINS       2000
#define PENDING     65536
#define OUTQ        65536
#define BURST_CODE_MIN 2           /* queries with the code needed to name a burst (temporary) */
#define PAIR_S      0.0003

struct analysis_state AS;

static struct interrog pend[PENDING];
static uint64_t ph, pt;                 /* head (write), tail (oldest) */
static struct interrog outq[OUTQ];
static uint64_t oh, ot;
static float bins[NBINS];
static long cur_bin = -1;
static int nbins_filled, bins_since_bg;
static int next_burst_id;
static bool (*rx_ready)(double t);
static void (*rx_match)(struct interrog *x);

void analysis_set_rx(bool (*ready)(double t), void (*match)(struct interrog *x))
{
    rx_ready = ready;
    rx_match = match;
}
static void (*sinkf)(const struct interrog *it);
static double last11_t = -1;
static int last11_code = -1;
static unsigned target_ver = ~0u;

static const char *group_labels[G_COUNT] = {
    "AA/CC", "A/C", "ACSall", "S", "alt", "sqk", "acid", "selAlt", "trk", "hdgSpd", "oBDS",
    "msg", "other"
};

int tok_group(uint16_t tok)
{
    int k = tok & 31, x = tok >> 5;
    if (x) {
        switch (x - 1) {
        case 0x20: return G_CAL;
        case 0x40: return G_SETALT;
        case 0x50: return G_TRK;
        case 0x60: return G_HDGSPD;
        default:   return G_OBDS;
        }
    }
    switch (k) {
    case K_MODEA: case K_MODEC: case K_MODE1: case K_MODE2: case K_MODEB: case K_MODED: return G_AC;
    case K_ACONLY_A: case K_ACONLY_C: return G_ACONLY;
    case K_ACS_A: case K_ACS_C: return G_ACS;
    case K_UF11: return G_SALL;
    case K_UF4: return G_SALT;
    case K_UF5: return G_SSQK;
    case K_UF20: case K_UF21: return G_MSG;
    default: return G_SOTHER;
    }
}

const char *group_label(int g)
{
    static const char *const reg[G_COUNT] = { [G_SALT] = "UF4", [G_SSQK] = "UF5", [G_CAL] = "b20", [G_SETALT] = "b40",
                                              [G_TRK] = "b50", [G_HDGSPD] = "b60" };
    return !bds_human && reg[g] ? reg[g] : group_labels[g];
}

static int cmpd(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
static int cmpfl(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

double median_d(const double *v, int n)
{
    if (n <= 0) return NAN;
    double t[64];
    if (n > 64) n = 64;
    memcpy(t, v, n * sizeof(double));
    qsort(t, n, sizeof(double), cmpd);
    return (n & 1) ? t[n / 2] : 0.5 * (t[n / 2 - 1] + t[n / 2]);
}

void analysis_init(const struct analysis_cfg *cfg)
{
    memset(&AS, 0, sizeof(AS));
    AS.cfg = *cfg;
    AS.floor_db = -20;
    AS.bg_db = -20;
    AS.prev_tok = -1;
    AS.last_burst_end = -1;
    AS.gap_end = -1e9;
    for (int k = 0; k < AC_MAX; k++) AS.tgt[k] = k < cfg->nac ? cfg->ac[k].addr : -1;
    target_ver = ~0u;
}

void analysis_set_floor(double thr_dbfs) { AS.floor_db = thr_dbfs; }
void analysis_set_sink(void (*sink)(const struct interrog *it)) { sinkf = sink; }

/* ---------------- background level ---------------- */

static void bins_advance(long b)
{
    if (cur_bin < 0) cur_bin = b - 1;
    while (cur_bin < b) {
        cur_bin++;
        bins[cur_bin % NBINS] = (float)AS.floor_db;
        if (nbins_filled < NBINS) nbins_filled++;
        if (++bins_since_bg >= 50) {
            static float tmp[NBINS];
            for (int k = 0; k < nbins_filled; k++) tmp[k] = bins[(cur_bin - k) % NBINS];
            qsort(tmp, nbins_filled, sizeof(float), cmpfl);
            AS.bg_db = tmp[nbins_filled / 2];
            bins_since_bg = 0;
        }
    }
}

void analysis_add(struct interrog *it)
{
    if (!it->tcas) {
        long b = (long)(it->t / BIN_S);
        if (b >= cur_bin) {
            bins_advance(b);
            float *v = &bins[b % NBINS];
            if (it->amp_db > *v) *v = it->amp_db;
        }
    }
    if (ph - pt >= PENDING) analysis_tick(1e300);
    it->radar = -1;
    pend[ph++ % PENDING] = *it;
}

/* ---------------- radars ---------------- */

static struct radar *radar_get(int code)
{
    for (int r = 0; r < AS.nradars; r++)
        if (AS.radars[r].code == code) return &AS.radars[r];
    if (AS.nradars == MAX_RADARS) return &AS.radars[MAX_RADARS - 1];
    struct radar *r = &AS.radars[AS.nradars++];
    memset(r, 0, sizeof(*r));
    r->code = code;
    radar_code_str(code, r->key, sizeof(r->key));
    for (int g = 0; g < G_COUNT; g++) r->cnt_min[g] = 1 << 30;
    return r;
}

double analysis_rotation(void)
{
    struct radar *best = NULL;
    for (int r = 0; r < AS.nradars; r++) {
        struct radar *X = &AS.radars[r];
        if (AS.cfg.radar_sel >= 0 && X->code != AS.cfg.radar_sel) continue;
        if (X->code < 0 && AS.cfg.radar_sel < 0 && best) continue;
        if (!best || X->nbursts > best->nbursts || (best->code < 0 && X->code >= 0)) best = X;
    }
    if (!best || best->np == 0) return 0;
    return median_d(best->periods, best->np < 64 ? best->np : 64);
}

double analysis_rotation_of(const char *key)
{
    for (int r = 0; r < AS.nradars; r++) {
        struct radar *X = &AS.radars[r];
        if (strcmp(X->key, key) || X->np == 0) continue;
        return median_d(X->periods, X->np < 64 ? X->np : 64);
    }
    return 0;
}

double analysis_last_burst_end(void) { return AS.last_burst_end; }

int analysis_main_code(void)
{
    if (AS.cfg.radar_sel >= 0) return AS.cfg.radar_sel;
    struct radar *best = NULL;
    for (int r = 0; r < AS.nradars; r++) {
        struct radar *X = &AS.radars[r];
        if (!best || (best->code < 0 && X->code >= 0) ||
            ((X->code >= 0) == (best->code >= 0) && X->nbursts > best->nbursts)) best = X;
    }
    return best ? best->code : -2;
}

void analysis_set_gain(float g) { AS.gain_db = g; }

void analysis_gain_shift(float delta)
{
    for (int k = 0; k < NBINS; k++) bins[k] += delta;
    AS.bg_db += delta;
}

/* ---------------- filters ---------------- */

static bool aircraft_filter_on(void) { return AS.cfg.nac > 0; }

static void resolve_target(void)
{
    if (!AS.cfg.nac || cs_version() == target_ver) return;
    target_ver = cs_version();
    for (int k = 0; k < AS.cfg.nac; k++)
        AS.tgt[k] = AS.cfg.ac[k].addr >= 0 ? AS.cfg.ac[k].addr : cs_find(AS.cfg.ac[k].cs);
}

int ac_list_parse(const char *s, struct ac_sel *out, char *bad, size_t badlen)
{
    char buf[256];
    snprintf(buf, sizeof(buf), "%s", s);
    int n = 0;
    for (char *t = strtok(buf, " ,;\t"); t; t = strtok(NULL, " ,;\t")) {
        if (n == AC_MAX) { snprintf(bad, badlen, "more than %d aircraft", AC_MAX); return -1; }
        struct ac_sel *e = &out[n];
        e->addr = -1;
        e->cs[0] = 0;
        const char *cs = NULL;
        if (!strncasecmp(t, "cs:", 3)) cs = t + 3;
        else {
            const char *h = !strncasecmp(t, "0x", 2) ? t + 2 : t;
            bool hex = strlen(h) == 6;
            for (int k = 0; hex && k < 6; k++) hex = isxdigit((unsigned char)h[k]);
            if (hex) e->addr = strtol(h, NULL, 16);
            else cs = t;
        }
        if (cs) {
            if (!*cs || strlen(cs) > 8) { snprintf(bad, badlen, "'%s'", t); return -1; }
            for (int k = 0; cs[k]; k++) e->cs[k] = toupper((unsigned char)cs[k]);
            e->cs[strlen(cs)] = 0;
        }
        n++;
    }
    return n;
}

void ac_list_str(const struct analysis_cfg *c, bool resolved, char *buf, size_t n)
{
    size_t o = 0;
    buf[0] = 0;
    for (int k = 0; k < c->nac && o + 1 < n; k++) {
        const struct ac_sel *e = &c->ac[k];
        const char *sep = k ? "," : "";
        if (e->addr >= 0) o += snprintf(buf + o, n - o, "%s%06llx", sep, (unsigned long long)e->addr);
        else if (!resolved) o += snprintf(buf + o, n - o, "%s%s", sep, e->cs);
        else {
            int64_t a = cs_find(e->cs);
            if (a >= 0) o += snprintf(buf + o, n - o, "%s%s=%06llx", sep, e->cs, (unsigned long long)a);
            else o += snprintf(buf + o, n - o, "%s%s(?)", sep, e->cs);
        }
    }
}

static bool directed(int k)
{
    return k == K_UF0 || k == K_UF4 || k == K_UF5 || k == K_UF16 || k == K_UF20 || k == K_UF21 ||
           k == K_UF24 || k == K_UF19;
}

static bool to_target(const struct interrog *x)
{
    if (!directed(x->kind)) return false;
    for (int k = 0; k < AS.cfg.nac; k++)
        if (AS.tgt[k] >= 0 && x->addr == (uint32_t)AS.tgt[k]) return true;
    return false;
}

/* II0 query: UF11 with II 0, or a roll-call with IIS 0 that is not part of
 * another radar's burst (SI radars send IIS 0 too) */
static bool is_ii0(const struct interrog *x)
{
    if (x->kind == K_UF11) return modes_own_code(x) == 0;
    return modes_iis0(x) && x->radar <= 0;
}

static bool passes(const struct interrog *x)
{
    if (AS.cfg.radar_sel >= 0 && x->radar != AS.cfg.radar_sel) { AS.filt_radar++; return false; }
    if (AS.cfg.modes_only && !kind_is_modes(x->kind)) { AS.filt_modes++; return false; }
    if (AS.cfg.ac_only && x->kind != K_MODEA && x->kind != K_MODEC) { AS.filt_ac++; return false; }
    if (AS.cfg.no_ac && (x->kind == K_MODEA || x->kind == K_MODEC)) { AS.filt_noac++; return false; }
    if (AS.cfg.ii0_only && !is_ii0(x)) { AS.filt_ii0++; return false; }
    if (AS.cfg.no_ii0 && is_ii0(x)) { AS.filt_ii0++; return false; }
    if (AS.cfg.rollcall_only && !directed(x->kind)) { AS.filt_rollcall++; return false; }
    if (aircraft_filter_on() && !to_target(x)) { AS.filt_aircraft++; return false; }
    return true;
}

/* ---------------- per-interrogation statistics (stage 2) ---------------- */

static struct addr_entry *addr_slot(uint32_t a, bool tcas, bool create)
{
    uint32_t h = (a * 2654435761u) >> 19;
    for (int k = 0; k < ADDR_HASH; k++) {
        struct addr_entry *e = &AS.addrs[(h + k) % ADDR_HASH];
        if (e->n == 0) {
            if (!create) return NULL;
            e->addr = a; e->tcas = tcas;
            if (tcas) AS.naddr_tcas++; else AS.naddr_radar++;
            return e;
        }
        if (e->addr == a && e->tcas == tcas) return e;
    }
    return NULL;
}

static void addr_note(uint32_t a, double t, bool tcas)
{
    struct addr_entry *e = addr_slot(a, tcas, true);
    if (e) { e->n++; e->last = t; }
}

bool analysis_addr_known(uint32_t a)
{
    return addr_slot(a, false, false) != NULL || addr_slot(a, true, false) != NULL;
}

static void attribute(struct interrog *x, int burst_code)
{
    int own = modes_own_code(x);
    if (own >= 0) x->radar = own;
    else if (x->burst >= 0 && burst_code >= 0) x->radar = burst_code;
    else if (x->kind != K_UF11 && kind_is_allcall(x->kind) && last11_t >= 0 && x->t - last11_t < PAIR_S)
        x->radar = last11_code;
    else x->radar = -1;
    if (x->kind == K_UF11 && own >= 0) { last11_t = x->t; last11_code = own; }
}

static struct trans *trans_get(uint16_t a, uint16_t b)
{
    for (int i = 0; i < AS.ntr; i++)
        if (AS.tr[i].a == a && AS.tr[i].b == b) return &AS.tr[i];
    if (AS.ntr == TR_N) return NULL;
    struct trans *t = &AS.tr[AS.ntr++];
    t->a = a; t->b = b; t->min = 1e9;
    return t;
}

static int cmp_trans(const void *x, const void *y)
{
    long a = ((const struct trans *)x)->n, b = ((const struct trans *)y)->n;
    return (b > a) - (b < a);
}

int analysis_top_transitions(struct trans *out, int max)
{
    static struct trans tmp[TR_N];
    memcpy(tmp, AS.tr, sizeof(struct trans) * AS.ntr);
    qsort(tmp, AS.ntr, sizeof(struct trans), cmp_trans);
    int n = AS.ntr < max ? AS.ntr : max;
    memcpy(out, tmp, sizeof(struct trans) * n);
    return n;
}

static void graph_note(const struct interrog *y);

/* returns true when the interrogation passed the filters */
static bool stage2(struct interrog *x, int burst_code)
{
    AS.n_all++;
    resolve_target();
    attribute(x, burst_code);
    graph_note(x);

    if (x->tcas) {
        AS.tcas_total++;
        if (x->kind == K_UF0) AS.tcas_uf0++;
        else if (x->kind == K_UF16) AS.tcas_uf16++;
        else AS.tcas_ws++;
        if ((x->kind == K_UF0 || x->kind == K_UF16) && x->addr != 0xFFFFFF) addr_note(x->addr, x->t, true);
        if (sinkf && (!aircraft_filter_on() || to_target(x))) sinkf(x);
        return false;
    }
    if (!passes(x)) return false;
    AS.n_pass++;

    int k = x->kind;
    AS.total[k]++;
    if (x->burst >= 0) AS.inbeam[k]++;
    if (x->sidelobe) AS.side[k]++;

    if (x->burst >= 0) {
        int tok = seq_token(x);
        if (AS.prev_tok >= 0 && AS.prev_burst == x->burst) {
            struct trans *t = trans_get((uint16_t)AS.prev_tok, (uint16_t)tok);
            if (t) {
                double dt = x->t - AS.prev_t;
                t->n++; t->sum += dt;
                if (dt < t->min) t->min = dt;
                if (dt > t->max) t->max = dt;
            }
        }
        AS.prev_tok = tok;
        AS.prev_t = x->t;
        AS.prev_burst = x->burst;
    }

    if (k == K_UF11) {
        int si, c = modes_iisi(x->msg, &si);
        if (c >= 0 && c < 64) AS.uf11_code[si][c][getbits(x->msg, 6, 9)]++;
    }
    if (k == K_UF4 || k == K_UF5 || k == K_UF20 || k == K_UF21) {
        int b1, b2;
        if (modes_bds(x->msg, &b1, &b2)) {
            int b = b1 * 16 + b2;
            AS.bds[b]++;
            if (x->burst >= 0) AS.bds_inbeam[b]++;
            if (x->sidelobe) AS.bds_side[b]++;
            AS.commb[k]++;
            if (x->burst >= 0) AS.commb_inbeam[k]++;
            if (x->sidelobe) AS.commb_side[k]++;
        }
        AS.pc[getbits(x->msg, 6, 8)]++;
        AS.di[getbits(x->msg, 14, 16)]++;
    }
    if (k == K_UF4 || k == K_UF5 || k == K_UF20 || k == K_UF21 || k == K_UF24)
        addr_note(x->addr, x->t, false);

    AS.recent[AS.nrecent % RECENT_N] = *x;
    AS.nrecent++;
    if (sinkf) sinkf(x);
    return true;
}

/* ---------------- bursts (stage 1) ---------------- */

static bool is_member(const struct interrog *x)
{
    return !x->tcas && !x->sidelobe && x->amp_db >= AS.bg_db + AS.cfg.member_db;
}
static bool is_core(const struct interrog *x)
{
    return is_member(x) && x->amp_db >= AS.bg_db + AS.cfg.core_db;
}

static int rec_add(struct burst *b, const struct interrog *x)
{
    if (b->n == 0) { b->t0 = x->t; b->peak = -99; }
    b->n++;
    b->cnt[x->kind]++;
    b->t1 = x->t;
    if (x->amp_db > b->peak) { b->peak = x->amp_db; b->tpeak = x->t; }
    uint16_t tok = seq_token(x);
    b->gcnt[tok_group(tok)]++;
    int k = NO_TARGET;
    if (directed(x->kind)) {
        for (k = 0; k < b->naddr; k++) if (b->addrs[k] == x->addr) break;
        if (k == b->naddr) {
            if (b->naddr < 64) b->addrs[b->naddr++] = x->addr;
            else k = NO_TARGET;
        }
    }
    if (b->nks < SEQ_MAX) { b->ks[b->nks] = tok; b->ka[b->nks] = (uint8_t)k; b->kans[b->nks] = x->answered; b->nks++; }
    return k;
}

/* non-member queries shown in the burst graph: same type/aircraft filters,
 * radar filter only when the query could be attributed */
static bool graph_extra_ok(const struct interrog *y)
{
    if (y->tcas) return false;
    if (AS.cfg.radar_sel >= 0 && y->radar >= 0 && y->radar != AS.cfg.radar_sel) return false;
    if (AS.cfg.modes_only && !kind_is_modes(y->kind)) return false;
    if (AS.cfg.ac_only && y->kind != K_MODEA && y->kind != K_MODEC) return false;
    if (AS.cfg.no_ac && (y->kind == K_MODEA || y->kind == K_MODEC)) return false;
    if (AS.cfg.ii0_only && !is_ii0(y)) return false;
    if (AS.cfg.no_ii0 && is_ii0(y)) return false;
    if (AS.cfg.rollcall_only && !directed(y->kind)) return false;
    if (aircraft_filter_on() && !to_target(y)) return false;
    return true;
}

/* II/SI code as sent; GP_IIS0 for a roll-call with IIS 0, which SI radars
 * send too (see modes_iis0) */
static int radio_code(const struct interrog *y)
{
    int c = modes_own_code(y);
    return c >= 0 ? c : modes_iis0(y) ? GP_IIS0 : -1;
}

static void graph_put(struct burst *b, double t, float amp, float p2, float p5, float sq, float noise, uint32_t addr, uint16_t tok,
                      int ka, int cls, bool ans, uint64_t rep_idx, int rep_n, int code)
{
    if (b->ngp >= SEQ_MAX) return;
    struct gpoint *g = &b->gp[b->ngp++];
    g->t = t;
    g->addr = addr;
    g->amp = amp;
    g->p2 = p2;
    g->p5 = p5;
    g->sq = sq;
    g->noise = noise;
    g->tok = tok;
    g->ka = (uint8_t)ka;
    g->cls = (uint8_t)cls;
    g->ans = ans;
    g->rep_idx = rep_idx;
    g->rep_n = (uint8_t)rep_n;
    g->code = (int8_t)code;
    if (cls == GP_SIDE) b->nside++;
    if (cls == GP_WEAK) b->nweak++;
}

/* absolute SLS pulse levels: P2 only for A/C (in Mode S it is the preamble) */
static float abs_p2(const struct interrog *y)
{
    return !kind_is_modes(y->kind) && !isnan(y->p2_db) ? y->amp_db + y->p2_db : NAN;
}
static float abs_p5(const struct interrog *y)
{
    return kind_is_modes(y->kind) && !isnan(y->p5_db) ? y->amp_db + y->p5_db : NAN;
}

static uint32_t gaddr(const struct interrog *y) { return directed(y->kind) ? y->addr : NO_ADDR; }

static void graph_add(struct burst *b, const struct interrog *y, int cls, int ka)
{
    graph_put(b, y->t, y->amp_db, abs_p2(y), abs_p5(y), y->sls_floor_db, y->noise_db, gaddr(y), seq_token(y), ka, cls, y->answered, y->rep_idx, y->rep_n, radio_code(y));
}

/* Non-beam queries for the graph. Sidelobes are kept from 60 % of a rotation
 * before to 60 % of a rotation after a burst, other weak queries only within
 * one dwell. A ring holds the recent past; open post-windows collect the
 * queries that arrive after a burst was committed. */
#define GRING 16384
#define NPOST 8
static struct {
    double t; float amp, p2, p5, sq; float noise; uint32_t addr; uint16_t tok; uint8_t cls; bool ans; uint8_t rep_n; int8_t code; uint64_t rep_idx;
} gring[GRING];
static uint64_t gring_n;
static struct { int idx; double until_side, until_weak; } posts[NPOST];
static int nposts;

static void graph_note(const struct interrog *y)
{
    if (y->burst >= 0 || !graph_extra_ok(y)) return;
    int cls = y->sidelobe ? GP_SIDE : GP_WEAK;
    uint16_t tok = seq_token(y);
    for (int p = 0; p < nposts; p++) {
        double until = cls == GP_SIDE ? posts[p].until_side : posts[p].until_weak;
        if (y->t <= until && AS.nhist - posts[p].idx <= BURST_HIST)
            graph_put(&AS.hist[posts[p].idx % BURST_HIST], y->t, y->amp_db, abs_p2(y), abs_p5(y), y->sls_floor_db, y->noise_db, gaddr(y), tok,
                      NO_TARGET, cls, y->answered, y->rep_idx, y->rep_n, radio_code(y));
    }
    int i = gring_n++ % GRING;
    gring[i].t = y->t; gring[i].amp = y->amp_db; gring[i].tok = tok; gring[i].cls = (uint8_t)cls;
    gring[i].p2 = abs_p2(y); gring[i].p5 = abs_p5(y); gring[i].sq = y->sls_floor_db; gring[i].noise = y->noise_db; gring[i].addr = gaddr(y); gring[i].ans = y->answered;
    gring[i].rep_idx = y->rep_idx; gring[i].rep_n = y->rep_n; gring[i].code = (int8_t)radio_code(y);
}

static int cmp_gp(const void *x, const void *y)
{
    double a = ((const struct gpoint *)x)->t, b = ((const struct gpoint *)y)->t;
    return (a > b) - (a < b);
}

/* sidelobe shoulder shown on each side of a burst: 60 % of the rotation */
#define SIDE_FRAC 0.6
static double side_window(int code)
{
    for (int r = 0; r < AS.nradars; r++)
        if (AS.radars[r].code == code && AS.radars[r].np)
            return SIDE_FRAC * median_d(AS.radars[r].periods, AS.radars[r].np < 64 ? AS.radars[r].np : 64);
    double rot = analysis_rotation();
    return rot > 0 ? SIDE_FRAC * rot : 2.4;
}

/* before commit: add buffered non-beam queries preceding / inside the burst */
static void graph_context(struct burst *b, double t0, double t1, double margin, double half)
{
    uint64_t first = gring_n > GRING ? gring_n - GRING : 0;
    for (uint64_t k = first; k < gring_n; k++) {
        int i = k % GRING;
        double lo = gring[i].cls == GP_SIDE ? t0 - half : t0 - margin;
        double hi = gring[i].cls == GP_SIDE ? t1 + half : t1 + margin;
        if (gring[i].t < lo || gring[i].t > hi) continue;
        graph_put(b, gring[i].t, gring[i].amp, gring[i].p2, gring[i].p5, gring[i].sq, gring[i].noise, gring[i].addr, gring[i].tok, NO_TARGET,
                  gring[i].cls, gring[i].ans, gring[i].rep_idx, gring[i].rep_n, gring[i].code);
    }
    qsort(b->gp, b->ngp, sizeof(b->gp[0]), cmp_gp);
}

const struct krow krows[] = {
    { NULL, K_MODE1, 0, false }, { NULL, K_MODE2, 0, false }, { NULL, K_MODEA, 0, false },
    { NULL, K_MODEB, 0, false }, { NULL, K_MODEC, 0, false }, { NULL, K_MODED, 0, false },
    { NULL, K_ACONLY_A, 0, true }, { NULL, K_ACONLY_C, 0, true },
    { NULL, K_ACS_A, 0, false }, { NULL, K_ACS_C, 0, false }, { NULL, K_INTERMODE_OTHER, 0, false },
    { NULL, K_UF4, 1, true }, { NULL, K_UF5, 1, true }, { NULL, K_UF11, 0, true },
    { "UF4  Comm-B req -> DF20", K_UF4, 2, true }, { "UF5  Comm-B req -> DF21", K_UF5, 2, true },
    { "UF other", K_UF19, 3, true }, { NULL, K_S_NOSYNC, 0, true }, { NULL, K_S_GARBLE, 0, true },
};
const int nkrows = sizeof(krows) / sizeof(krows[0]);

/* one counter of a table row, from per-kind counters (total, in-beam, sidelobe or a snapshot) */
long krow_count(const struct krow *r, const long *by_kind, const long *commb_by_kind)
{
    switch (r->part) {
    case 1: return by_kind[r->k] - commb_by_kind[r->k];        /* surveillance only: DF4/DF5 reply */
    case 2: return commb_by_kind[r->k];                         /* Comm-B request */
    case 3: return by_kind[K_UF19] + by_kind[K_UF20] + by_kind[K_UF21] + by_kind[K_UF24] + by_kind[K_UF_OTHER];
    default: return by_kind[r->k];
    }
}

long krow_all(void)
{
    long n = 0;
    for (int i = 0; i < nkrows; i++) n += krow_count(&krows[i], AS.total, AS.commb);
    return n;
}

bool krow_shown(const struct krow *r, long tot, long all)
{
    return r->always || tot >= KROW_MIN_SHARE * all;
}

const struct burst *analysis_graph_burst(void)
{
    if (!AS.nhist) return NULL;
    int code = analysis_main_code();
    int lim = AS.nhist < BURST_HIST ? AS.nhist : BURST_HIST;
    for (int i = 0; i < lim; i++) {
        const struct burst *b = &AS.hist[(AS.nhist - 1 - i) % BURST_HIST];
        if (b->code == code) return b;
    }
    return &AS.hist[(AS.nhist - 1) % BURST_HIST];
}

static void commit(struct burst *b)
{
    struct radar *r = radar_get(b->code);
    snprintf(b->key, sizeof(b->key), "%s", r->key);
    b->period = 0;
    if (b->partial) {
        r->npartial++;
        r->gap_since = true;
    } else {
        if (r->nbursts > 0) {
            double p = b->tpeak - r->last_peak;
            int np = r->np < 64 ? r->np : 64;
            double med = np ? median_d(r->periods, np) : 0;
            bool ok = p > 0.5 && p < 60;
            if (r->gap_since) ok = ok && np >= 3 && fabs(p - med) < 0.25 * med;  /* no burst missed */
            if (ok) {
                b->period = p;
                r->periods[r->np % 64] = p;
                r->np++;
            }
        }
        r->gap_since = false;
        r->last_peak = b->tpeak;
        r->nbursts++;
        r->dwell_sum += b->t1 - b->t0;
        r->total_sum += b->n;
        const int *gc = b->gcnt;
        for (int g = 0; g < G_COUNT; g++) {
            r->cnt_sum[g] += gc[g];
            if (gc[g] < r->cnt_min[g]) r->cnt_min[g] = gc[g];
            if (gc[g] > r->cnt_max[g]) r->cnt_max[g] = gc[g];
        }
    }
    b->gain = AS.gain_db;
    AS.hist[AS.nhist % BURST_HIST] = *b;
    AS.nhist++;
}

/* ---------------- burst rendering helpers ---------------- */

static const struct burst *rank_b;
static int cmp_target(const void *x, const void *y)
{
    char a[12], b[12];
    snprintf(a, sizeof(a), "%s", ac_name(rank_b->addrs[*(const int *)x]));
    snprintf(b, sizeof(b), "%s", ac_name(rank_b->addrs[*(const int *)y]));
    return strcmp(a, b);
}

void burst_ranks(const struct burst *b, int *rank)
{
    int idx[64];
    for (int k = 0; k < b->naddr; k++) idx[k] = k;
    rank_b = b;
    qsort(idx, b->naddr, sizeof(int), cmp_target);
    for (int r = 0; r < b->naddr; r++) rank[idx[r]] = r + 1;
}

int burst_item(const struct burst *b, int s, const int *rank, char *out, size_t n)
{
    int e = s;
    while (e < b->nks && b->ks[e] == b->ks[s] && b->ka[e] == b->ka[s] && b->kans[e] == b->kans[s]) e++;
    char ts[24], tg[8] = "";
    seq_token_str(b->ks[s], ts, sizeof(ts));
    if (b->ka[s] != NO_TARGET) snprintf(tg, sizeof(tg), "/%d", rank[b->ka[s]]);
    if (b->kans[s]) strcat(tg, "*");                /* answered on 1090 */
    if (e - s > 1) snprintf(out, n, "%s%sx%d ", ts, tg, e - s);
    else snprintf(out, n, "%s%s ", ts, tg);
    return e;
}

void burst_targets(const struct burst *b, const int *rank, char *out, size_t n)
{
    int byrank[64];
    for (int k = 0; k < b->naddr; k++) byrank[rank[k] - 1] = k;
    size_t o = 0;
    out[0] = 0;
    for (int r = 0; r < b->naddr && o + 1 < n; r++)
        o += snprintf(out + o, n - o, "%s%d/%s", r ? " " : "", r + 1, ac_name(b->addrs[byrank[r]]));
}

/* drain the output queue through stage 2; members of burst `id` are
 * collected into rec when valid */
static void drain(int id, int code, bool valid, struct burst *rec)
{
    while (ot < oh) {
        struct interrog *y = &outq[ot++ % OUTQ];
        bool member = id >= 0 && y->burst == id;
        /* plain A/C (no P4) never comes from a Mode S radar: those send A/C-only
         * all-calls instead, so a plain A/C query in a coded burst is another radar */
        if (member && (!valid || (code >= 0 && y->kind <= K_MODED))) { y->burst = -1; member = false; }
        bool pass = stage2(y, member ? code : -1);
        if (!valid || !rec) continue;
        if (pass && member) {
            int ka = rec_add(rec, y);
            graph_add(rec, y, GP_BEAM, ka);
        }
    }
}

static void close_burst(void)
{
    struct open_burst *o = &AS.ob;
    if (!o->open) return;
    o->open = false;
    int code = -1, bn = 0;
    for (int c = 0; c < NCODES; c++) if (o->votes[c] > bn) { bn = o->votes[c]; code = c; }
    /* temporary: a code needs BURST_CODE_MIN queries carrying it, so one stray
     * query of another radar heard during a burst does not name it */
    if (bn < BURST_CODE_MIN) code = -1;
    /* IIS 0 alone: an II0 interrogator (SI radars send their code in UF11 and DI 3 too) */
    if (bn == 0 && o->iis0 >= BURST_CODE_MIN) code = 0;
    bool valid = o->n_raw >= AS.cfg.burst_min;

    static struct burst rec;
    memset(&rec, 0, sizeof(rec));
    rec.id = o->id;
    rec.code = code;
    rec.partial = o->partial;
    drain(o->id, code, valid, &rec);
    if (valid && (AS.cfg.radar_sel < 0 || code == AS.cfg.radar_sel)) AS.last_burst_end = o->t1;
    /* radar filter: only that radar's bursts. Another radar's burst can hold a
     * query that passed it (a UF11 of the filtered radar heard in between, or
     * an all-call paired with one), but the burst is not the filtered radar's */
    if (valid && rec.n > 0 && (AS.cfg.radar_sel < 0 || code == AS.cfg.radar_sel)) {
        double margin = o->t1 - o->t0;
        if (margin < 0.03) margin = 0.03;
        if (margin > 0.2) margin = 0.2;
        double half = side_window(code);
        graph_context(&rec, o->t0, o->t1, margin, half);
        commit(&rec);
        /* keep collecting late queries for this burst; drop expired windows */
        int w = 0;
        for (int p = 0; p < nposts; p++) if (posts[p].until_side >= o->t1) posts[w++] = posts[p];
        nposts = w;
        if (nposts == NPOST) { memmove(posts, posts + 1, sizeof(posts[0]) * (NPOST - 1)); nposts--; }
        posts[nposts].idx = AS.nhist - 1;
        posts[nposts].until_side = o->t1 + half;
        posts[nposts].until_weak = o->t1 + margin;
        nposts++;
    }
}

static void open_burst(const struct interrog *x)
{
    struct open_burst *o = &AS.ob;
    memset(o, 0, sizeof(*o));
    o->open = true;
    o->id = next_burst_id++;
    o->t0 = o->t1 = x->t;
    o->partial = x->t - AS.gap_end < 0.1;
}

/* does a chain of members starting at pending index i reach a core one? */
static bool chain_reaches_core(uint64_t i)
{
    double last = pend[i % PENDING].t;
    for (uint64_t j = i; j < ph; j++) {
        struct interrog *p = &pend[j % PENDING];
        if (p->t - last > AS.cfg.burst_gap) return false;
        if (!is_member(p)) continue;
        if (is_core(p)) return true;
        last = p->t;
    }
    return false;
}

static void stage1(struct interrog *x, uint64_t idx)
{
    struct open_burst *o = &AS.ob;
    if (rx_match) rx_match(x);
    if (!x->tcas) {
        if (o->open && x->t - o->t1 > AS.cfg.burst_gap) close_burst();
        if (is_member(x)) {
            if (!o->open && chain_reaches_core(idx)) open_burst(x);
            if (o->open) {
                x->burst = o->id;
                o->n_raw++;
                o->t1 = x->t;
                int c = modes_own_code(x);
                if (c >= 0 && c < NCODES) o->votes[c]++;
                else if (modes_iis0(x)) o->iis0++;
            }
        }
    }
    if (o->open) {
        outq[oh++ % OUTQ] = *x;
        if (oh - ot >= OUTQ - 1) close_burst();
    } else {
        drain(-1, -1, false, NULL);
        stage2(x, -1);
    }
}

void analysis_tick(double now)
{
    if (now < 1e299) AS.now = now;
    bool flush = now >= 1e299;
    /* with 1090 replies, a query waits until its replies are in */
    while (pt < ph && pend[pt % PENDING].t < now - LOOKAHEAD &&
           (flush || !rx_ready || rx_ready(pend[pt % PENDING].t))) {
        uint64_t idx = pt++;
        stage1(&pend[idx % PENDING], idx);
    }
    double front = now - LOOKAHEAD;            /* everything before this went through stage 1 */
    if (pt < ph && pend[pt % PENDING].t < front) front = pend[pt % PENDING].t;
    if (AS.ob.open && front - AS.ob.t1 > AS.cfg.burst_gap) close_burst();
    if (now >= 1e299) { close_burst(); drain(-1, -1, false, NULL); }
}

void analysis_gap(double t0, double t1)
{
    analysis_tick(t0 + LOOKAHEAD + 1e-6);
    if (AS.ob.open) { AS.ob.partial = true; close_burst(); }
    drain(-1, -1, false, NULL);
    for (int r = 0; r < AS.nradars; r++) AS.radars[r].gap_since = true;
    AS.gap_end = t1;
}

/* ---------------- text summary ---------------- */

void analysis_print_summary(FILE *f)
{
    fprintf(f, "\n=== ssr1030 summary: %.1f s analysed, %ld interrogations, %ld passed filters ===\n",
            AS.now, AS.n_all, AS.n_pass);
    if (AS.filt_radar || AS.filt_modes || AS.filt_ac || AS.filt_noac || AS.filt_rollcall || AS.filt_aircraft || AS.filt_ii0)
        fprintf(f, "filtered out: other/unknown radar %ld, non-Mode S %ld, not plain A/C %ld, plain A/C %ld, "
                   "all-call/unaddressed %ld, other aircraft %ld, II0 %ld\n",
                AS.filt_radar, AS.filt_modes, AS.filt_ac, AS.filt_noac, AS.filt_rollcall, AS.filt_aircraft, AS.filt_ii0);
    if (aircraft_filter_on()) {
        char l[512];
        ac_list_str(&AS.cfg, true, l, sizeof(l));
        fprintf(f, "aircraft filter: %s%s\n", l, strstr(l, "(?)") ? "   ((?) = callsign not resolved yet)" : "");
    }
    fprintf(f, "background (sidelobe) level %.1f dBFS, burst member >= %.1f dBFS, core >= %.1f dBFS\n\n",
            AS.bg_db, AS.bg_db + AS.cfg.member_db, AS.bg_db + AS.cfg.core_db);
    fprintf(f, "%-26s %9s %9s %9s\n", "Radar interrogations", "total", "in-beam", "sidelobe");
    long all = krow_all(), oth[3] = { 0 };
    for (int i = 0; i < nkrows; i++) {
        const struct krow *r = &krows[i];
        long c[3] = { krow_count(r, AS.total, AS.commb), krow_count(r, AS.inbeam, AS.commb_inbeam),
                      krow_count(r, AS.side, AS.commb_side) };
        if (!c[0]) continue;
        if (!krow_shown(r, c[0], all)) { for (int j = 0; j < 3; j++) oth[j] += c[j]; continue; }
        fprintf(f, "%-26s %9ld %9ld %9ld\n", r->label ? r->label : kind_label(r->k), c[0], c[1], c[2]);
    }
    if (oth[0]) fprintf(f, "%-26s %9ld %9ld %9ld\n", KROW_OTHER, oth[0], oth[1], oth[2]);
    fprintf(f, "\nTCAS/ACAS (excluded): %ld  (UF0 %ld, UF16 %ld, whisper-shout %ld; %d aircraft)\n",
            AS.tcas_total, AS.tcas_uf0, AS.tcas_uf16, AS.tcas_ws, AS.naddr_tcas);

    fprintf(f, "\nRadars (by II/SI code):\n");
    for (int r = 0; r < AS.nradars; r++) {
        struct radar *R = &AS.radars[r];
        int np = R->np < 64 ? R->np : 64;
        fprintf(f, "  %-6s bursts %d (+%d partial)  rotation %.3f s  mean dwell %.1f ms  mean %.1f interrog/burst\n",
                R->key, R->nbursts, R->npartial, np ? median_d(R->periods, np) : 0.0,
                R->nbursts ? 1e3 * R->dwell_sum / R->nbursts : 0,
                R->nbursts ? (double)R->total_sum / R->nbursts : 0);
        if (!R->nbursts) continue;
        fprintf(f, "         per burst (avg [min-max]):");
        for (int g = 0; g < G_COUNT; g++) {
            if (!R->cnt_max[g]) continue;
            fprintf(f, " %s %.1f [%d-%d]", group_label(g), (double)R->cnt_sum[g] / R->nbursts,
                    R->cnt_min[g], R->cnt_max[g]);
        }
        fputc('\n', f);
    }

    fprintf(f, "\nBursts (last %d):\n", AS.nhist < BURST_HIST ? AS.nhist : BURST_HIST);
    int first = AS.nhist > BURST_HIST ? AS.nhist - BURST_HIST : 0;
    for (int i = first; i < AS.nhist; i++) {
        struct burst *b = &AS.hist[i % BURST_HIST];
        fprintf(f, "  #%-4d t=%9.3f %-5s%s dur %6.1f ms peak %5.1f dBFS n=%-4d", b->id, b->t0, b->key,
                b->partial ? "*" : " ", 1e3 * (b->t1 - b->t0), b->peak, b->n);
        if (b->period > 0) fprintf(f, " P=%.3fs", b->period);
        for (int g = 0; g < G_COUNT; g++) if (b->gcnt[g]) fprintf(f, " %s:%d", group_label(g), b->gcnt[g]);
        fprintf(f, "\n     ");
        int col = 5, rank[64];
        burst_ranks(b, rank);
        for (int s = 0; s < b->nks;) {
            char tok[48];
            s = burst_item(b, s, rank, tok, sizeof(tok));
            if (col + (int)strlen(tok) > 110) { fprintf(f, "\n     "); col = 5; }
            fputs(tok, f);
            col += strlen(tok);
        }
        if (b->naddr) {
            char tl[1024];
            burst_targets(b, rank, tl, sizeof(tl));
            fprintf(f, "\n     targets: %s", tl);
        }
        fputc('\n', f);
    }
    if (AS.nhist) {
        fprintf(f, "  (* = burst touched an ADS-B window; excluded from per-burst statistics)\n");
        fprintf(f, "  tokens: A C A/C-only all-call | AA CC plain Mode A/C | As Cs A/C/S all-call | S UF11 all-call\n"
                   "          %s | S-msg Comm-A (UF20/21) | S-elm UF24\n",
                bds_human ? "alt UF4 altitude | sqk UF5 Mode A code" : "UF4 altitude roll-call | UF5 Mode A code roll-call");
        if (bds_human)
            fprintf(f, "          Comm-B: acid 2,0 aircraft identification (callsign) | selAlt 4,0 selected alt | trk 5,0 track/turn |\n"
                       "          hdgSpd 6,0 heading/speed | cap 1,0 | gicb 1,7 | pos 5,1/5,2 | met 4,4 | S-bXY other\n");
        else
            fprintf(f, "          Comm-B: UFnbXY = UFn roll-call requesting BDS X,Y (UF4b20 callsign, b40 selected alt,\n"
                       "          b50 track/turn, b60 heading/speed); -h shows short names\n");
    }

    fprintf(f, "\nIn-beam transitions (top 16):\n");
    {
        struct trans top[16];
        int nt = analysis_top_transitions(top, 16);
        for (int i = 0; i < nt; i++) {
            char ta[24], tb[24];
            fprintf(f, "  %-10s -> %-10s %7ld  dt avg %8.3f ms  min %8.3f  max %8.3f\n",
                    seq_token_str(top[i].a, ta, sizeof(ta)), seq_token_str(top[i].b, tb, sizeof(tb)), top[i].n,
                    1e3 * top[i].sum / top[i].n, 1e3 * top[i].min, 1e3 * top[i].max);
        }
    }

    fprintf(f, "\nUF11 all-call codes:");
    for (int si = 0; si < 2; si++)
        for (int c = 0; c < 64; c++)
            for (int p = 0; p < 16; p++)
                if (AS.uf11_code[si][c][p]) fprintf(f, " %s%d/PR%d:%ld", si ? "SI" : "II", c, p, AS.uf11_code[si][c][p]);
    fprintf(f, "\nComm-B requests (BDS):");
    for (int b = 0; b < 256; b++) if (AS.bds[b]) fprintf(f, " %X,%X:%ld", b >> 4, b & 15, AS.bds[b]);
    fprintf(f, "\nDI values:");
    for (int k = 0; k < 8; k++) if (AS.di[k]) fprintf(f, " DI%d:%ld", k, AS.di[k]);
    fprintf(f, "\nPC values:");
    for (int k = 0; k < 8; k++) if (AS.pc[k]) fprintf(f, " PC%d:%ld", k, AS.pc[k]);
    fprintf(f, "\nRoll-called aircraft: %d distinct addresses (callsigns known: %d)\n", AS.naddr_radar, cs_count());
    int shown = 0;
    for (int k = 0; k < ADDR_HASH && shown < 60; k++) {
        struct addr_entry *e = &AS.addrs[k];
        if (!e->n || e->tcas) continue;
        fprintf(f, "  %-8s %5d", ac_name(e->addr), e->n);
        if (++shown % 5 == 0) fputc('\n', f);
    }
    fputc('\n', f);
}
