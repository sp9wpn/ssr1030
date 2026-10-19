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
#include "analysis.h"

#define LOOKAHEAD   0.3
#define BIN_S       0.01
#define NBINS       2000
#define PENDING     65536
#define OUTQ        65536
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
static void (*sinkf)(const struct interrog *it);
static double last11_t = -1;
static int last11_code = -1;
static unsigned target_ver = ~0u;

static const char *group_labels[G_COUNT] = {
    "AA/CC", "A/C", "ACSall", "S", "alt", "sqk", "acid", "setAlt", "trk", "hdgSpd", "oBDS",
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

const char *group_label(int g) { return group_labels[g]; }

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
    AS.target = cfg->ac_addr;
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

double analysis_last_burst_end(void) { return AS.last_burst_end; }
int64_t analysis_target_addr(void) { return AS.target; }

/* ---------------- filters ---------------- */

static bool aircraft_filter_on(void) { return AS.cfg.ac_addr >= 0 || AS.cfg.ac_callsign[0]; }

static void resolve_target(void)
{
    if (AS.cfg.ac_addr >= 0 || !AS.cfg.ac_callsign[0]) return;
    if (cs_version() == target_ver) return;
    target_ver = cs_version();
    AS.target = cs_find(AS.cfg.ac_callsign);
}

static bool directed(int k)
{
    return k == K_UF0 || k == K_UF4 || k == K_UF5 || k == K_UF16 || k == K_UF20 || k == K_UF21 ||
           k == K_UF24 || k == K_UF19;
}

static bool to_target(const struct interrog *x)
{
    return AS.target >= 0 && directed(x->kind) && x->addr == (uint32_t)AS.target;
}

static bool passes(const struct interrog *x)
{
    if (AS.cfg.radar_sel >= 0 && x->radar != AS.cfg.radar_sel) { AS.filt_radar++; return false; }
    if (AS.cfg.modes_only && !kind_is_modes(x->kind)) { AS.filt_modes++; return false; }
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

/* returns true when the interrogation passed the filters */
static bool stage2(struct interrog *x, int burst_code)
{
    AS.n_all++;
    resolve_target();
    attribute(x, burst_code);

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
        if (modes_bds(x->msg, &b1, &b2)) AS.bds[b1 * 16 + b2]++;
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

static void rec_add(struct burst *b, const struct interrog *x)
{
    if (b->n == 0) { b->t0 = x->t; b->peak = -99; }
    b->n++;
    b->cnt[x->kind]++;
    b->t1 = x->t;
    if (x->amp_db > b->peak) { b->peak = x->amp_db; b->tpeak = x->t; }
    uint16_t tok = seq_token(x);
    b->gcnt[tok_group(tok)]++;
    if (b->nks < SEQ_MAX) b->ks[b->nks++] = tok;
    if (directed(x->kind)) {
        int k;
        for (k = 0; k < b->naddr; k++) if (b->addrs[k] == x->addr) break;
        if (k == b->naddr && b->naddr < 64) b->addrs[b->naddr++] = x->addr;
    }
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
    AS.hist[AS.nhist % BURST_HIST] = *b;
    AS.nhist++;
}

/* drain the output queue through stage 2; members of burst `id` are
 * collected into rec when valid */
static void drain(int id, int code, bool valid, struct burst *rec)
{
    while (ot < oh) {
        struct interrog *y = &outq[ot++ % OUTQ];
        bool member = id >= 0 && y->burst == id;
        if (member && !valid) y->burst = -1;
        bool pass = stage2(y, member && valid ? code : -1);
        if (pass && member && valid) rec_add(rec, y);
    }
}

static void close_burst(void)
{
    struct open_burst *o = &AS.ob;
    if (!o->open) return;
    o->open = false;
    int code = -1, bn = 0;
    for (int c = 0; c < NCODES; c++) if (o->votes[c] > bn) { bn = o->votes[c]; code = c; }
    bool valid = o->n_raw >= AS.cfg.burst_min;

    static struct burst rec;
    memset(&rec, 0, sizeof(rec));
    rec.id = o->id;
    rec.code = code;
    rec.partial = o->partial;
    drain(o->id, code, valid, &rec);
    if (valid && (AS.cfg.radar_sel < 0 || code == AS.cfg.radar_sel)) AS.last_burst_end = o->t1;
    if (valid && rec.n > 0) commit(&rec);
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
    while (pt < ph && pend[pt % PENDING].t < now - LOOKAHEAD) {
        uint64_t idx = pt++;
        stage1(&pend[idx % PENDING], idx);
    }
    if (AS.ob.open && now - LOOKAHEAD - AS.ob.t1 > AS.cfg.burst_gap) close_burst();
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
    if (AS.filt_radar || AS.filt_modes || AS.filt_aircraft)
        fprintf(f, "filtered out: other/unknown radar %ld, non-Mode S %ld, other aircraft %ld\n",
                AS.filt_radar, AS.filt_modes, AS.filt_aircraft);
    if (aircraft_filter_on()) {
        if (AS.target >= 0) fprintf(f, "aircraft filter: %06llX %s\n", (unsigned long long)AS.target,
                                    cs_lookup((uint32_t)AS.target) ? cs_lookup((uint32_t)AS.target) : "");
        else fprintf(f, "aircraft filter: callsign %s not resolved (no ADS-B identification seen)\n",
                     AS.cfg.ac_callsign);
    }
    fprintf(f, "background (sidelobe) level %.1f dBFS, burst member >= %.1f dBFS, core >= %.1f dBFS\n\n",
            AS.bg_db, AS.bg_db + AS.cfg.member_db, AS.bg_db + AS.cfg.core_db);
    fprintf(f, "%-26s %9s %9s %9s\n", "Radar interrogations", "total", "in-beam", "sidelobe");
    for (int k = 0; k < K_COUNT; k++) {
        if (!AS.total[k]) continue;
        fprintf(f, "%-26s %9ld %9ld %9ld\n", kind_label(k), AS.total[k], AS.inbeam[k], AS.side[k]);
    }
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
        int col = 5;
        for (int s = 0; s < b->nks;) {
            int e = s;
            while (e < b->nks && b->ks[e] == b->ks[s]) e++;
            char tok[48], ts[24];
            seq_token_str(b->ks[s], ts, sizeof(ts));
            if (e - s > 1) snprintf(tok, sizeof(tok), "%s#%d ", ts, e - s);
            else snprintf(tok, sizeof(tok), "%s ", ts);
            if (col + (int)strlen(tok) > 110) { fprintf(f, "\n     "); col = 5; }
            fputs(tok, f);
            col += strlen(tok);
            s = e;
        }
        if (b->naddr) {
            fprintf(f, "\n     targets:");
            for (int a = 0; a < b->naddr; a++) {
                fprintf(f, " %s", ac_name(b->addrs[a]));
            }
        }
        fputc('\n', f);
    }
    if (AS.nhist) {
        fprintf(f, "  (* = burst touched an ADS-B window; excluded from per-burst statistics)\n");
        fprintf(f, "  tokens: A C A/C-only all-call | AA CC plain Mode A/C | As Cs A/C/S all-call | S UF11 all-call\n"
                   "          alt UF4 altitude | sqk UF5 Mode A code | S-msg Comm-A (UF20/21) | S-elm UF24\n"
                   "          Comm-B: acid 2,0 aircraft identification (callsign) | setAlt 4,0 selected alt | trk 5,0 track/turn |\n"
                   "          hdgSpd 6,0 heading/speed | cap 1,0 | gicb 1,7 | pos 5,1/5,2 | met 4,4 | S-bXY other\n");
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
