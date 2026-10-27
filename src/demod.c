/* Pulse detection, interrogation classification and Mode S DPSK demodulation.
 *
 * Signal chain per sample: DC removal -> magnitude -> threshold segmentation
 * (with short-gap merging so DPSK phase-reversal dips do not split P6) ->
 * half-amplitude edge refinement -> pulse queue -> pattern matcher.
 *
 * Two segmenters: the main one (threshold --threshold) gives the "strong"
 * pulses; a second one with the lower Mode S threshold adds pulses that never
 * reach the main threshold, so interrogations heard through the radar's
 * sidelobes are decoded too (DPSK with 56+ bits and parity is hard to fake
 * from noise). Mode A/C, intermode, S1 and the pulse density use only strong
 * pulses: a P1-P3 pair of weak noise pulses would be too common. Pulses above
 * both thresholds appear twice (one copy per segmenter); Mode S takes either
 * and then marks both used. Near the main threshold the main segmenter cuts
 * a pulse short or breaks P6 into fragments, the weak copy is whole. The weak
 * segmenter also ends a segment below SEG_REL of its peak, so the edges of
 * strong pulses do not bridge into one weak segment.
 *
 * Interrogation timing (ICAO Annex 10 Vol IV, from P1 leading edge):
 *   Mode A/C           P3 at 8/21 us, P2 (SLS) at 2 us
 *   intermode          P4 2 us after P3: 0.8 us = A/C-only, 1.6 us = A/C/S all-call
 *   whisper-shout      S1 2 us before P1
 *   Mode S             P2 at 2 us, P6 at 3.5 us (16.25 / 30.25 us long),
 *                      sync phase reversal 1.25 us into P6; data reversals at
 *                      sync + N x 0.25 us, N >= 2; a reversal = 1. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ssr1030.h"

#define BUSY_GAP_US 21.0           /* longest gap inside an interrogation: Mode C P1-P3 is 20.2 us */
#define BUSY_MIN_US 0.6            /* busy pulses: interrogation pulses are 0.8 us (+-0.1) */
#define BUSY_MAX_HZ 0.5e6          /* ... and within this of 1030 MHz */
#define RING_BITS   20
#define RING        (1u << RING_BITS)
#define RMASK       (RING - 1)
#define QN          4096
#define QMASK       (QN - 1)
#define MAX_PULSE_US 40.0
#define SEG_REL     0.3f

struct demod {
    double fs, spus;            /* samples per second / per microsecond */
    float complex *z;           /* ring of DC-corrected samples */
    float *m;                   /* ring of magnitudes */
    uint64_t n;                 /* global index of next sample */

    float complex dc;
    float noise;                /* mean noise magnitude */
    float thr_db, thr_abs, thr;
    float modes_db, thr_lo;     /* Mode S pulse threshold (<= thr) */
    float thr_sls;              /* P2/P5 presence at their known positions: noise + 6 dB */
    double nacc; uint64_t ncnt;

    struct seg { int on; uint64_t ps, plast; float peak; } hi, lo;   /* main / Mode S segmenter */
    double sleads[QN];          /* leads of recent strong pulses, for the density */
    uint64_t ns;
    int merge_gap, min_len, max_len;

    struct pulse q[QN];
    uint64_t qhead, qtail;      /* qhead: next write; qtail: oldest unexamined */
    uint64_t qdens;             /* oldest strong pulse inside the 1 ms density window */
    int dens_max;
    uint64_t garbled;

    interrog_cb cb;
    void *ctx;
    uint64_t seq;

    uint64_t clipped, pulses, stray;
    uint64_t rx, busy;          /* samples received on 1030 MHz; of them busy (see busy_gap) */
    uint64_t busy_end;          /* end (sample index) of the last pulse counted as busy */
    uint64_t busy_gap;          /* samples: a gap this short between pulses is busy too */
    float peak;
    float *tmp;
};

struct demod *demod_new(double fs, float thr_db, float thr_abs, float modes_db, int dens_max, interrog_cb cb, void *ctx)
{
    struct demod *d = calloc(1, sizeof(*d));
    d->fs = fs;
    d->modes_db = modes_db > 0 && modes_db < thr_db ? modes_db : thr_db;
    d->spus = fs / 1e6;
    d->z = calloc(RING, sizeof(float complex));
    d->m = calloc(RING, sizeof(float));
    d->thr_db = thr_db;
    d->thr_abs = thr_abs;
    d->dens_max = dens_max;
    d->noise = -1;
    d->busy_gap = (uint64_t)lround(BUSY_GAP_US * d->spus);
    d->merge_gap = (int)lround(0.25 * d->spus);
    if (d->merge_gap < 2) d->merge_gap = 2;
    d->min_len = (int)lround(0.3 * d->spus);
    d->max_len = (int)lround(MAX_PULSE_US * d->spus);
    d->tmp = malloc(sizeof(float) * (d->max_len + 16));
    d->cb = cb;
    d->ctx = ctx;
    return d;
}

static inline float M(struct demod *d, uint64_t i) { return d->m[i & RMASK]; }
static inline float complex Z(struct demod *d, uint64_t i) { return d->z[i & RMASK]; }

static float dbfs(float a) { return a > 0 ? 20.0f * log10f(a / 127.0f) : -99.0f; }

static int cmpf(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

/* ---------------- Mode S DPSK ---------------- */

/* Normalised correlation of the phase before and after boundary position t
 * (in samples). -1 means a clean 180 degree reversal, +1 no reversal. */
static float reversal_corr(struct demod *d, double t, int h, float complex *prod)
{
    uint64_t c = (uint64_t)llround(t);
    float complex a = 0, b = 0;
    for (int k = 1; k <= h; k++) {
        a += Z(d, c + k);
        b += Z(d, c - k);
    }
    float complex p = a * conjf(b);
    if (prod) *prod = p;
    float mag = cabsf(p);
    return mag > 0 ? crealf(p) / mag : 0;
}

/* Around the sync reversal the signal is P5 (constant, if present) plus P6
 * whose sign flips once. Least-squares fit over +-0.35 us: for each split,
 * c = (mean_before + mean_after) / 2 is P5 and (mean_before - mean_after) / 2
 * is P6; the split with the smallest residual is the reversal. Samples next
 * to the split (the band-limited transition) are left out. */
static void p5_fit(struct demod *d, double center, float *p5, float *p6, double *tau)
{
    int w = (int)lround(0.35 * d->spus), g = (int)lround(0.08 * d->spus);
    if (g < 1) g = 1;
    uint64_t c0 = (uint64_t)llround(center);
    float best = 1e30f;
    *p5 = 0; *p6 = 0; *tau = center;
    for (int s = -w + 2 * g + 1; s <= w - 2 * g - 1; s++) {
        float complex mb = 0, ma = 0;
        int nb = 0, na = 0;
        for (int k = -w; k <= s - g; k++) { mb += Z(d, c0 + k); nb++; }
        for (int k = s + g; k <= w; k++) { ma += Z(d, c0 + k); na++; }
        if (nb < 2 || na < 2) continue;
        mb /= nb; ma /= na;
        float r = 0;
        for (int k = -w; k <= s - g; k++) { float complex e = Z(d, c0 + k) - mb; r += crealf(e * conjf(e)); }
        for (int k = s + g; k <= w; k++) { float complex e = Z(d, c0 + k) - ma; r += crealf(e * conjf(e)); }
        r /= nb + na;
        if (r < best) {
            best = r;
            *p5 = cabsf(mb + ma) / 2;
            *p6 = cabsf(mb - ma) / 2;
            *tau = (double)c0 + s;
        }
    }
}

/* P6 level: 70th percentile of the plateau magnitude, leaving out the
 * edges and the P5 window around the sync (centre ts); the percentile
 * skips the dips at phase reversals. */
static float p6_level(struct demod *d, double p6lead, double p6width, double ts)
{
    uint64_t a = (uint64_t)llround(p6lead + 0.3 * d->spus);
    uint64_t b = (uint64_t)llround(p6lead + (p6width - 0.3) * d->spus);
    double x0 = ts - 0.45 * d->spus, x1 = ts + 0.45 * d->spus;
    int n = 0;
    for (uint64_t k = a; k < b && n < d->max_len; k++) {
        if ((double)k > x0 && (double)k < x1) continue;
        d->tmp[n++] = M(d, k);
    }
    if (n < 8) return 0;
    qsort(d->tmp, n, sizeof(float), cmpf);
    return d->tmp[(n * 7) / 10];
}

static void decode_modes(struct demod *d, struct interrog *it, double p6lead, double p6width)
{
    double chip = 0.25 * d->spus;
    int h = (int)(chip / 2);
    if (h < 1) h = 1;

    /* Locate the sync phase reversal near P6 + 1.25 us. It is the first
     * reversal in P6 and data reversals follow at >= 0.5 us, so within
     * +-0.35 us the earliest clear reversal (local minimum below -0.5) is
     * the sync; otherwise fall back to the deepest point. */
    double exp_sync = p6lead + 1.25 * d->spus;
    double best_t = exp_sync;
    float best = 2;
    int span = (int)lround(0.35 * d->spus);
    float cprev = 2;
    bool found = false;
    for (int o = -span; o <= span + 1; o++) {
        float c = reversal_corr(d, exp_sync + o, h, NULL);
        if (o > -span && cprev < -0.5f && cprev <= c) {     /* local minimum at o - 1 */
            best = cprev; best_t = exp_sync + o - 1; found = true;
            break;
        }
        cprev = c;
    }
    if (!found) {
        for (int o = -span; o <= span; o++) {
            float c = reversal_corr(d, exp_sync + o, h, NULL);
            if (c < best) { best = c; best_t = exp_sync + o; }
        }
    }
    it->sync = best;

    /* P5 (Mode S SLS) overlaps the sync reversal: fitted around the nominal
     * sync, see p5_fit(). The fitted split is the reversal itself, so when P6
     * is clearly present and stronger than P5 it also gives the sync, even if
     * the reversal falls between samples and correlates weakly. */
    {
        float p5, p6;
        double tau;
        p5_fit(d, exp_sync, &p5, &p6, &tau);
        float p6lvl = p6_level(d, p6lead, p6width, tau);
        it->p6_db = p6lvl > 0 ? dbfs(p6lvl) : NAN;
        it->p5_db = p5 >= d->thr_sls && p6lvl > 0 ? dbfs(p5) - it->p6_db : NAN;   /* relative to P6 */
        if (isnan(it->p5_db) && p5 > 0) it->sls_floor_db = dbfs(p5);
        if (!isnan(it->p5_db) && p5 >= p6) {
            /* transponder rule: P5 at least as strong as P6 = sidelobe. A
             * transponder would not reply, but the data chips after P5 are intact */
            it->sidelobe = true;
            best = 2;
            exp_sync = tau;
        } else if (p6 >= d->thr && best > -0.5f) {
            best_t = tau;
            best = -0.5f - 1e-3f;                     /* sync taken from the fit */
            it->sync = best;
        }
    }

    int len_class = 0;
    if (p6width > 14.5 && p6width < 18.0) len_class = 56;
    else if (p6width > 28.5 && p6width < 32.5) len_class = 112;
    int maxbits = (p6width > 18.0) ? 112 : 56;

    if (best > -0.5f) {
        if (isnan(it->p5_db)) {             /* no sync and no P5: garbled */
            it->kind = K_S_GARBLE;
            return;
        }
        /* Sync masked by P5: recover the bit timing from the data chips. P5
         * covers only sync +-0.4 us. A shift by a whole chip decodes just as
         * cleanly, so search only +-half a chip around the nominal sync (the
         * P6 leading edge, before P5 starts, fixes it that well) and keep the
         * offset whose chips (from the 2nd bit on) give the most confident
         * decisions. */
        it->sidelobe = true;
        int nb = maxbits < 56 ? maxbits : 56;
        int half = (int)(chip / 2);
        float score = 0;
        for (int o = -half; o <= half; o++) {
            float sum = 0;
            for (int k = 1; k < nb; k++) sum += fabsf(reversal_corr(d, exp_sync + o + (k + 2) * chip, h, NULL));
            sum /= nb - 1;
            if (sum > score) { score = sum; best_t = exp_sync + o; }
        }
        if (score < 0.6f) {                 /* chips not clean enough to trust */
            it->kind = K_S_NOSYNC;
            return;
        }
    }

    /* Data reversals occur N x 0.25 us after sync with N >= 2 (Annex 10
     * Vol IV 3.1.2.1.5.2.3), so data bit k (0-based) sits at sync + (k + 2) chips. */
    float conf[MODES_MAX_BITS];
    memset(it->msg, 0, sizeof(it->msg));
    for (int k = 0; k < maxbits; k++) {
        float c = reversal_corr(d, best_t + (k + 2) * chip, h, NULL);
        if (c < 0) it->msg[k / 8] |= 0x80 >> (k % 8);
        conf[k] = fabsf(c);
    }

    int uf = modes_uf(it->msg);
    int nb = modes_uf_bits(uf);
    if (len_class && nb != len_class) {
        it->len_mismatch = true;
        nb = len_class;
    }
    if (nb > maxbits) nb = maxbits;
    if (nb == 56) memset(it->msg + 7, 0, 7);
    float qmin = 1;
    for (int k = 0; k < nb; k++) if (conf[k] < qmin) qmin = conf[k];
    it->nbits = nb;
    it->uf = uf;
    it->quality = qmin;
    it->addr = modes_ap_to_addr(it->msg, nb);

    switch (uf) {
    case 0:  it->kind = K_UF0;  it->tcas = true; break;
    case 4:  it->kind = K_UF4;  break;
    case 5:  it->kind = K_UF5;  break;
    case 11:
        /* all-calls carry address FFFFFF in AP: anything else is a corrupted decode */
        it->kind = it->addr == 0xFFFFFF ? K_UF11 : K_S_GARBLE;
        break;
    case 16: it->kind = K_UF16; it->tcas = true; break;
    case 19: it->kind = K_UF19; break;
    case 20: it->kind = K_UF20; break;
    case 21: it->kind = K_UF21; break;
    case 24: it->kind = K_UF24; break;
    default: it->kind = K_UF_OTHER; break;
    }
}

/* ---------------- pattern matcher ---------------- */

/* median magnitude between two sample positions */
static float slot_level(struct demod *d, double a, double b)
{
    int n = 0;
    for (uint64_t k = (uint64_t)llround(a); k <= (uint64_t)llround(b) && n < 64; k++) d->tmp[n++] = M(d, k);
    if (!n) return 0;
    qsort(d->tmp, n, sizeof(float), cmpf);
    return d->tmp[n / 2];
}

static inline double rel_us(struct demod *d, const struct pulse *a, const struct pulse *b)
{
    return (b->lead - a->lead) / d->spus;
}

static bool is_short(const struct pulse *p) { return p->width > 0.6 && p->width < 1.05; }

/* find unused pulse after index i at offset off +- tol us; strong pulses only
 * for Mode A/C (weak ones are Mode S candidates) */
static int find_at(struct demod *d, uint64_t i, double off, double tol, bool need_short, bool need_strong)
{
    struct pulse *p1 = &d->q[i & QMASK];
    for (uint64_t j = i + 1; j < d->qhead; j++) {
        struct pulse *p = &d->q[j & QMASK];
        double r = rel_us(d, p1, p);
        if (r > off + tol) break;
        if (p->used || r < off - tol) continue;
        if (need_short && !is_short(p)) continue;
        if (need_strong && !p->strong) continue;
        return (int)(j - i);
    }
    return -1;
}

static void emit(struct demod *d, struct interrog *it, const struct pulse *p1)
{
    it->seq = d->seq++;
    it->t = p1->lead / d->fs;
    it->p1_db = dbfs(p1->amp);
    it->amp_db = isnan(it->p6_db) ? it->p1_db : it->p6_db;   /* Mode S: P6 (longer, less noisy) */
    it->burst = -1;
    d->cb(it, d->ctx);
}

static void match(struct demod *d, uint64_t i)
{
    struct pulse *p1 = &d->q[i & QMASK];
    if (p1->used) return;
    if (!is_short(p1)) { if (p1->strong) d->stray++; return; }

    struct interrog it;
    memset(&it, 0, sizeof(it));
    it.p2_db = NAN;
    it.p5_db = NAN;
    it.sls_floor_db = NAN;
    it.p6_db = NAN;
    it.noise_db = d->noise > 0 ? dbfs(d->noise) : NAN;

    /* --- Mode S: P1, P2 @2.0, P6 @3.5 --- */
    int j2 = find_at(d, i, 2.0, 0.2, true, false);
    if (j2 > 0) {
        int j6 = find_at(d, i, 3.5, 0.25, false, true);    /* the main segmenter's P6 when there is one */
        if (j6 < 0) j6 = find_at(d, i, 3.5, 0.25, false, false);
        if (j6 > 0) {
            struct pulse *p6 = &d->q[(i + j6) & QMASK];
            /* chain fragments of a broken-up P6 */
            double end = p6->lead + p6->width * d->spus;
            uint64_t last = i + j6;
            for (uint64_t k = i + j6 + 1; k < d->qhead; k++) {
                struct pulse *p = &d->q[k & QMASK];
                if (p->lead - end > 0.5 * d->spus) break;
                if (p6->strong && !p->strong && p->lead > end) continue;   /* echo or ringing after a strong P6 */
                double e = p->lead + p->width * d->spus;
                if (e > end) end = e;
                last = k;
            }
            double w6 = (end - p6->lead) / d->spus;
            if (w6 > 10.0) {
                struct pulse *p2 = &d->q[(i + j2) & QMASK];
                it.p2_db = 20.0f * log10f(p2->amp / p1->amp);
                decode_modes(d, &it, p6->lead, w6);
                p1->used = p2->used = true;
                for (uint64_t k = i + j6; k <= last; k++) d->q[k & QMASK].used = true;
                /* the other segmenter's copies of P1 and P2 */
                for (uint64_t k = i > d->qtail + 8 ? i - 8 : d->qtail; k < i + j6; k++) {
                    struct pulse *p = &d->q[k & QMASK];
                    if (fabs(p->lead - p1->lead) < 0.3 * d->spus || fabs(p->lead - p2->lead) < 0.3 * d->spus) p->used = true;
                }
                emit(d, &it, p1);
                return;
            }
        }
    }

    /* --- Mode A/C with optional P2, P4, S1. Modes 1/2/B/D are not detected:
     * on air they were almost all coincidences of unrelated pulses. Strong
     * pulses only (see the top of the file). --- */
    if (!p1->strong) return;
    static const double off3[] = { 8, 21 };
    static const int base[] = { K_MODEA, K_MODEC };
    for (int m = 0; m < 2; m++) {
        int j3 = find_at(d, i, off3[m], 0.15, true, true);
        if (j3 < 0) continue;
        struct pulse *p3 = &d->q[(i + j3) & QMASK];
        if (fabsf(20.0f * log10f(p3->amp / p1->amp)) > 3.0f) continue;
        if (p1->density > d->dens_max) { d->garbled++; return; }
        it.kind = base[m];
        int j2s = find_at(d, i, 2.0, 0.2, true, true);
        if (j2s > 0 && j2s < j3) j2 = j2s;
        if (j2 > 0 && j2 < j3 && d->q[(i + j2) & QMASK].strong) {
            struct pulse *p2 = &d->q[(i + j2) & QMASK];
            it.p2_db = 20.0f * log10f(p2->amp / p1->amp);
            if (it.p2_db >= 0) it.sidelobe = true;
            p2->used = true;
        } else {
            /* no P2 pulse detected: read the P2 slot directly (P2 is usually
             * well below the pulse threshold in the main beam) */
            float a2 = slot_level(d, p1->lead + 2.15 * d->spus, p1->lead + 2.65 * d->spus);
            if (a2 >= d->thr_sls) it.p2_db = 20.0f * log10f(a2 / p1->amp);
            else if (a2 > 0) it.sls_floor_db = dbfs(a2);
        }
        int j4 = find_at(d, i, off3[m] + 2.0, 0.15, false, true);
        if (j4 > 0) {
            struct pulse *p4 = &d->q[(i + j4) & QMASK];
            if (p4->width > 0.5 && p4->width < 2.3 &&
                fabsf(20.0f * log10f(p4->amp / p3->amp)) <= 3.0f) {
                it.p4_width = p4->width;
                bool lng = p4->width >= 1.25;
                if (m == 0) it.kind = lng ? K_ACS_A : K_ACONLY_A;
                else it.kind = lng ? K_ACS_C : K_ACONLY_C;
                p4->used = true;
            }
        }
        /* S1 (whisper-shout) 2 us before P1: look back */
        for (uint64_t k = i; k > 0 && k + QN > d->qhead + 8 && i - k < 8; k--) {
            struct pulse *p = &d->q[(k - 1) & QMASK];
            double r = rel_us(d, p, p1);
            if (r > 2.2) break;
            if (r > 1.8 && is_short(p) && p->strong) { it.s1 = true; break; }
        }
        if (it.s1) it.tcas = true;
        p1->used = p3->used = true;
        emit(d, &it, p1);
        return;
    }
    d->stray++;
}

static void assemble(struct demod *d, bool all)
{
    double lim = MAX_PULSE_US * d->spus;
    while (d->qtail < d->qhead) {
        struct pulse *p = &d->q[d->qtail & QMASK];
        if (!all && (double)d->n - p->lead < lim) break;
        match(d, d->qtail);
        d->qtail++;
    }
}

/* ---------------- pulse finalisation ---------------- */

static void finish_pulse(struct demod *d, const struct seg *g, bool strong)
{
    uint64_t s = g->ps, e = g->plast;
    int len = (int)(e - s + 1);
    if (len < d->min_len || len > d->max_len) return;

    /* plateau amplitude: median of the segment, ignoring the edges */
    int a0 = len > 6 ? 2 : 0, n = 0;
    for (int k = a0; k < len - a0; k++) d->tmp[n++] = M(d, s + k);
    qsort(d->tmp, n, sizeof(float), cmpf);
    float amp = d->tmp[n / 2];
    if (len > 24) amp = d->tmp[(n * 3) / 4]; /* long DPSK pulse: skip reversal dips */
    if (amp < 0.8f * (strong ? d->thr : d->thr_lo)) return;   /* merged noise blips, not a pulse */
    float hl = amp * 0.5f;

    /* leading edge at half amplitude */
    uint64_t i = s;
    if (M(d, i) >= hl) {
        int lim = (int)(1.0 * d->spus);
        while (lim-- > 0 && M(d, i - 1) >= hl) i--;
    } else {
        while (i < e && M(d, i) < hl) i++;
    }
    float m0 = M(d, i - 1), m1 = M(d, i);
    double lead = (double)i - 1 + (m1 > m0 ? (hl - m0) / (m1 - m0) : 1.0);

    /* trailing edge */
    uint64_t j = e;
    if (M(d, j) >= hl) {
        while (j + 1 < d->n && M(d, j + 1) >= hl) j++;
    } else {
        while (j > i && M(d, j) < hl) j--;
    }
    m0 = M(d, j); m1 = M(d, j + 1);
    double trail = (double)j + (m0 > m1 ? (m0 - hl) / (m0 - m1) : 0.0);

    struct pulse np = { .lead = lead, .width = (trail - lead) / d->spus, .amp = amp, .strong = strong };
    if (strong) {
        /* channel busy: each pulse, plus the gap from the previous pulse when it
         * is up to busy_gap (inside an interrogation, also a garbled or partial
         * one). Only pulses of interrogation width on 1030 MHz count: on air,
         * shorter blips (noise near the threshold) and pulses elsewhere in the
         * received band were most of the rest */
        uint64_t pe = e + 1;
        float complex rot = 0;
        for (uint64_t k = s + 1; k <= e; k++) rot += Z(d, k) * conjf(Z(d, k - 1));
        double off_hz = cargf(rot) * d->fs / (2 * M_PI);
        if (np.width >= BUSY_MIN_US && fabs(off_hz) < BUSY_MAX_HZ) {
            if (d->busy_end && s <= d->busy_end + d->busy_gap) { if (pe > d->busy_end) d->busy += pe - d->busy_end; }
            else d->busy += pe - s;
            if (pe > d->busy_end) d->busy_end = pe;
        }
        /* strong pulses in the preceding 1 ms */
        while (d->qdens < d->ns && lead - d->sleads[d->qdens % QN] > 1000.0 * d->spus) d->qdens++;
        if (d->qdens + QN <= d->ns) d->qdens = d->ns - QN + 1;
        np.density = (int)(d->ns - d->qdens);
        d->sleads[d->ns++ % QN] = lead;
    }
    /* the two segmenters finish pulses slightly out of order: keep the queue sorted */
    uint64_t k = d->qhead;
    while (k > d->qtail && d->q[(k - 1) & QMASK].lead > lead) { d->q[k & QMASK] = d->q[(k - 1) & QMASK]; k--; }
    d->q[k & QMASK] = np;
#ifdef DEBUG_PULSES
    fprintf(stderr, "P %.3f w=%.2f a=%.1f len=%d %s\n", lead / d->spus, np.width, amp, len, strong ? "S" : "w");
#endif
    d->qhead++;
    if (strong) d->pulses++;
    if (d->qhead - d->qtail >= QN - 64) assemble(d, true); /* overflow guard */
    assemble(d, false);
}

/* one sample into a segmenter; rel > 0 also ends the segment below rel x its peak */
static inline void segment(struct demod *d, struct seg *g, float a, uint64_t idx, float thr, float rel, bool strong)
{
    if (a >= thr && (!g->on || a >= rel * g->peak)) {
        if (!g->on) { g->on = 1; g->ps = idx; g->peak = a; }
        if (a > g->peak) g->peak = a;
        g->plast = idx;
    } else if (g->on && idx - g->plast > (uint64_t)d->merge_gap) {
        g->on = 0;
        finish_pulse(d, g, strong);
    } else if (g->on && idx - g->ps > (uint64_t)d->max_len) {
        g->on = 0;     /* CW / interference: drop */
    }
}

void demod_process(struct demod *d, const int8_t *iq, size_t ns)
{
    /* DC estimate from block mean, smoothed */
    double si = 0, sq = 0;
    for (size_t k = 0; k < ns; k++) { si += iq[2 * k]; sq += iq[2 * k + 1]; }
    float complex bm = (float)(si / ns) + I * (float)(sq / ns);
    d->dc = (d->noise < 0) ? bm : 0.9f * d->dc + 0.1f * bm;

    if (d->noise < 0) {
        double acc = 0;
        for (size_t k = 0; k < ns; k++) {
            float complex v = (iq[2 * k] + I * iq[2 * k + 1]) - d->dc;
            acc += cabsf(v);
        }
        d->noise = (float)(acc / ns);
    }
    d->thr = d->thr_abs > 0 ? d->thr_abs : d->noise * powf(10.0f, d->thr_db / 20.0f);
    d->thr_lo = d->thr * powf(10.0f, (d->modes_db - d->thr_db) / 20.0f);
    d->thr_sls = d->noise * 2.0f;
    float noise_cut = d->noise * 2.5f;

    for (size_t k = 0; k < ns; k++) {
        int8_t ii = iq[2 * k], qq = iq[2 * k + 1];
        if (ii == 127 || ii == -128 || qq == 127 || qq == -128) d->clipped++;
        float complex v = (ii + I * qq) - d->dc;
        float a = cabsf(v);
        uint64_t idx = d->n;
        d->z[idx & RMASK] = v;
        d->m[idx & RMASK] = a;
        d->n++;
        if (a > d->peak) d->peak = a;
        if (a < noise_cut) { d->nacc += a; d->ncnt++; }

        segment(d, &d->hi, a, idx, d->thr, 0, true);
        if (d->thr_lo < d->thr) segment(d, &d->lo, a, idx, d->thr_lo, SEG_REL, false);
    }
    d->rx += ns;
    if (d->ncnt > 1000) {
        /* mean of a Rayleigh distribution truncated at 2.5x mean is ~0.93 x mean */
        float nb = (float)(d->nacc / d->ncnt) / 0.93f;
        d->noise = 0.95f * d->noise + 0.05f * nb;
        d->nacc = 0; d->ncnt = 0;
    }
    assemble(d, false);
}

void demod_flush(struct demod *d)
{
    if (d->hi.on) { d->hi.on = 0; finish_pulse(d, &d->hi, true); }
    if (d->lo.on) { d->lo.on = 0; if (d->lo.peak < d->thr) finish_pulse(d, &d->lo, false); }
    assemble(d, true);
}

void demod_set_busy_gap(struct demod *d, double us) { d->busy_gap = (uint64_t)lround(us * d->spus); }

void demod_get_stats(struct demod *d, struct demod_stats *s, bool reset_peak)
{
    s->samples = d->n;
    s->clipped = d->clipped;
    s->rx = d->rx;
    s->busy = d->busy;
    s->pulses = d->pulses;
    s->stray_pulses = d->stray;
    s->garbled = d->garbled;
    s->noise_dbfs = dbfs(d->noise);
    s->thr_dbfs = dbfs(d->thr);
    s->peak_dbfs = dbfs(d->peak);
    if (reset_peak) d->peak = 0;
}

/* Advance the time base over samples that were not received on 1030 MHz
 * (ADS-B windows, retune settling). Pending pulses are finalised first. */
void demod_skip(struct demod *d, size_t ns)
{
    demod_flush(d);
    size_t z = ns < RING ? ns : RING;
    for (size_t k = 0; k < z; k++) {
        d->z[(d->n + k) & RMASK] = 0;
        d->m[(d->n + k) & RMASK] = 0;
    }
    d->n += ns;
}
