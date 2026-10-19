/* Pulse detection, interrogation classification and Mode S DPSK demodulation.
 *
 * Signal chain per sample: DC removal -> magnitude -> threshold segmentation
 * (with short-gap merging so DPSK phase-reversal dips do not split P6) ->
 * half-amplitude edge refinement -> pulse queue -> pattern matcher.
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

#define RING_BITS   20
#define RING        (1u << RING_BITS)
#define RMASK       (RING - 1)
#define QN          4096
#define QMASK       (QN - 1)
#define MAX_PULSE_US 40.0

struct demod {
    double fs, spus;            /* samples per second / per microsecond */
    float complex *z;           /* ring of DC-corrected samples */
    float *m;                   /* ring of magnitudes */
    uint64_t n;                 /* global index of next sample */

    float complex dc;
    float noise;                /* mean noise magnitude */
    float thr_db, thr_abs, thr;
    double nacc; uint64_t ncnt;

    int in_pulse;
    uint64_t ps, plast;
    int merge_gap, min_len, max_len;

    struct pulse q[QN];
    uint64_t qhead, qtail;      /* qhead: next write; qtail: oldest unexamined */
    uint64_t qdens;             /* oldest pulse inside the 1 ms density window */
    int dens_max;
    uint64_t garbled;

    interrog_cb cb;
    void *ctx;
    uint64_t seq;

    uint64_t clipped, pulses, stray;
    float peak;
    float *tmp;
};

struct demod *demod_new(double fs, float thr_db, float thr_abs, int dens_max, interrog_cb cb, void *ctx)
{
    struct demod *d = calloc(1, sizeof(*d));
    d->fs = fs;
    d->spus = fs / 1e6;
    d->z = calloc(RING, sizeof(float complex));
    d->m = calloc(RING, sizeof(float));
    d->thr_db = thr_db;
    d->thr_abs = thr_abs;
    d->dens_max = dens_max;
    d->noise = -1;
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

static void decode_modes(struct demod *d, struct interrog *it, double p6lead, double p6width)
{
    double chip = 0.25 * d->spus;
    int h = (int)(chip / 2);
    if (h < 1) h = 1;

    /* locate sync phase reversal near P6 + 1.25 us */
    double exp_sync = p6lead + 1.25 * d->spus;
    double best_t = exp_sync;
    float best = 2;
    int span = (int)lround(0.125 * d->spus);
    for (int o = -span; o <= span; o++) {
        float c = reversal_corr(d, exp_sync + o, h, NULL);
        if (c < best) { best = c; best_t = exp_sync + o; }
    }
    it->sync = best;

    /* P5 (Mode S SLS) overlaps the sync reversal. P6 flips phase there, P5
     * does not: with z1/z2 the mean samples just before/after the reversal,
     * P5 = |z1 + z2| / 2 and P6 = |z1 - z2| / 2. Uses the nominal sync when
     * a strong P5 masks the reversal. */
    {
        double ts = best > -0.5f ? exp_sync : best_t;
        int a0 = (int)lround(0.125 * d->spus), a1 = (int)lround(0.375 * d->spus);
        if (a1 <= a0) a1 = a0 + 1;
        float complex z1 = 0, z2 = 0;
        uint64_t c = (uint64_t)llround(ts);
        for (int k = a0; k < a1; k++) { z1 += Z(d, c - k); z2 += Z(d, c + k); }
        z1 /= (a1 - a0); z2 /= (a1 - a0);
        float p5 = cabsf(z1 + z2) / 2;
        it->p5_db = p5 >= d->thr ? dbfs(p5) : NAN;     /* absolute here, made relative to P1 by the caller */
    }
    if (best > -0.5f) {
        it->kind = K_S_NOSYNC;
        it->sidelobe = true;
        return;
    }

    int len_class = 0;
    if (p6width > 14.5 && p6width < 18.0) len_class = 56;
    else if (p6width > 28.5 && p6width < 32.5) len_class = 112;

    /* Data reversals occur N x 0.25 us after sync with N >= 2 (Annex 10
     * Vol IV 3.1.2.1.5.2.3), so data bit k sits at sync + (k + 1) chips. */
    int maxbits = (p6width > 18.0) ? 112 : 56;
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
    case 11: it->kind = K_UF11; break;
    case 16: it->kind = K_UF16; it->tcas = true; break;
    case 19: it->kind = K_UF19; break;
    case 20: it->kind = K_UF20; break;
    case 21: it->kind = K_UF21; break;
    case 24: it->kind = K_UF24; break;
    default: it->kind = K_UF_OTHER; break;
    }
}

/* ---------------- pattern matcher ---------------- */

static inline double rel_us(struct demod *d, const struct pulse *a, const struct pulse *b)
{
    return (b->lead - a->lead) / d->spus;
}

static bool is_short(const struct pulse *p) { return p->width > 0.6 && p->width < 1.05; }

/* find unused pulse after index i at offset off +- tol us */
static int find_at(struct demod *d, uint64_t i, double off, double tol, bool need_short)
{
    struct pulse *p1 = &d->q[i & QMASK];
    for (uint64_t j = i + 1; j < d->qhead; j++) {
        struct pulse *p = &d->q[j & QMASK];
        double r = rel_us(d, p1, p);
        if (r > off + tol) break;
        if (p->used || r < off - tol) continue;
        if (need_short && !is_short(p)) continue;
        return (int)(j - i);
    }
    return -1;
}

static void emit(struct demod *d, struct interrog *it, const struct pulse *p1)
{
    it->seq = d->seq++;
    it->t = p1->lead / d->fs;
    it->amp_db = dbfs(p1->amp);
    it->burst = -1;
    d->cb(it, d->ctx);
}

static void match(struct demod *d, uint64_t i)
{
    struct pulse *p1 = &d->q[i & QMASK];
    if (p1->used) return;
    if (!is_short(p1)) { d->stray++; return; }

    struct interrog it;
    memset(&it, 0, sizeof(it));
    it.p2_db = NAN;
    it.p5_db = NAN;

    /* --- Mode S: P1, P2 @2.0, P6 @3.5 --- */
    int j2 = find_at(d, i, 2.0, 0.2, true);
    if (j2 > 0) {
        int j6 = find_at(d, i, 3.5, 0.25, false);
        if (j6 > 0) {
            struct pulse *p6 = &d->q[(i + j6) & QMASK];
            /* chain fragments of a broken-up P6 */
            double end = p6->lead + p6->width * d->spus;
            uint64_t last = i + j6;
            for (uint64_t k = i + j6 + 1; k < d->qhead; k++) {
                struct pulse *p = &d->q[k & QMASK];
                if (p->lead - end > 0.5 * d->spus) break;
                double e = p->lead + p->width * d->spus;
                if (e > end) end = e;
                last = k;
            }
            double w6 = (end - p6->lead) / d->spus;
            if (w6 > 10.0) {
                struct pulse *p2 = &d->q[(i + j2) & QMASK];
                it.p2_db = 20.0f * log10f(p2->amp / p1->amp);
                decode_modes(d, &it, p6->lead, w6);
                if (!isnan(it.p5_db)) it.p5_db -= dbfs(p1->amp);
                p1->used = p2->used = true;
                for (uint64_t k = i + j6; k <= last; k++) d->q[k & QMASK].used = true;
                emit(d, &it, p1);
                return;
            }
        }
    }

    /* --- Mode A/C with optional P2, P4, S1. Modes 1/2/B/D are not detected:
     * on air they were almost all coincidences of unrelated pulses. --- */
    static const double off3[] = { 8, 21 };
    static const int base[] = { K_MODEA, K_MODEC };
    for (int m = 0; m < 2; m++) {
        int j3 = find_at(d, i, off3[m], 0.15, true);
        if (j3 < 0) continue;
        struct pulse *p3 = &d->q[(i + j3) & QMASK];
        if (fabsf(20.0f * log10f(p3->amp / p1->amp)) > 3.0f) continue;
        if (p1->density > d->dens_max) { d->garbled++; return; }
        it.kind = base[m];
        if (j2 > 0 && j2 < j3) {
            struct pulse *p2 = &d->q[(i + j2) & QMASK];
            it.p2_db = 20.0f * log10f(p2->amp / p1->amp);
            if (it.p2_db >= 0) it.sidelobe = true;
            p2->used = true;
        }
        int j4 = find_at(d, i, off3[m] + 2.0, 0.15, false);
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
            if (r > 1.8 && is_short(p)) { it.s1 = true; break; }
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

static void finish_pulse(struct demod *d)
{
    uint64_t s = d->ps, e = d->plast;
    int len = (int)(e - s + 1);
    if (len < d->min_len || len > d->max_len) return;

    /* plateau amplitude: median of the segment, ignoring the edges */
    int a0 = len > 6 ? 2 : 0, n = 0;
    for (int k = a0; k < len - a0; k++) d->tmp[n++] = M(d, s + k);
    qsort(d->tmp, n, sizeof(float), cmpf);
    float amp = d->tmp[n / 2];
    if (len > 24) amp = d->tmp[(n * 3) / 4]; /* long DPSK pulse: skip reversal dips */
    if (amp < 0.8f * d->thr) return;   /* merged noise blips, not a pulse */
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

    struct pulse *p = &d->q[d->qhead & QMASK];
    p->lead = lead;
    p->width = (trail - lead) / d->spus;
    p->amp = amp;
    p->used = false;
    while (d->qdens < d->qhead && lead - d->q[d->qdens & QMASK].lead > 1000.0 * d->spus) d->qdens++;
    if (d->qdens + QN <= d->qhead) d->qdens = d->qhead - QN + 1;
    p->density = (int)(d->qhead - d->qdens);
#ifdef DEBUG_PULSES
    fprintf(stderr, "P %.3f w=%.2f a=%.1f len=%d\n", lead / d->spus, p->width, amp, len);
#endif
    d->qhead++;
    d->pulses++;
    if (d->qhead - d->qtail >= QN - 64) assemble(d, true); /* overflow guard */
    assemble(d, false);
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

        if (a >= d->thr) {
            if (!d->in_pulse) { d->in_pulse = 1; d->ps = idx; }
            d->plast = idx;
        } else if (d->in_pulse && idx - d->plast > (uint64_t)d->merge_gap) {
            d->in_pulse = 0;
            finish_pulse(d);
        } else if (d->in_pulse && idx - d->ps > (uint64_t)d->max_len) {
            d->in_pulse = 0;     /* CW / interference: drop */
        }
    }
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
    if (d->in_pulse) { d->in_pulse = 0; finish_pulse(d); }
    assemble(d, true);
}

void demod_get_stats(struct demod *d, struct demod_stats *s, bool reset_peak)
{
    s->samples = d->n;
    s->clipped = d->clipped;
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
