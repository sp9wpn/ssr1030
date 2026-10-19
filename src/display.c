/* Live ANSI terminal dashboard. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include "analysis.h"

#define B   "\033[1m"
#define R   "\033[0m"
#define RED "\033[31m"
#define YEL "\033[33m"
#define GRN "\033[32m"
#define DIM "\033[2m"
#define EL  "\033[K\n"

static uint64_t last_samples, last_clip, last_pulses, last_clip_prev, last_samples_prev;
static long last_total[K_COUNT];
static double last_now = -1;
static float rate[K_COUNT];

/* signal level history, so warnings cover several antenna rotations */
#define LVL_N 128
static struct { double t; float peak; uint64_t clip, samples; } lvl[LVL_N];
static int nlvl;
static double last_clip_t = -1e9;
static int last_gain = -1;

void display_init(void)
{
    fputs("\033[?1049h\033[?25l\033[2J", stdout);
    fflush(stdout);
}

void display_restore(void)
{
    fputs("\033[?25h\033[?1049l", stdout);
    fflush(stdout);
}

static int term_cols(void)
{
    struct winsize w;
    if (ioctl(1, TIOCGWINSZ, &w) == 0 && w.ws_col > 40) return w.ws_col;
    return 120;
}

static void print_seq(const struct burst *b, int cols, int maxlines)
{
    int col = 2, lines = 0, rank[64];
    burst_ranks(b, rank);
    fputs("  ", stdout);
    for (int s = 0; s < b->nks;) {
        char tok[48];
        s = burst_item(b, s, rank, tok, sizeof(tok));
        if (col + (int)strlen(tok) >= cols - 1) {
            if (++lines >= maxlines) { fputs("...", stdout); break; }
            fputs(EL "  ", stdout);
            col = 2;
        }
        fputs(tok, stdout);
        col += strlen(tok);
    }
    fputs(EL, stdout);
    char tl[1024];
    burst_targets(b, rank, tl, sizeof(tl));
    if (b->naddr) {
        /* wrap the target list like the sequence */
        fputs("  targets:", stdout);
        col = 10;
        lines = 0;
        for (char *t = strtok(tl, " "); t; t = strtok(NULL, " ")) {
            if (col + (int)strlen(t) + 1 >= cols - 1) {
                if (++lines >= 2) { fputs(" ...", stdout); break; }
                fputs(EL "          ", stdout);
                col = 10;
            }
            printf(" %s", t);
            col += strlen(t) + 1;
        }
    }
    fputs(EL, stdout);
}

void display_render(const struct run_info *ri, double now, struct demod *d, uint64_t drops)
{
    struct demod_stats st;
    demod_get_stats(d, &st, true);
    int cols = term_cols();

    double dt = (last_now >= 0) ? now - last_now : 0;
    double clip_pct = 0, pps = 0;
    if (dt > 0) {
        uint64_t ds = st.samples - last_samples;
        clip_pct = ds ? 100.0 * (st.clipped - last_clip) / ds : 0;
        pps = (st.pulses - last_pulses) / dt;
        for (int k = 0; k < K_COUNT; k++) {
            float r = (AS.total[k] - last_total[k]) / dt;
            rate[k] = rate[k] * 0.7f + r * 0.3f;
        }
    }
    last_now = now; last_samples = st.samples; last_clip = st.clipped; last_pulses = st.pulses;
    memcpy(last_total, AS.total, sizeof(last_total));

    fputs("\033[H", stdout);
    int h = (int)now / 3600, m = ((int)now / 60) % 60, s = (int)now % 60;
    printf(B "ssr1030" R "  %s  %.3f MHz  %.1f MS/s", ri->source, ri->freq / 1e6, ri->fs / 1e6);
    printf("  t=%02d:%02d:%02d" EL, h, m, s);

    /* level hold window: max(20 s, 3 rotations), cleared when the gain changes */
    int gain_sig = ri->amp * 10000 + ri->lna * 100 + ri->vga;
    if (gain_sig != last_gain) { nlvl = 0; last_gain = gain_sig; last_clip_t = -1e9; }
    if (dt > 0) {
        lvl[nlvl % LVL_N].t = now;
        lvl[nlvl % LVL_N].peak = st.peak_dbfs;
        lvl[nlvl % LVL_N].clip = st.clipped - last_clip_prev;
        lvl[nlvl % LVL_N].samples = st.samples - last_samples_prev;
        nlvl++;
        if (clip_pct > 0.0005) last_clip_t = now;
    }
    last_clip_prev = st.clipped; last_samples_prev = st.samples;
    double rot = analysis_rotation();
    double win = rot * 3 > 20 ? rot * 3 : 20;
    float peak_w = -99; uint64_t clip_w = 0, samp_w = 0;
    double span = 0;
    for (int i = 0; i < nlvl && i < LVL_N; i++) {
        int k = (nlvl - 1 - i) % LVL_N;
        if (now - lvl[k].t > win) break;
        if (lvl[k].peak > peak_w) peak_w = lvl[k].peak;
        clip_w += lvl[k].clip; samp_w += lvl[k].samples;
        span = now - lvl[k].t;
    }
    double clip_w_pct = samp_w ? 100.0 * clip_w / samp_w : 0;
    bool clipping = now - last_clip_t < win;
    bool weak = !clipping && span > 0.8 * win && peak_w - st.noise_dbfs < 25;

    if (ri->live) {
        printf("HackRF: %s  " B "AMP %s  LNA %d dB  VGA %d dB" R "  (nominal %d dB)  %s%s" R,
               ri->hw, ri->amp ? "on" : "off", ri->lna, ri->vga, (ri->amp ? 11 : 0) + ri->lna + ri->vga,
               ri->agc ? GRN "AUTO GAIN" : YEL "MANUAL GAIN ('g' = auto)", "");
        if (ri->gain_age >= 0)
            printf(DIM "  last change %+d dB %.0f s ago (%s)" R, ri->gain_delta, ri->gain_age, ri->gain_why);
        fputs(EL, stdout);
    }
    printf("Signal: noise %.1f dBFS  thr %.1f  peak %s%.1f" R " (last %.0f s: %s%.1f" R ")  clip %s%.4f%%" R
           " (last %.0f s: %.4f%%)  pulses %.0f/s  interf-rejected %lu",
           st.noise_dbfs, st.thr_dbfs, clip_pct > 0.0005 ? RED : "", st.peak_dbfs, win,
           clipping ? RED : "", peak_w, clip_pct > 0.0005 ? RED : "", clip_pct, win, clip_w_pct, pps,
           (unsigned long)st.garbled);
    if (drops) printf("  " RED "USB drops %lu" R, (unsigned long)drops);
    fputs(EL, stdout);
    if (ri->agc && ri->live)
        printf("  %s" R DIM "  auto gain keeps the main radar's beam peak near -6 dBFS (peak %.1f, clip %.4f%% over %.0f s)" R EL,
               clipping ? RED "clipping seen" : (weak ? YEL "weak" : GRN "level OK"), peak_w, clip_w_pct, win);
    else if (clipping)
        printf(RED B "  CLIPPING %.0f s ago" R RED " - lower gain: '-' VGA -2 dB, 'l' LNA -8 dB, 'a' amp off"
               "  (warning holds %.0f s)" R EL, now - last_clip_t, win);
    else if (weak)
        printf(YEL B "  WEAK SIGNAL" R YEL " - strongest pulse in %.0f s only %.1f dB above noise; raise gain: "
               "'+' VGA +2 dB, 'L' LNA +8 dB, 'a' amp on" R EL, win, peak_w - st.noise_dbfs);
    else if (span < 0.8 * win)
        printf(DIM "  level check: collecting %.0f/%.0f s since start or last gain change" R EL, span, win);
    else printf(GRN "  level OK" R DIM " - peak %.1f dB above noise over %.0f s, no clipping" R EL,
                peak_w - st.noise_dbfs, win);
    {
        char fl[200];
        int o = snprintf(fl, sizeof(fl), "Filters: %s", ri->filters[0] ? ri->filters : "none");
        if (AS.filt_radar || AS.filt_modes || AS.filt_rollcall || AS.filt_aircraft)
            snprintf(fl + o, sizeof(fl) - o, "  (dropped: radar %ld, non-S %ld, all-call %ld, aircraft %ld)",
                     AS.filt_radar, AS.filt_modes, AS.filt_rollcall, AS.filt_aircraft);
        printf("%s", fl);
        if (ri->adsb) {
            struct adsb_stats as;
            adsb_get_stats(&as);
            printf("   ADS-B: %s%s" R " %d callsigns, %ld windows, %ld frames (%ld id)",
                   ri->adsb_now ? GRN : "", ri->adsb_now ? "LISTENING 1090" :
                   (ri->adsb_next < 0 ? "waiting for gap after burst" : "idle"),
                   cs_count(), as.windows, as.frames, as.idents + as.commb_idents);
            if (!ri->adsb_now && ri->adsb_next >= 0) printf(", next in %.0f s", ri->adsb_next);
        }
        fputs(EL, stdout);
    }
    printf("Background (sidelobe) %.1f dBFS   beam member >= %.1f   core >= %.1f   TCAS/ACAS (excluded): "
           GRN "%ld" R " (UF0 %ld, UF16 %ld, WS %ld; %d aircraft)" EL,
           AS.bg_db, AS.bg_db + AS.cfg.member_db, AS.bg_db + AS.cfg.core_db,
           AS.tcas_total, AS.tcas_uf0, AS.tcas_uf16, AS.tcas_ws, AS.naddr_tcas);
    fputs(EL, stdout);

    /* left: per kind table; right: radars */
    char left[24][80], right[24][160];
    int nl = 0, nr = 0;
    snprintf(left[nl++], 80, B "%-24s %8s %7s %7s %6s" R, "Interrogation", "total", "in-beam", "sidelb", "/s");
    for (int k = 0; k < K_COUNT && nl < 24; k++) {
        if (!AS.total[k] || k == K_UF0 || k == K_UF16) continue;
        snprintf(left[nl++], 80, "%-24s %8ld %7ld %7ld %6.1f", kind_label(k), AS.total[k], AS.inbeam[k],
                 AS.side[k], rate[k]);
    }
    if (nl == 1) snprintf(left[nl++], 80, DIM "(no radar interrogations yet)" R);

    snprintf(right[nr++], 160, B "Radar   bursts  rotation  dwell  interrog/burst avg [min-max]" R);
    for (int r = 0; r < AS.nradars && nr < 22; r++) {
        struct radar *X = &AS.radars[r];
        int np = X->np < 64 ? X->np : 64;
        char rot[16] = "      -";
        if (np) snprintf(rot, sizeof(rot), "%7.3fs", median_d(X->periods, np));
        snprintf(right[nr++], 160, "%-6s %7d  %s %5.1fms  total %.1f", X->key, X->nbursts + X->npartial,
                 rot, X->nbursts ? 1e3 * X->dwell_sum / X->nbursts : 0,
                 X->nbursts ? (double)X->total_sum / X->nbursts : 0);
        char line[160]; int o = 0;
        o += snprintf(line + o, sizeof(line) - o, "      ");
        for (int g = 0; g < G_COUNT && o < 150; g++) {
            if (!X->cnt_max[g]) continue;
            o += snprintf(line + o, sizeof(line) - o, " %s %.1f[%d-%d]", group_label(g),
                          (double)X->cnt_sum[g] / X->nbursts, X->cnt_min[g], X->cnt_max[g]);
        }
        snprintf(right[nr++], 160, "%s", line);
    }
    if (AS.nradars == 0) snprintf(right[nr++], 160, DIM "(waiting for a main-beam burst)" R);
    nr++; right[nr - 1][0] = 0;
    {
        char line[160]; int o = 0;
        o += snprintf(line, sizeof(line), B "UF11 codes" R);
        for (int si = 0; si < 2; si++)
            for (int c = 0; c < 64; c++)
                for (int p = 0; p < 16; p++)
                    if (AS.uf11_code[si][c][p] && o < 150)
                        o += snprintf(line + o, sizeof(line) - o, " %s%d/PR%d:%ld", si ? "SI" : "II", c, p,
                                      AS.uf11_code[si][c][p]);
        snprintf(right[nr++], 160, "%s", line);
        o = snprintf(line, sizeof(line), B "Comm-B BDS" R);
        for (int b = 0; b < 256; b++)
            if (AS.bds[b] && o < 150) o += snprintf(line + o, sizeof(line) - o, " %X,%X:%ld", b >> 4, b & 15, AS.bds[b]);
        snprintf(right[nr++], 160, "%s", line);
        o = snprintf(line, sizeof(line), B "DI" R);
        for (int k = 0; k < 8; k++) if (AS.di[k]) o += snprintf(line + o, sizeof(line) - o, " %d:%ld", k, AS.di[k]);
        o += snprintf(line + o, sizeof(line) - o, "  " B "PC" R);
        for (int k = 0; k < 8; k++) if (AS.pc[k]) o += snprintf(line + o, sizeof(line) - o, " %d:%ld", k, AS.pc[k]);
        snprintf(right[nr++], 160, "%s", line);
        snprintf(right[nr++], 160, B "Roll-called aircraft" R " %d", AS.naddr_radar);
    }
    int rows = nl > nr ? nl : nr;
    for (int i = 0; i < rows; i++) {
        /* left column is padded by visible width (escape codes excluded) */
        const char *l = i < nl ? left[i] : "";
        int vis = 0;
        for (const char *p = l; *p; p++) {
            if (*p == '\033') { while (*p && *p != 'm') p++; continue; }
            vis++;
        }
        printf("%s%*s  %s" EL, l, 58 - vis > 0 ? 58 - vis : 0, "", i < nr ? right[i] : "");
    }
    fputs(EL, stdout);

    /* recent bursts */
    printf(B "%-5s %9s %-5s %7s %6s %4s %5s %7s", "Burst", "t[s]", "radar", "dur ms", "peak", "gain", "n", "P[s]");
    for (int g = 0; g < G_COUNT; g++) printf(" %6s", group_label(g));
    printf(R EL);
    int nb = AS.nhist < 8 ? AS.nhist : 8;
    for (int i = AS.nhist - nb; i < AS.nhist; i++) {
        struct burst *b = &AS.hist[i % BURST_HIST];
        printf("#%-4d %9.3f %-4s%s %7.1f %6.1f %4.0f %5d", b->id, b->t0, b->key, b->partial ? "*" : " ",
               1e3 * (b->t1 - b->t0), b->peak, b->gain, b->n);
        if (b->period > 0) printf(" %7.3f", b->period); else printf(" %7s", "-");
        for (int g = 0; g < G_COUNT; g++) {
            if (b->gcnt[g]) printf(" %6d", b->gcnt[g]);
            else printf(" " DIM "%6s" R, ".");
        }
        fputs(EL, stdout);
    }
    for (int i = nb; i < 8; i++) fputs(EL, stdout);

    if (AS.nhist > 0) {
        struct burst *b = &AS.hist[(AS.nhist - 1) % BURST_HIST];
        printf(B "Last burst #%d sequence" R DIM
               "  (A C A/C-only all-call, AA CC plain A/C, S UF11 all-call, alt/sqk UF4/5, acid/setAlt/trk/hdgSpd BDS 2,0/4,0/5,0/6,0, S-msg Comm-A, #N repeat)" R EL,
               b->id);
        print_seq(b, cols, 4);
    } else fputs(EL EL EL, stdout);

    /* transitions */
    printf(B "In-beam interleaving (most frequent transitions, dt avg [min-max] ms)" R EL);
    struct trans top[8];
    int nt = analysis_top_transitions(top, 8);
    for (int i = 0; i < 8; i += 2) {
        for (int j = i; j < i + 2; j++) {
            if (j >= nt) continue;
            /* fixed-width cell so the second column lines up */
            char ta[24], tb[24], cell[96];
            snprintf(cell, sizeof(cell), "%-12s -> %-12s %6ld %8.3f [%7.3f -%8.3f]",
                     seq_token_str(top[j].a, ta, sizeof(ta)), seq_token_str(top[j].b, tb, sizeof(tb)), top[j].n,
                     1e3 * top[j].sum / top[j].n, 1e3 * top[j].min, 1e3 * top[j].max);
            printf("  %-66s", cell);
        }
        fputs(EL, stdout);
    }

    printf(B "Recent radar interrogations" R EL);
    int nrc = AS.nrecent < RECENT_N ? AS.nrecent : RECENT_N;
    for (int i = AS.nrecent - nrc; i < AS.nrecent; i++) {
        char buf[300];
        FILE *mf = fmemopen(buf, sizeof(buf), "w");
        print_interrog_line(mf, &AS.recent[i % RECENT_N]);
        fclose(mf);
        buf[strcspn(buf, "\n")] = 0;
        if ((int)strlen(buf) > cols - 1) buf[cols - 1] = 0;
        printf("%s" EL, buf);
    }
    printf(DIM "keys: q quit  +/- VGA  L/l LNA  a amp (manual)  g auto gain  r reset stats" R EL);
    fputs("\033[J", stdout);
    fflush(stdout);
}
