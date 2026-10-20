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
/* snapshots of the counters for per-minute rates over the last 60 s */
#define SNAP_N 160
static struct { double t; long total[K_COUNT]; long bds[256]; } snap[SNAP_N];
static int nsnap;
static double last_now = -1;
static float rate[K_COUNT], bds_rate[256];   /* per minute */

/* spaces needed to pad s to w visible columns (escape sequences excluded) */
static int pad_to(const char *s, int w)
{
    int vis = 0;
    for (const char *p = s; *p; p++) {
        if (*p == '\033') { while (*p && *p != 'm') p++; continue; }
        if ((*p & 0xC0) != 0x80) vis++;
    }
    return w - vis > 0 ? w - vis : 0;
}

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

static int term_rows(void)
{
    struct winsize w;
    if (ioctl(1, TIOCGWINSZ, &w) == 0 && w.ws_row > 5) return w.ws_row;
    return 50;
}

/* burst table column headers (the UF11 column gets a longer name than its token) */
static const char *tbl_label(int g) { return g == G_SALL ? "UF11-all" : group_label(g); }
static int tbl_w(int g)
{
    int w = (int)strlen(tbl_label(g));
    return w > 6 ? w : 6;
}

static void draw_status_line(const struct run_info *ri, bool last_row);

/* print s cut to `width` visible columns (escape sequences and UTF-8
 * continuation bytes do not count); resets attributes at the end */
static void put_cut(const char *s, int width)
{
    int vis = 0;
    for (const char *p = s; *p; p++) {
        if (*p == '\033') {
            const char *q = p;
            while (*q && *q != 'm' && *q != 'K' && *q != 'H' && *q != 'J') q++;
            fwrite(p, 1, (size_t)(q - p) + (*q ? 1 : 0), stdout);
            if (!*q) break;
            p = q;
            continue;
        }
        if ((*p & 0xC0) != 0x80) {
            if (vis >= width) break;
            vis++;
        }
        fputc(*p, stdout);
    }
    fputs(R, stdout);
}

void display_status_only(const struct run_info *ri)
{
    printf("\033[%d;1H", term_rows());
    draw_status_line(ri, true);
    fflush(stdout);
}

/* bottom line: filter input, last command result, or key help. On the last
 * screen row it must not end with a newline, or the terminal scrolls. */
static void draw_status_line(const struct run_info *ri, bool last_row)
{
    char line[512];
    int o = 0;
    if (ri->frozen) o += snprintf(line + o, sizeof(line) - o, "\033[1;7m FROZEN \033[0m ");
    if (ri->cmd_mode)
        snprintf(line + o, sizeof(line) - o, B "%s: " R "%s\033[7m \033[0m" DIM "   (Enter apply, empty = no filter, Esc cancel, Ctrl-U clear)" R,
                 ri->cmd_kind == 'r' ? "Radar filter (SI1-63 or II0-15)" : "Aircraft filter (hex or callsign)", ri->cmd);
    else if (ri->cmd_msg[0])
        snprintf(line + o, sizeof(line) - o, YEL "%s" R DIM "   keys: q  +/-  L/l  p  g  r radar  a aircraft  S Mode S only  R roll-calls only  b graph  space freeze  x reset" R,
                 ri->cmd_msg);
    else if (ri->graph)
        snprintf(line + o, sizeof(line) - o, DIM "keys: b dashboard  \u2191\u2193 PgUp PgDn Home End scroll  r radar filter  a aircraft filter  "
                 "S Mode S only  R roll-calls only  space freeze  q quit" R);
    else snprintf(line + o, sizeof(line) - o, DIM "keys: q quit  +/- VGA  L/l LNA (manual)  p RF amp  g auto gain  r radar filter  a aircraft filter  "
                  "S Mode S only  R roll-calls only  b burst graph  space freeze  x reset stats" R);
    /* on the last row one column stays free so the terminal never wraps/scrolls */
    put_cut(line, term_cols() - 1);
    fputs(last_row ? "\033[K" : EL, stdout);
}

/* ---------- burst graph screen ('b') ----------
 * One row per query, in time order: time relative to the middle of the beam
 * pass, query/target, amplitude and a horizontal bar whose length is the
 * amplitude. The bar character tells the class: beam, sidelobe, weak. */
static const char *gp_char(int cls) { return cls == GP_BEAM ? "\u2588" : cls == GP_SIDE ? "\u2591" : "\u2592"; }
static const char *gp_color(int cls) { return cls == GP_BEAM ? GRN : cls == GP_SIDE ? YEL : DIM; }

#define GRAPH_BOT_DBFS (-40.0f)
#define GRAPH_TOP_DBFS (3.0f)

static int g_scroll, g_burst = -1, g_page = 20, g_key;
static bool frozen_view;

void display_graph_key(int key)
{
    g_key = key;            /* applied at the next draw, where the page size is known */
}

static void draw_graph(int cols, int rows)
{
    /* frozen: keep showing the pinned burst while it is still in the history */
    static int pin = -1;
    const struct burst *b = NULL;
    if (frozen_view && pin >= 0) {
        int lim = AS.nhist < BURST_HIST ? AS.nhist : BURST_HIST;
        for (int i = 0; i < lim && !b; i++)
            if (AS.hist[(AS.nhist - 1 - i) % BURST_HIST].id == pin) b = &AS.hist[(AS.nhist - 1 - i) % BURST_HIST];
    }
    if (!b) b = analysis_graph_burst();
    pin = frozen_view && b ? b->id : -1;
    if (!b || !b->ngp) {
        printf(B "Burst graph" R DIM "  waiting for a main-beam burst   b back" R EL);
        return;
    }
    double mid = 0.5 * (b->t0 + b->t1);
    int n = b->ngp, first_beam = 0;
    float amax = -99, amin = 99;
    for (int i = n - 1; i >= 0; i--) {
        if (b->gp[i].amp > amax) amax = b->gp[i].amp;
        if (!isnan(b->gp[i].p2) && b->gp[i].p2 > amax) amax = b->gp[i].p2;
        if (!isnan(b->gp[i].p5) && b->gp[i].p5 > amax) amax = b->gp[i].p5;
        if (b->gp[i].amp < amin) amin = b->gp[i].amp;
        /* the scale must also cover the SLS pulses, which are often weaker than every bar */
        if (!isnan(b->gp[i].p2) && b->gp[i].p2 < amin) amin = b->gp[i].p2;
        if (!isnan(b->gp[i].p5) && b->gp[i].p5 < amin) amin = b->gp[i].p5;
        if (b->gp[i].cls == GP_BEAM) first_beam = i;
    }
    /* fixed scale: left edge just below the noise floor, right edge where
     * I/Q clipping starts, so bars compare between bursts */
    float ytop = GRAPH_TOP_DBFS, ybot = GRAPH_BOT_DBFS;
    (void)amax; (void)amin;

    int page = rows - 3;
    if (page < 5) page = 5;
    g_page = page;
    if (b->id != g_burst) {                 /* new burst: start just before the beam */
        g_burst = b->id;
        g_scroll = first_beam - page / 4;
    }
    switch (g_key) {
    case 'U': g_scroll--; break;
    case 'D': g_scroll++; break;
    case 'P': g_scroll -= page - 1; break;
    case 'N': g_scroll += page - 1; break;
    case 'H': g_scroll = 0; break;
    case 'E': g_scroll = n; break;
    }
    g_key = 0;
    if (g_scroll > n - page) g_scroll = n - page;
    if (g_scroll < 0) g_scroll = 0;
    int last = g_scroll + page < n ? g_scroll + page : n;

    {
        char hdr[256];
        snprintf(hdr, sizeof(hdr), B "Burst #%d %s%s" R "  dur %.1f ms  peak %.1f dBFS  gain %.0f dB  beam %d  sidelobe %d  weak %d"
                 DIM "  rows %d-%d/%d  bars %.0f..%+.0f dBFS" R,
                 b->id, b->key, b->partial ? "*" : "", 1e3 * (b->t1 - b->t0), b->peak, b->gain, b->n, b->nside,
                 b->nweak, g_scroll + 1, last, n, GRAPH_BOT_DBFS, GRAPH_TOP_DBFS);
        put_cut(hdr, cols - 1);
        fputs(EL, stdout);
    }

    /* bar width: screen minus the row prefix, 5 columns spare */
    char probe[96];
    int plen = snprintf(probe, sizeof(probe), "%+10.2f ms  %-12.12s %-8.8s %-3s %6.1f ", -9999.99, "", "", "", -99.9);
    int W = cols - plen - 5;
    if (W < 10) W = 10;
    for (int i = g_scroll; i < last; i++) {
        const struct gpoint *g = &b->gp[i];
        char ts[24];
        seq_token_str(g->tok, ts, sizeof(ts));
        const char *tg = g->addr != NO_ADDR ? ac_name(g->addr) : "";
        int len = (int)lroundf((g->amp - ybot) / (ytop - ybot) * W);
        if (len < 1) len = 1;
        if (len > W) len = W;
        /* SLS pulse marker "P2"/"P5" ending at the column of its level */
        const char *mk = NULL;
        float mv = NAN;
        if (!isnan(g->p2)) { mk = "P2"; mv = g->p2; }
        else if (!isnan(g->p5)) { mk = "P5"; mv = g->p5; }
        int mpos = -1;
        if (mk) {
            mpos = (int)lroundf((mv - ybot) / (ytop - ybot) * W);
            if (mpos < 2) mpos = 2;
            if (mpos > W) mpos = W;
        }
        /* fixed-width fields: query 12 (longest token, e.g. S-msg/hdgSpd), target 8 (callsign or hex),
         * "sls" when the query reached us through a sidelobe */
        printf("%s%+10.2f ms  %-12.12s %-8.8s %-3s %6.1f ", g->cls == GP_BEAM ? "" : DIM, 1e3 * (g->t - mid), ts, tg,
               g->cls == GP_SIDE ? "sls" : "", g->amp);
        fputs(R, stdout);
        int end = mpos > len ? mpos : len;
        for (int k = 0; k < end; k++) {
            if (mk && k >= mpos - 2 && k < mpos) {
                /* inside the bar: inverse, beyond it: plain bold */
                printf(k < len ? R "\033[1;7;31m%c" R : R B RED "%c" R, mk[k - (mpos - 2)]);
            } else if (k < len) {
                printf("%s%s", gp_color(g->cls), gp_char(g->cls));
            } else {
                fputs(R " ", stdout);
            }
        }
        fputs(R EL, stdout);
    }
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
    }
    last_now = now; last_samples = st.samples; last_clip = st.clipped; last_pulses = st.pulses;

    /* per-minute rates: counts over the last 60 s (or since start / reset) */
    {
        int cur = nsnap % SNAP_N;
        snap[cur].t = now;
        memcpy(snap[cur].total, AS.total, sizeof(AS.total));
        memcpy(snap[cur].bds, AS.bds, sizeof(AS.bds));
        nsnap++;
        int old = cur;
        for (int i = 1; i < nsnap && i < SNAP_N; i++) {
            int j = (cur - i + SNAP_N) % SNAP_N;
            if (now - snap[j].t > 60) break;
            bool reset = false;              /* counters went down: statistics were reset */
            for (int k = 0; k < K_COUNT && !reset; k++) if (snap[j].total[k] > AS.total[k]) reset = true;
            if (reset) break;
            old = j;
        }
        double span = now - snap[old].t;
        for (int k = 0; k < K_COUNT; k++)
            rate[k] = span > 0 ? (float)((AS.total[k] - snap[old].total[k]) * 60.0 / span) : 0;
        for (int b = 0; b < 256; b++)
            bds_rate[b] = span > 0 ? (float)((AS.bds[b] - snap[old].bds[b]) * 60.0 / span) : 0;
    }

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

    frozen_view = ri->frozen;
    if (ri->graph) {
        int rows = term_rows();
        draw_graph(cols, rows - 1);           /* last row: status / filter input */
        fputs("\033[J", stdout);               /* clear below the chart, then the last row */
        printf("\033[%d;1H", rows);
        draw_status_line(ri, true);
        fflush(stdout);
        return;
    }

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
        printf(RED B "  CLIPPING %.0f s ago" R RED " - lower gain: '-' VGA -2 dB, 'l' LNA -8 dB, 'p' amp off" R EL,
               now - last_clip_t);
    else if (weak)
        printf(YEL B "  WEAK SIGNAL" R YEL " - strongest pulse in %.0f s only %.1f dB above noise; raise gain: "
               "'+' VGA +2 dB, 'L' LNA +8 dB, 'p' amp on" R EL, win, peak_w - st.noise_dbfs);
    else if (span < 0.8 * win)
        printf(DIM "  level check: collecting %.0f/%.0f s since start or last gain change" R EL, span, win);
    else printf(GRN "  level OK" R DIM " - peak %.1f dB above noise over %.0f s, no clipping" R EL,
                peak_w - st.noise_dbfs, win);
    {
        char fl[360];
        int o = snprintf(fl, sizeof(fl), "Filters: %s", ri->filters[0] ? ri->filters : "none");
        if (AS.filt_radar || AS.filt_modes || AS.filt_ac || AS.filt_rollcall || AS.filt_aircraft)
            snprintf(fl + o, sizeof(fl) - o, "  (dropped: radar %ld, non-S %ld, not plain A/C %ld, all-call %ld, aircraft %ld)",
                     AS.filt_radar, AS.filt_modes, AS.filt_ac, AS.filt_rollcall, AS.filt_aircraft);
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

    /* left: per kind table; middle: per BDS register table; right: radars and misc */
    char left[24][80], mid[24][96], right[24][160];
    int nl = 0, nm = 0, nr = 0;
    snprintf(left[nl++], 80, B "%-24s %8s %7s %7s %6s" R, "Interrogation", "total", "in-beam", "sidelb", "/min");
    for (int k = 0; k < K_COUNT && nl < 24; k++) {
        if (!AS.total[k] || k == K_UF0 || k == K_UF16) continue;
        snprintf(left[nl++], 80, "%-24s %8ld %7ld %7ld %6.0f", kind_label(k), AS.total[k], AS.inbeam[k],
                 AS.side[k], rate[k]);
    }
    if (nl == 1) snprintf(left[nl++], 80, DIM "(no radar interrogations yet)" R);

    snprintf(mid[nm++], 96, B "%-14s %7s %7s %6s %5s" R, "Comm-B BDS", "total", "in-beam", "sidelb", "/min");
    for (int b = 0; b < 256 && nm < 24; b++) {
        if (!AS.bds[b]) continue;
        const char *ab = bds_abbrev(b);
        char lab[16];
        snprintf(lab, sizeof(lab), "%X,%X %s", b >> 4, b & 15, ab ? ab : "");
        snprintf(mid[nm++], 96, "%-14s %7ld %7ld %6ld %5.0f", lab, AS.bds[b], AS.bds_inbeam[b], AS.bds_side[b], bds_rate[b]);
    }
    if (nm == 1) snprintf(mid[nm++], 96, DIM "(no Comm-B requests yet)" R);

    snprintf(right[nr++], 160, B "Radar  bursts rotation  dwell  per burst" R);
    for (int r = 0; r < AS.nradars && nr < 12; r++) {
        struct radar *X = &AS.radars[r];
        int np = X->np < 64 ? X->np : 64;
        char rot[16] = "      -";
        if (np) snprintf(rot, sizeof(rot), "%7.3fs", median_d(X->periods, np));
        snprintf(right[nr++], 160, "%-6s %6d %s %5.1fms %6.1f", X->key, X->nbursts + X->npartial, rot,
                 X->nbursts ? 1e3 * X->dwell_sum / X->nbursts : 0,
                 X->nbursts ? (double)X->total_sum / X->nbursts : 0);
    }
    if (AS.nradars == 0) snprintf(right[nr++], 160, DIM "(waiting for a main-beam burst)" R);
    int rows = nl > nr ? nl : nr;
    if (nm > rows) rows = nm;
    for (int i = 0; i < rows; i++) {
        /* columns are padded by visible width (escape codes excluded) */
        const char *l = i < nl ? left[i] : "", *mm = i < nm ? mid[i] : "";
        printf("%s%*s  %s%*s  %s" EL, l, pad_to(l, 58), "", mm, pad_to(mm, 46), "", i < nr ? right[i] : "");
    }
    fputs(EL, stdout);

    /* recent bursts */
    printf(B "%-5s %9s %-5s %7s %6s %4s %5s %7s", "Burst", "t[s]", "radar", "dur ms", "peak", "gain", "n", "P[s]");
    for (int g = 0; g < G_COUNT; g++) printf(" %*s", tbl_w(g), tbl_label(g));
    printf(R EL);
    int nb = AS.nhist < 8 ? AS.nhist : 8;
    for (int i = AS.nhist - nb; i < AS.nhist; i++) {
        struct burst *b = &AS.hist[i % BURST_HIST];
        printf("#%-4d %9.3f %-4s%s %7.1f %6.1f %4.0f %5d", b->id, b->t0, b->key, b->partial ? "*" : " ",
               1e3 * (b->t1 - b->t0), b->peak, b->gain, b->n);
        if (b->period > 0) printf(" %7.3f", b->period); else printf(" %7s", "-");
        for (int g = 0; g < G_COUNT; g++) {
            if (b->gcnt[g]) printf(" %*d", tbl_w(g), b->gcnt[g]);
            else printf(" " DIM "%*s" R, tbl_w(g), ".");
        }
        fputs(EL, stdout);
    }
    for (int i = nb; i < 8; i++) fputs(EL, stdout);

    if (AS.nhist > 0) {
        struct burst *b = &AS.hist[(AS.nhist - 1) % BURST_HIST];
        printf(B "Last burst #%d sequence" R DIM
               "  (A C A/C-only all-call, AA CC plain A/C, S UF11 all-call, alt/sqk UF4/5, acid/selAlt/trk/hdgSpd BDS 2,0/4,0/5,0/6,0, S-msg Comm-A, #N repeat)" R EL,
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
    fputs("\033[J", stdout);                   /* clear below the dashboard, then the last row */
    printf("\033[%d;1H", term_rows());
    draw_status_line(ri, true);
    fflush(stdout);
}
