/* ssr1030 - live SSR / Mode S 1030 MHz interrogation analyser for HackRF. */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <getopt.h>
#include <pthread.h>
#include <termios.h>
#include <time.h>
#include <sys/stat.h>
#include <libhackrf/hackrf.h>
#include "analysis.h"

#define RB_SIZE   (128u << 20)       /* 4 s of 16 MS/s int8 IQ */
#define CHUNK     (256 * 1024)
#define GUARD_S   0.06               /* retune settling + samples in flight in USB */
#define ADSB_FREQ 1090e6

static volatile sig_atomic_t stop;
static void on_sig(int s) { (void)s; stop = 1; }

static struct {
    double freq, fs;
    int lna, vga, amp;
    float thr_db, thr_abs;
    int dens_max;
    const char *file, *record, *log, *serial, *adsb_db;
    char sep;
    bool list, dashboard, quiet, show_tcas;
    double duration;
    bool adsb;
    double adsb_every, adsb_dwell;
    bool lol;
    double lol_lat, lol_lon;
    int lol_radius;
    bool agc;
    struct analysis_cfg acfg;
} opt = {
    .freq = 1030e6, .fs = 16e6, .lna = 32, .vga = 24, .amp = 0,
    .thr_db = 12, .dens_max = 50, .sep = ',',
    .adsb_every = 30, .adsb_dwell = 2.9, .lol_radius = 150,
    .agc = true,
    .acfg = { .burst_gap = 0.030, .burst_min = 3, .member_db = 3, .core_db = 6,
              .radar_sel = -1 },
};

/* ---------- ring buffer between libhackrf thread and main thread ---------- */
static uint8_t *rb;
static volatile uint64_t rb_w, rb_r;
static volatile uint64_t drops;
static pthread_mutex_t rb_mx = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t rb_cv = PTHREAD_COND_INITIALIZER;
static FILE *recf;

static int rx_cb(hackrf_transfer *t)
{
    size_t n = t->valid_length;
    pthread_mutex_lock(&rb_mx);
    if (rb_w - rb_r + n > RB_SIZE) {
        drops++;
    } else {
        size_t off = rb_w % RB_SIZE, first = n < RB_SIZE - off ? n : RB_SIZE - off;
        memcpy(rb + off, t->buffer, first);
        memcpy(rb, t->buffer + first, n - first);
        rb_w += n;
        pthread_cond_signal(&rb_cv);
    }
    pthread_mutex_unlock(&rb_mx);
    return stop ? -1 : 0;
}

static uint64_t captured(void)
{
    pthread_mutex_lock(&rb_mx);
    uint64_t w = rb_w;
    pthread_mutex_unlock(&rb_mx);
    return w;
}

/* ---------- frequency switching and stream routing ---------- */
enum { SEG_1030, SEG_1090, SEG_SKIP };
struct mark { uint64_t pos; int to1090; };
static struct mark marks[64];
static int marks_w, marks_r;
static int route_mode = SEG_1030;
static uint64_t guard_end;
static bool in_gap, in_adsb;
static double gap_t0;
static uint8_t zeros[CHUNK];

/* ---------- gain control ----------
 * Nominal gain = amp (11 dB) + LNA + VGA. A single number is mapped onto the
 * three stages so that VGA stays >= 16 dB while LNA can still move, and the
 * RF amp is used only above 80 dB (it overloads first near a radar). */
#define AMP_DB      11
#define GAIN_MAX    (AMP_DB + 40 + 62)
#define AGC_TARGET  -6.0f      /* beam peak, dBFS */
#define AGC_HIGH    -2.0f
#define AGC_LOW     -14.0f

struct gain { int amp, lna, vga; };
static struct gain cur_gain;
static bool adsb_active;
static int gain_changes, agc_low_count;
static double gain_change_t = -1, now_t;
static int gain_delta;
static char gain_why[64];
static float win_peak = -99;          /* strongest radar interrogation since last AGC step */
static uint64_t agc_clip0, agc_samp0;

static int gain_total(struct gain g) { return g.amp * AMP_DB + g.lna + g.vga; }

static struct gain gain_from_total(int t)
{
    if (t < 0) t = 0;
    if (t > GAIN_MAX) t = GAIN_MAX;
    struct gain g;
    g.amp = t > 80;
    int r = t - g.amp * AMP_DB;
    g.lna = r > 16 ? ((r - 16) / 8) * 8 : 0;
    if (g.lna > 40) g.lna = 40;
    g.vga = r - g.lna;
    if (g.vga > 62) g.vga = 62;
    g.vga &= ~1;
    return g;
}

static void hw_gain(hackrf_device *dev, struct gain g)
{
    if (!dev) return;
    hackrf_set_amp_enable(dev, g.amp);
    hackrf_set_lna_gain(dev, g.lna);
    hackrf_set_vga_gain(dev, g.vga);
}

static void set_gain(hackrf_device *dev, struct gain g, const char *why)
{
    int delta = gain_total(g) - gain_total(cur_gain);
    cur_gain = g;
    opt.amp = g.amp; opt.lna = g.lna; opt.vga = g.vga;
    if (!adsb_active) hw_gain(dev, g);
    analysis_set_gain(gain_total(g));
    if (delta) {
        analysis_gain_shift((float)delta);
        gain_changes++;
        gain_change_t = now_t;
        gain_delta = delta;
        snprintf(gain_why, sizeof(gain_why), "%s", why);
    }
}

/* called when a burst of the main radar was committed */
static double last_clip_seen = -1e9;   /* no gain increase within 30 s of clipping */

static void agc_burst(hackrf_device *dev, const struct burst *b, struct demod *d)
{
    struct demod_stats st;
    demod_get_stats(d, &st, false);
    uint64_t ds = st.samples - agc_samp0, dc = st.clipped - agc_clip0;
    agc_samp0 = st.samples; agc_clip0 = st.clipped;
    if (!opt.agc || b->partial) return;
    bool clip = ds && (double)dc / ds > 2e-6;
    if (clip) last_clip_seen = now_t;
    int t = gain_total(cur_gain);
    char why[64];
    if (b->peak > AGC_HIGH || clip) {
        int step = (int)lroundf((b->peak - AGC_TARGET) / 2) * 2;
        if (step < 2) step = clip ? 4 : 2;
        if (step > 16) step = 16;
        snprintf(why, sizeof(why), "beam peak %.1f dBFS%s", b->peak, clip ? ", clipping" : "");
        set_gain(dev, gain_from_total(t - step), why);
        agc_low_count = 0;
    } else if (b->peak < AGC_LOW && now_t - last_clip_seen > 30) {
        if (++agc_low_count >= 2) {
            int step = (int)lroundf((AGC_TARGET - b->peak) / 2) * 2;
            if (step > 12) step = 12;
            snprintf(why, sizeof(why), "beam peak %.1f dBFS", b->peak);
            if (t < GAIN_MAX) set_gain(dev, gain_from_total(t + step), why);
            agc_low_count = 0;
        }
    } else agc_low_count = 0;
    win_peak = -99;
}

/* Clipping from any source: if 2 consecutive 1 s windows clip, lower the
 * gain by 4 dB, but only while no beam pass is open and at most every 2 s. */

static void agc_clip_guard(hackrf_device *dev, struct demod *d)
{
    static double win_t0 = -1, last_step = -1e9;
    static uint64_t c0, s0;
    static int hot;
    struct demod_stats st;
    demod_get_stats(d, &st, false);
    if (win_t0 < 0) { win_t0 = now_t; c0 = st.clipped; s0 = st.samples; return; }
    if (now_t - win_t0 < 1.0) return;
    uint64_t ds = st.samples - s0, dc = st.clipped - c0;
    win_t0 = now_t; c0 = st.clipped; s0 = st.samples;
    bool clip = ds && (double)dc / ds > 1e-6;
    if (clip) last_clip_seen = now_t;
    hot = clip ? hot + 1 : 0;
    if (!opt.agc || hot < 2 || AS.ob.open || now_t - last_step < 2.0 || gain_total(cur_gain) <= 0) return;
    char why[64];
    snprintf(why, sizeof(why), "clipping %.4f%% for %d s", 100.0 * dc / ds, hot);
    set_gain(dev, gain_from_total(gain_total(cur_gain) - 4), why);
    last_step = now_t;
    hot = 0;
}

/* no burst of the main radar for a while: raise gain while nothing is strong */
static void agc_idle(hackrf_device *dev, double since_burst, double rot)
{
    static double last_idle;
    double lim = rot > 0 ? 2 * rot : 10;
    if (lim < 10) lim = 10;
    if (!opt.agc || since_burst < lim || now_t - last_idle < lim || now_t - last_clip_seen < 30) return;
    last_idle = now_t;
    if (win_peak < -20 && gain_total(cur_gain) < GAIN_MAX) {
        int step = win_peak < -90 ? 12 : (int)lroundf((AGC_TARGET - win_peak) / 2) * 2;
        if (step < 4) step = 4;
        if (step > 12) step = 12;
        char why[64];
        if (win_peak < -90) snprintf(why, sizeof(why), "nothing heard for %.0f s", since_burst);
        else snprintf(why, sizeof(why), "no beam pass for %.0f s, peak %.1f dBFS", since_burst, win_peak);
        set_gain(dev, gain_from_total(gain_total(cur_gain) + step), why);
    }
    win_peak = -99;
}

static void set_mark(hackrf_device *dev, bool to1090)
{
    hackrf_set_freq(dev, (uint64_t)(to1090 ? ADSB_FREQ : opt.freq));
    marks[marks_w % 64].pos = captured() & ~(uint64_t)1;
    marks[marks_w % 64].to1090 = to1090;
    marks_w++;
}

static void route(struct demod *d, int type, const uint8_t *p, size_t n)
{
    struct demod_stats st;
    if (type == SEG_1030) {
        if (in_gap) {
            demod_get_stats(d, &st, false);
            analysis_gap(gap_t0, st.samples / opt.fs);
            in_gap = false;
        }
        in_adsb = false;
        demod_process(d, (const int8_t *)p, n / 2);
        if (recf) fwrite(p, 1, n, recf);
        return;
    }
    if (!in_gap) {
        demod_get_stats(d, &st, false);
        gap_t0 = st.samples / opt.fs;
        in_gap = true;
    }
    demod_skip(d, n / 2);
    if (recf) for (size_t k = 0; k < n; k += CHUNK) fwrite(zeros, 1, n - k < CHUNK ? n - k : CHUNK, recf);
    if (type == SEG_1090) {
        if (!in_adsb) { adsb_window_start(); in_adsb = true; }
        adsb_process((const int8_t *)p, n / 2);
    }
}

static void consume(struct demod *d, const uint8_t *buf, size_t n, uint64_t pos0)
{
    uint64_t p = pos0, end = pos0 + n;
    uint64_t guard = (uint64_t)(GUARD_S * opt.fs) * 2;
    while (p < end) {
        while (marks_r < marks_w && marks[marks_r % 64].pos <= p) {
            route_mode = marks[marks_r % 64].to1090 ? SEG_1090 : SEG_1030;
            guard_end = marks[marks_r % 64].pos + guard;
            marks_r++;
        }
        uint64_t seg = end;
        if (marks_r < marks_w && marks[marks_r % 64].pos < seg) seg = marks[marks_r % 64].pos;
        int type = route_mode;
        if (p < guard_end) { type = SEG_SKIP; if (guard_end < seg) seg = guard_end; }
        route(d, type, buf + (p - pos0), seg - p);
        p = seg;
    }
}

/* ---------- interrogation sink ---------- */
static double start_epoch;

static void sink(const struct interrog *it)
{
    if (opt.log) output_write(it, start_epoch);
    if (opt.list && (!it->tcas || opt.show_tcas)) print_interrog_line(stdout, it);
}

static void on_interrog(struct interrog *it, void *ctx)
{
    (void)ctx;
    it->gain_db = (float)gain_total(cur_gain);
    if (!it->tcas && it->amp_db > win_peak) win_peak = it->amp_db;
    analysis_add(it);
}

/* ---------- keyboard ---------- */
static struct termios tio_saved;
static bool tio_raw;

static void kb_raw(void)
{
    if (!isatty(0)) return;
    struct termios t;
    tcgetattr(0, &tio_saved);
    t = tio_saved;
    t.c_lflag &= ~(ICANON | ECHO);
    t.c_cc[VMIN] = 0;
    t.c_cc[VTIME] = 0;
    tcsetattr(0, TCSANOW, &t);
    tio_raw = true;
}

static void kb_restore(void)
{
    if (tio_raw) tcsetattr(0, TCSANOW, &tio_saved);
}

static bool parse_aircraft(const char *s, char *bad, size_t badlen)
{
    struct ac_sel tmp[AC_MAX];
    int n = ac_list_parse(s, tmp, bad, badlen);
    if (n < 0) return false;
    memcpy(opt.acfg.ac, tmp, sizeof(tmp));
    opt.acfg.nac = n;
    return true;
}

static void build_filters(char *buf, size_t n)
{
    int o = 0;
    char rc[8];
    buf[0] = 0;
    if (opt.acfg.radar_sel >= 0) {
        radar_code_str(opt.acfg.radar_sel, rc, sizeof(rc));
        o += snprintf(buf + o, n - o, "radar %s  ", rc);
    }
    if (opt.acfg.modes_only) o += snprintf(buf + o, n - o, "Mode S only  ");
    if (opt.acfg.ac_only) o += snprintf(buf + o, n - o, "plain A/C only  ");
    if (opt.acfg.rollcall_only) o += snprintf(buf + o, n - o, "roll-calls only  ");
    if (opt.acfg.nac) {
        char l[160];
        ac_list_str(&opt.acfg, false, l, sizeof(l));
        o += snprintf(buf + o, n - o, "aircraft %s  ", l);
    }
}

/* ---------- filter input line ('r' radar, 'a' aircraft; empty clears) ---------- */
static bool show_graph;               /* 'b': burst graph screen */
static bool frozen;                   /* space: freeze the screen */
static bool force_full;               /* one full redraw (screen switched while frozen) */
static char edit_kind;                /* 'r', 'a' or 0 when not editing */
static char cmd_buf[160];
static int cmd_len;
static char cmd_msg[160];

static void apply_filters(struct run_info *ri, const char *what)
{
    analysis_init(&opt.acfg);                 /* statistics restart with the new filter */
    analysis_set_gain(gain_total(cur_gain));
    build_filters(ri->filters, sizeof(ri->filters));
    snprintf(cmd_msg, sizeof(cmd_msg), "%s - statistics reset", what);
}

static void apply_input(struct run_info *ri, char kind, char *arg)
{
    while (*arg == ' ') arg++;
    char *e = arg + strlen(arg);
    while (e > arg && e[-1] == ' ') *--e = 0;
    char m[160];
    if (kind == 'r') {
        if (!*arg) { opt.acfg.radar_sel = -1; apply_filters(ri, "radar filter cleared"); return; }
        int code = radar_code_parse(arg);
        if (code < 0) { snprintf(cmd_msg, sizeof(cmd_msg), "bad radar code '%s' (SI1-63 or II0-15), filter unchanged", arg); return; }
        opt.acfg.radar_sel = code;
        char rc[8];
        radar_code_str(code, rc, sizeof(rc));
        snprintf(m, sizeof(m), "radar filter %s", rc);
        apply_filters(ri, m);
    } else {
        if (!*arg) {
            opt.acfg.nac = 0;
            apply_filters(ri, "aircraft filter cleared");
            return;
        }
        char bad[64];
        if (!parse_aircraft(arg, bad, sizeof(bad))) {
            snprintf(cmd_msg, sizeof(cmd_msg), "bad aircraft %s (hex or callsign list), filter unchanged", bad);
            return;
        }
        char l[120];
        ac_list_str(&opt.acfg, true, l, sizeof(l));
        snprintf(m, sizeof(m), "aircraft filter %s", l);
        apply_filters(ri, m);
    }
}

/* read the rest of an escape sequence; returns 0 for a bare Esc, else a
 * key code: 'U' up, 'D' down, 'P' PgUp, 'N' PgDn, 'H' Home, 'E' End, '?' other */
static int eat_escape(void)
{
    char c, num = 0;
    if (read(0, &c, 1) != 1) return 0;
    if (c != '[' && c != 'O') return '?';
    while (read(0, &c, 1) == 1) {
        if (c >= '0' && c <= '9') { num = c; continue; }
        if (c == 'A') return 'U';
        if (c == 'B') return 'D';
        if (c == 'H') return 'H';
        if (c == 'F') return 'E';
        if (c == '~') return num == '5' ? 'P' : num == '6' ? 'N' : (num == '1' || num == '7') ? 'H' :
                              (num == '4' || num == '8') ? 'E' : '?';
        if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) return '?';
    }
    return '?';
}

/* returns true when the screen should be redrawn now */
static bool handle_keys(hackrf_device *dev, struct run_info *ri)
{
    char c;
    bool redraw = false;
    while (read(0, &c, 1) == 1) {
        redraw = true;
        if (edit_kind) {
            if (c == '\r' || c == '\n') {
                cmd_buf[cmd_len] = 0;
                char k = edit_kind;
                edit_kind = 0;
                apply_input(ri, k, cmd_buf);
            } else if (c == 27) {
                if (eat_escape() == 0) {           /* bare Esc cancels; arrow keys are ignored */
                    edit_kind = 0;
                    snprintf(cmd_msg, sizeof(cmd_msg), "cancelled, filter unchanged");
                }
            } else if (c == 127 || c == 8) {
                if (cmd_len > 0) cmd_len--;
            } else if (c == 21) {
                cmd_len = 0;                       /* Ctrl-U clears the line */
            } else if (c >= 32 && c < 127 && cmd_len < (int)sizeof(cmd_buf) - 1) {
                cmd_buf[cmd_len++] = c;
            }
            continue;
        }
        if (c == 27) {
            int k = eat_escape();
            if (show_graph && k && k != '?') display_graph_key(k);
            continue;
        }
        if (show_graph && (c == 'j' || c == 'k')) { display_graph_key(c == 'j' ? 'D' : 'U'); continue; }
        if (c == 'r' || c == 'a') {
            edit_kind = c;
            /* start from the current filter so it can be edited */
            cmd_len = 0;
            if (c == 'r' && opt.acfg.radar_sel >= 0) {
                radar_code_str(opt.acfg.radar_sel, cmd_buf, sizeof(cmd_buf));
                cmd_len = strlen(cmd_buf);
            } else if (c == 'a' && opt.acfg.nac) {
                ac_list_str(&opt.acfg, false, cmd_buf, sizeof(cmd_buf));
                cmd_len = strlen(cmd_buf);
            }
            cmd_msg[0] = 0;
            continue;
        }
        struct gain g = cur_gain;
        bool manual = true;
        switch (c) {
        case 'q': case 'Q': stop = 1; manual = false; break;
        case '+': case '=': if (g.vga < 62) g.vga += 2; break;
        case '-': case '_': if (g.vga > 0) g.vga -= 2; break;
        case 'L': if (g.lna < 40) g.lna += 8; break;
        case 'l': if (g.lna > 0) g.lna -= 8; break;
        case 'p': case 'P': g.amp = !g.amp; break;
        case 'g': case 'G': opt.agc = !opt.agc; agc_low_count = 0; manual = false; break;
        case 'S':
            opt.acfg.modes_only = !opt.acfg.modes_only;
            apply_filters(ri, opt.acfg.modes_only ? "Mode S only ON (non-Mode S dropped)" : "Mode S only OFF");
            manual = false;
            break;
        case 'R':
            opt.acfg.rollcall_only = !opt.acfg.rollcall_only;
            apply_filters(ri, opt.acfg.rollcall_only ? "roll-calls only ON (all-calls dropped)" : "roll-calls only OFF");
            manual = false;
            break;
        case 'b': case 'B': show_graph = !show_graph; force_full = true; manual = false; break;
        case ' ':
            frozen = !frozen;
            snprintf(cmd_msg, sizeof(cmd_msg), frozen ? "screen frozen, data collection continues" : "screen live");
            manual = false;
            break;
        case 'x': case 'X':
            analysis_init(&opt.acfg);
            analysis_set_gain(gain_total(cur_gain));
            snprintf(cmd_msg, sizeof(cmd_msg), "statistics reset");
            manual = false;
            break;
        default: manual = false;
        }
        if (manual) {
            opt.agc = false;               /* a manual setting overrides automatic gain */
            set_gain(dev, g, "manual");
        }
    }
    ri->graph = show_graph;
    ri->frozen = frozen;
    ri->cmd_mode = edit_kind != 0;
    ri->cmd_kind = edit_kind;
    snprintf(ri->cmd, sizeof(ri->cmd), "%.*s", cmd_len, cmd_buf);
    snprintf(ri->cmd_msg, sizeof(ri->cmd_msg), "%s", cmd_msg);
    return redraw;
}

static void usage(void)
{
    printf(
"ssr1030 - SSR Mode A/C/S 1030 MHz interrogation analyser (HackRF)\n\n"
"Source:\n"
"  -f, --file PATH        replay int8 IQ (hackrf_transfer -r format) instead of HackRF\n"
"  -w, --record PATH      also save raw 1030 IQ to PATH (32 MB/s; ADS-B windows as zeros)\n"
"  -s, --rate HZ          sample rate (default 16e6; 8e6 minimum for Mode S)\n"
"      --freq HZ          centre frequency (default 1030e6)\n"
"  -l, --lna DB           LNA gain 0-40 step 8 (default 32)\n"
"  -g, --vga DB           VGA gain 0-62 step 2 (default 24)\n"
"  -a, --amp 0|1          RF amp +11 dB (default 0 = off); -a 1 for a distant radar\n"
"      --no-auto-gain     fixed gain; by default the gain follows the main radar's beam peak\n"
"                         (target -6 dBFS, changed only between beam passes; 'g' toggles)\n"
"      --serial SN        HackRF serial number\n"
"  -t, --duration S       stop after S seconds\n\n"
"Filters (apply to statistics, bursts, logs and list):\n"
"  -r, --radar CODE       keep only one radar: SI1-63 or II0-15; A/C interrogations are\n"
"                         attributed via the beam burst or the preceding UF11, else dropped\n"
"  -S, --modes-only       drop Mode A/C and intermode (A/C-only, A/C/S all-call) interrogations\n"
"  -C, --ac-only          keep only plain Mode A/C (tokens AA / CC, no P4) as sent by older SSR\n"
"                         radars; drops Mode S and all intermode all-calls (not with -S, -R or -A)\n"
"  -R, --rollcall-only    keep only interrogations addressed to an aircraft (UF4/5/20/21/24/19):\n"
"                         drops UF11 and A/C/intermode all-calls and plain Mode A/C\n"
"  -A, --aircraft LIST    keep only interrogations addressed to these aircraft: up to 16 hex\n"
"                         addresses (48af01) or callsigns (LOT3KM; needs --adsb/--adsblol),\n"
"                         separated by commas or spaces: -A 48af01,LOT3KM,DLH123\n"
"                         Prefix with cs: to force a callsign that looks like hex\n\n"
"ADS-B callsigns (1090 MHz, HackRF retunes for short windows):\n"
"  -B, --adsb             listen on 1090 periodically to map addresses to callsigns\n"
"      --adsb-every S     minimum time between windows (default 30)\n"
"      --adsb-dwell S     window length (default 2.9); windows are placed right after a\n"
"                         main-beam burst when the rotation is known, so no burst is lost\n"
"      --adsb-db PATH     load/save the address-callsign table (works with --file too)\n"
"      --adsblol LAT,LON  at startup, fetch hex/callsign pairs of all aircraft around this\n"
"                         position from api.adsb.lol (sends the coordinates to adsb.lol)\n"
"      --adsblol-radius NM  radius for --adsblol (default 150, max 250)\n\n"
"Detection:\n"
"      --threshold DB     pulse threshold above noise (default 12)\n"
"      --threshold-abs A  absolute pulse threshold, magnitude 0-180 (overrides)\n"
"      --dens-max N       max pulses per ms before A/C detections are rejected (default 50)\n\n"
"Bursts (antenna main beam):\n"
"      --burst-gap MS     gap that ends a burst (default 30)\n"
"      --burst-min N      minimum interrogations per burst (default 3)\n"
"      --member-db DB     in-beam if this far above background (default 3)\n"
"      --core-db DB       a burst needs one this far above background (default 6)\n\n"
"Output:\n"
"  -o, --log PATH         write every interrogation to PATH\n"
"      --space            space-delimited log instead of CSV\n"
"  -L, --list             scrolling one-line-per-interrogation output instead of dashboard\n"
"      --tcas             include TCAS/ACAS in --list output\n"
"  -d, --dashboard        dashboard also in --file mode (default: summary only)\n"
"  -q, --quiet            no live output; summary at exit\n\n"
"Keys (dashboard): q quit, +/- VGA, L/l LNA, p RF amp (any of these: manual gain), g auto gain,\n"
"                  r radar filter, a aircraft filter (type it, Enter; empty clears),\n"
"                  S Mode S only on/off, R roll-calls only on/off, b burst graph, space freeze screen,\n"
"                  x reset statistics\n");
}

int main(int argc, char **argv)
{
    static const struct option lo[] = {
        { "file", 1, 0, 'f' }, { "record", 1, 0, 'w' }, { "rate", 1, 0, 's' }, { "freq", 1, 0, 1 },
        { "lna", 1, 0, 'l' }, { "vga", 1, 0, 'g' }, { "amp", 1, 0, 'a' }, { "serial", 1, 0, 2 },
        { "duration", 1, 0, 't' }, { "threshold", 1, 0, 3 }, { "threshold-abs", 1, 0, 4 },
        { "dens-max", 1, 0, 5 }, { "burst-gap", 1, 0, 6 }, { "burst-min", 1, 0, 7 },
        { "member-db", 1, 0, 8 }, { "core-db", 1, 0, 9 }, { "log", 1, 0, 'o' }, { "space", 0, 0, 10 },
        { "list", 0, 0, 'L' }, { "tcas", 0, 0, 11 }, { "dashboard", 0, 0, 'd' }, { "quiet", 0, 0, 'q' },
        { "radar", 1, 0, 'r' }, { "modes-only", 0, 0, 'S' }, { "aircraft", 1, 0, 'A' },
        { "rollcall-only", 0, 0, 'R' }, { "ac-only", 0, 0, 'C' },
        { "adsb", 0, 0, 'B' }, { "adsb-every", 1, 0, 12 }, { "adsb-dwell", 1, 0, 13 },
        { "adsb-db", 1, 0, 14 }, { "adsblol", 1, 0, 15 }, { "adsblol-radius", 1, 0, 16 },
        { "no-auto-gain", 0, 0, 17 },
        { "help", 0, 0, 'h' }, { 0, 0, 0, 0 }
    };
    int c;
    while ((c = getopt_long(argc, argv, "f:w:s:l:g:a:t:o:Ldqr:SCRA:Bh", lo, NULL)) != -1) {
        switch (c) {
        case 'f': opt.file = optarg; break;
        case 'w': opt.record = optarg; break;
        case 's': opt.fs = atof(optarg); break;
        case 1: opt.freq = atof(optarg); break;
        case 'l': opt.lna = atoi(optarg); break;
        case 'g': opt.vga = atoi(optarg); break;
        case 'a': opt.amp = atoi(optarg) != 0; break;
        case 2: opt.serial = optarg; break;
        case 't': opt.duration = atof(optarg); break;
        case 3: opt.thr_db = atof(optarg); break;
        case 4: opt.thr_abs = atof(optarg); break;
        case 5: opt.dens_max = atoi(optarg); break;
        case 6: opt.acfg.burst_gap = atof(optarg) / 1e3; break;
        case 7: opt.acfg.burst_min = atoi(optarg); break;
        case 8: opt.acfg.member_db = atof(optarg); break;
        case 9: opt.acfg.core_db = atof(optarg); break;
        case 'o': opt.log = optarg; break;
        case 10: opt.sep = ' '; break;
        case 'L': opt.list = true; break;
        case 11: opt.show_tcas = true; break;
        case 'd': opt.dashboard = true; break;
        case 'q': opt.quiet = true; break;
        case 'r':
            opt.acfg.radar_sel = radar_code_parse(optarg);
            if (opt.acfg.radar_sel < 0) { fprintf(stderr, "bad radar code '%s' (use SI1-63 or II0-15)\n", optarg); return 1; }
            break;
        case 'S': opt.acfg.modes_only = true; break;
        case 'R': opt.acfg.rollcall_only = true; break;
        case 'C': opt.acfg.ac_only = true; break;
        case 'A':
            {
                char bad[64];
                if (!parse_aircraft(optarg, bad, sizeof(bad))) { fprintf(stderr, "bad aircraft %s\n", bad); return 1; }
            }
            break;
        case 'B': opt.adsb = true; break;
        case 12: opt.adsb_every = atof(optarg); break;
        case 13: opt.adsb_dwell = atof(optarg); break;
        case 14: opt.adsb_db = optarg; break;
        case 15:
            if (sscanf(optarg, "%lf,%lf", &opt.lol_lat, &opt.lol_lon) != 2 ||
                fabs(opt.lol_lat) > 90 || fabs(opt.lol_lon) > 180) {
                fprintf(stderr, "bad coordinates '%s' (use LAT,LON in decimal degrees, e.g. 52.23,21.01)\n", optarg);
                return 1;
            }
            opt.lol = true;
            break;
        case 16: opt.lol_radius = atoi(optarg); break;
        case 17: opt.agc = false; break;
        default: usage(); return c == 'h' ? 0 : 1;
        }
    }
    if (opt.fs < 8e6) { fprintf(stderr, "sample rate must be >= 8e6 for Mode S DPSK\n"); return 1; }
    if (opt.acfg.ac_only && (opt.acfg.modes_only || opt.acfg.rollcall_only || opt.acfg.nac)) {
        fprintf(stderr, "-C (plain A/C only) cannot be combined with -S, -R or -A: nothing would pass\n");
        return 1;
    }
    bool live = opt.file == NULL;
    bool dash = !opt.list && !opt.quiet && (live || opt.dashboard);
    if (opt.adsb && !live) { fprintf(stderr, "--adsb needs the HackRF; ignored in --file mode\n"); opt.adsb = false; }
    bool any_cs = false;
    for (int k = 0; k < opt.acfg.nac; k++) if (opt.acfg.ac[k].addr < 0) any_cs = true;
    if (any_cs && !opt.adsb && !opt.adsb_db && !opt.lol)
        fprintf(stderr, "warning: callsign filter without --adsb, --adsblol or --adsb-db can never resolve\n");
    if (opt.lol_radius < 1) opt.lol_radius = 1;
    if (opt.lol_radius > 250) opt.lol_radius = 250;            /* API maximum */

    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);

    if (opt.file) opt.agc = false;
    cur_gain.amp = opt.amp; cur_gain.lna = opt.lna; cur_gain.vga = opt.vga;
    adsb_init(opt.fs);
    if (opt.adsb_db) {
        int n = cs_load(opt.adsb_db);
        if (n >= 0) fprintf(stderr, "loaded %d callsigns from %s\n", n, opt.adsb_db);
    }
    if (opt.lol) {
        int n = cs_fetch_adsblol(opt.lol_lat, opt.lol_lon, opt.lol_radius);
        if (n >= 0) fprintf(stderr, "adsb.lol: %d callsigns within %d NM of %.4f,%.4f\n", n, opt.lol_radius,
                            opt.lol_lat, opt.lol_lon);
        else fprintf(stderr, "adsb.lol: request failed (no network or API error), continuing without\n");
    }
    analysis_init(&opt.acfg);
    analysis_set_gain(gain_total(cur_gain));
    analysis_set_sink(sink);
    struct demod *d = demod_new(opt.fs, opt.thr_db, opt.thr_abs, opt.dens_max, on_interrog, NULL);

    int fd = -1;
    hackrf_device *dev = NULL;
    struct run_info ri = { .fs = opt.fs, .freq = opt.freq, .live = live, .adsb = opt.adsb };
    build_filters(ri.filters, sizeof(ri.filters));
    char srcname[300];

    if (live) {
        int r = hackrf_init();
        if (r == HACKRF_SUCCESS) r = hackrf_open_by_serial(opt.serial, &dev);
        if (r != HACKRF_SUCCESS) {
            fprintf(stderr, "HackRF open failed: %s\n", hackrf_error_name(r));
            return 1;
        }
        hackrf_set_sample_rate(dev, opt.fs);
        ri.bbw = hackrf_compute_baseband_filter_bw((uint32_t)(0.75 * opt.fs));
        hackrf_set_baseband_filter_bandwidth(dev, (uint32_t)ri.bbw);
        {
            uint8_t bid = BOARD_ID_UNDETECTED;
            char ver[256] = "?";
            read_partid_serialno_t ps;
            memset(&ps, 0, sizeof(ps));
            hackrf_board_id_read(dev, &bid);
            hackrf_version_string_read(dev, ver, sizeof(ver) - 1);
            hackrf_board_partid_serialno_read(dev, &ps);
            snprintf(ri.hw, sizeof(ri.hw), "%s fw %s SN ..%08x%08x", hackrf_board_id_name(bid), ver,
                     ps.serial_no[2], ps.serial_no[3]);
        }
        hackrf_set_freq(dev, (uint64_t)opt.freq);
        hackrf_set_lna_gain(dev, opt.lna);
        hackrf_set_vga_gain(dev, opt.vga);
        hackrf_set_amp_enable(dev, opt.amp);
        rb = malloc(RB_SIZE);
        if (opt.record && !(recf = fopen(opt.record, "wb"))) { perror(opt.record); return 1; }
        start_epoch = (double)time(NULL);
        r = hackrf_start_rx(dev, rx_cb, NULL);
        if (r != HACKRF_SUCCESS) { fprintf(stderr, "start_rx: %s\n", hackrf_error_name(r)); return 1; }
        snprintf(srcname, sizeof(srcname), "HackRF %s", opt.serial ? opt.serial : "");
    } else {
        fd = open(opt.file, O_RDONLY);
        if (fd < 0) { perror(opt.file); return 1; }
        struct stat sb;
        fstat(fd, &sb);
        start_epoch = (double)sb.st_mtime - sb.st_size / 2.0 / opt.fs;
        snprintf(srcname, sizeof(srcname), "file %s", opt.file);
    }
    ri.source = srcname;
    if (opt.log && output_open(opt.log, opt.sep) < 0) { perror(opt.log); return 1; }

    if (dash) { display_init(); kb_raw(); }
    else if (live && isatty(0)) kb_raw();

    uint8_t *buf = malloc(CHUNK);
    struct timespec last_draw = { 0 };
    double now = 0, last_main_burst = 0;
    uint64_t adsb_pos = 0;
    double adsb_last_end = -1e9;
    while (!stop) {
        size_t n = 0;
        uint64_t pos0 = 0;
        if (live) {
            pthread_mutex_lock(&rb_mx);
            while (rb_w == rb_r && !stop) {
                struct timespec ts;
                clock_gettime(CLOCK_REALTIME, &ts);
                ts.tv_nsec += 200000000;
                if (ts.tv_nsec >= 1000000000) { ts.tv_sec++; ts.tv_nsec -= 1000000000; }
                pthread_cond_timedwait(&rb_cv, &rb_mx, &ts);
                if (rb_w == rb_r) break;
            }
            n = rb_w - rb_r;
            if (n > CHUNK) n = CHUNK;
            size_t off = rb_r % RB_SIZE, first = n < RB_SIZE - off ? n : RB_SIZE - off;
            memcpy(buf, rb + off, first);
            memcpy(buf + first, rb, n - first);
            pthread_mutex_unlock(&rb_mx);
            pos0 = rb_r;
            rb_r += n;
            if (!hackrf_is_streaming(dev) && n == 0) { fprintf(stderr, "HackRF stopped streaming\n"); break; }
        } else {
            ssize_t r = read(fd, buf, CHUNK);
            if (r <= 0) break;
            n = (size_t)r & ~(size_t)1;
        }
        if (n) {
            if (live) consume(d, buf, n, pos0);
            else demod_process(d, (const int8_t *)buf, n / 2);
            struct demod_stats st;
            demod_get_stats(d, &st, false);
            now = st.samples / opt.fs;
            analysis_set_floor(st.thr_dbfs);
            analysis_tick(now);
            now_t = now;
            if (live) {
                static int seen_hist;
                int code = analysis_main_code();
                if (seen_hist > AS.nhist) seen_hist = AS.nhist;      /* statistics were reset */
                while (seen_hist < AS.nhist) {
                    const struct burst *b = &AS.hist[seen_hist % BURST_HIST];
                    if (code != -2 && b->code == code) { agc_burst(dev, b, d); last_main_burst = now; }
                    seen_hist++;
                }
                agc_idle(dev, now - last_main_burst, analysis_rotation());
                agc_clip_guard(dev, d);
                if (now - last_main_burst > 3600) last_main_burst = now;   /* keep agc_idle periodic */
            }
        }
        if (opt.duration > 0 && now >= opt.duration) break;

        /* ADS-B window scheduling, in capture time */
        if (opt.adsb) {
            double cap_t = captured() / 2.0 / opt.fs;
            double rot = analysis_rotation();
            ri.adsb_now = adsb_active;
            if (!adsb_active) {
                double wait = opt.adsb_every - (cap_t - adsb_last_end);
                ri.adsb_next = wait > 0 ? wait : -1;
                if (wait <= 0) {
                    bool go = true;
                    double lb = analysis_last_burst_end();
                    double win = opt.adsb_dwell + 2 * GUARD_S;
                    if (rot > win + 1.0 && lb >= 0) {
                        double ph = fmod(cap_t - lb, rot);
                        go = ph > 0.02 && ph < rot - win - 0.4;   /* next beam pass well after the window */
                    }
                    if (-wait > 3 * rot + 10) go = true;              /* never starve the table */
                    if (go) {
                        set_mark(dev, true);
                        adsb_active = true;
                        adsb_pos = captured();
                    }
                }
            } else if (captured() - adsb_pos >= (uint64_t)((opt.adsb_dwell + GUARD_S) * opt.fs) * 2) {
                set_mark(dev, false);
                adsb_active = false;
                hw_gain(dev, cur_gain);       /* apply a gain change deferred during the window */
                adsb_last_end = captured() / 2.0 / opt.fs;
                if (opt.adsb_db) cs_save(opt.adsb_db);
            }
        }

        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        double el = (ts.tv_sec - last_draw.tv_sec) + (ts.tv_nsec - last_draw.tv_nsec) / 1e9;
        bool kred = tio_raw && handle_keys(dev, &ri);
        if (el > 0.5 || kred) {
            last_draw = ts;
            ri.lna = opt.lna; ri.vga = opt.vga; ri.amp = opt.amp;
            ri.agc = opt.agc;
            ri.gain_age = gain_change_t >= 0 ? now - gain_change_t : -1;
            ri.gain_delta = gain_delta;
            snprintf(ri.gain_why, sizeof(ri.gain_why), "%s", gain_why);
            /* frozen: the chart (pinned burst) still redraws for scrolling; the
             * dashboard keeps its picture and only the bottom row is updated */
            if (dash && (!frozen || show_graph || force_full)) display_render(&ri, now, d, drops);
            else if (dash && kred) display_status_only(&ri);
            force_full = false;
            if (opt.list) fflush(stdout);
        }
    }

    if (live) {
        stop = 1;
        hackrf_stop_rx(dev);
        hackrf_close(dev);
        hackrf_exit();
        if (recf) fclose(recf);
    }
    demod_flush(d);
    analysis_tick(1e300);
    if (dash) display_restore();
    kb_restore();
    output_close();
    if (opt.adsb_db) cs_save(opt.adsb_db);

    struct demod_stats st;
    demod_get_stats(d, &st, false);
    FILE *sf = opt.list ? stderr : stdout;
    if (live)
        fprintf(sf, "\nHackRF: %s\n  %.3f MHz, %.1f MS/s, baseband filter %.2f MHz, AMP %s, LNA %d dB, VGA %d dB "
                    "(nominal gain %d dB)%s\n",
                ri.hw, opt.freq / 1e6, opt.fs / 1e6, ri.bbw / 1e6, opt.amp ? "on" : "off", opt.lna, opt.vga,
                (opt.amp ? 11 : 0) + opt.lna + opt.vga, "");
    if (live)
        fprintf(sf, "  gain: %s, %d changes during the run%s%s%s\n",
                opt.agc ? "automatic" : "manual", gain_changes, gain_changes ? " (last: " : "",
                gain_changes ? gain_why : "", gain_changes ? ")" : "");
    fprintf(sf, "\nnoise %.1f dBFS, threshold %.1f dBFS, %lu pulses, %lu unmatched, %lu interference-rejected, "
                "clipped samples %.4f%%, USB drops %lu\n",
            st.noise_dbfs, st.thr_dbfs, (unsigned long)st.pulses, (unsigned long)st.stray_pulses,
            (unsigned long)st.garbled, st.samples ? 100.0 * st.clipped / st.samples : 0, (unsigned long)drops);
    if (opt.adsb) {
        struct adsb_stats as;
        adsb_get_stats(&as);
        fprintf(sf, "ADS-B: %ld windows, %.1f s on 1090, %ld frames (DF17/18 %ld, DF11 %ld), "
                    "%ld identifications (+%ld Comm-B BDS2,0), %d callsigns known\n",
                as.windows, as.seconds, as.frames, as.df17, as.df11, as.idents, as.commb_idents, cs_count());
    }
    analysis_print_summary(sf);
    return 0;
}
