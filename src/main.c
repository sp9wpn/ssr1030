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
    struct analysis_cfg acfg;
} opt = {
    .freq = 1030e6, .fs = 16e6, .lna = 32, .vga = 24, .amp = 0,
    .thr_db = 12, .dens_max = 50, .sep = ',',
    .adsb_every = 30, .adsb_dwell = 2.0, .lol_radius = 150,
    .acfg = { .burst_gap = 0.030, .burst_min = 3, .member_db = 3, .core_db = 6,
              .radar_sel = -1, .ac_addr = -1 },
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

static int gain_changes;

static void handle_keys(hackrf_device *dev)
{
    char c;
    while (read(0, &c, 1) == 1) {
        if (c && strchr("+=-_LlaA", c)) gain_changes++;
        switch (c) {
        case 'q': case 'Q': stop = 1; break;
        case '+': case '=': if (opt.vga < 62) opt.vga += 2; if (dev) hackrf_set_vga_gain(dev, opt.vga); break;
        case '-': case '_': if (opt.vga > 0) opt.vga -= 2; if (dev) hackrf_set_vga_gain(dev, opt.vga); break;
        case 'L': if (opt.lna < 40) opt.lna += 8; if (dev) hackrf_set_lna_gain(dev, opt.lna); break;
        case 'l': if (opt.lna > 0) opt.lna -= 8; if (dev) hackrf_set_lna_gain(dev, opt.lna); break;
        case 'a': case 'A': opt.amp = !opt.amp; if (dev) hackrf_set_amp_enable(dev, opt.amp); break;
        case 'r': case 'R': analysis_init(&opt.acfg); break;
        }
    }
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
"      --serial SN        HackRF serial number\n"
"  -t, --duration S       stop after S seconds\n\n"
"Filters (apply to statistics, bursts, logs and list):\n"
"  -r, --radar CODE       keep only one radar, e.g. SI17 or II3; A/C interrogations are\n"
"                         attributed via the beam burst or the preceding UF11, else dropped\n"
"  -S, --modes-only       drop Mode A/C and intermode (A/C-only, A/C/S all-call) interrogations\n"
"  -A, --aircraft ID      keep only interrogations addressed to one aircraft: 6 hex digits\n"
"                         (48AF01) or a callsign (LOT3KM; needs --adsb or --adsb-db).\n"
"                         Prefix with cs: to force a callsign that looks like hex\n\n"
"ADS-B callsigns (1090 MHz, HackRF retunes for short windows):\n"
"  -B, --adsb             listen on 1090 periodically to map addresses to callsigns\n"
"      --adsb-every S     minimum time between windows (default 30)\n"
"      --adsb-dwell S     window length (default 2.0); windows are placed right after a\n"
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
"Keys (dashboard): q quit, +/- VGA, L/l LNA, a amp, r reset statistics\n");
}

static bool parse_aircraft(const char *s)
{
    if (!strncasecmp(s, "cs:", 3)) {
        snprintf(opt.acfg.ac_callsign, sizeof(opt.acfg.ac_callsign), "%s", s + 3);
    } else {
        const char *h = !strncasecmp(s, "0x", 2) ? s + 2 : s;
        bool hex = strlen(h) == 6;
        for (int k = 0; hex && k < 6; k++) hex = isxdigit((unsigned char)h[k]);
        if (hex) opt.acfg.ac_addr = strtol(h, NULL, 16);
        else if (strlen(s) <= 8) snprintf(opt.acfg.ac_callsign, sizeof(opt.acfg.ac_callsign), "%s", s);
        else return false;
    }
    for (char *p = opt.acfg.ac_callsign; *p; p++) *p = toupper((unsigned char)*p);
    return true;
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
        { "adsb", 0, 0, 'B' }, { "adsb-every", 1, 0, 12 }, { "adsb-dwell", 1, 0, 13 },
        { "adsb-db", 1, 0, 14 }, { "adsblol", 1, 0, 15 }, { "adsblol-radius", 1, 0, 16 },
        { "help", 0, 0, 'h' }, { 0, 0, 0, 0 }
    };
    int c;
    while ((c = getopt_long(argc, argv, "f:w:s:l:g:a:t:o:Ldqr:SA:Bh", lo, NULL)) != -1) {
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
        case 'A':
            if (!parse_aircraft(optarg)) { fprintf(stderr, "bad aircraft '%s'\n", optarg); return 1; }
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
        default: usage(); return c == 'h' ? 0 : 1;
        }
    }
    if (opt.fs < 8e6) { fprintf(stderr, "sample rate must be >= 8e6 for Mode S DPSK\n"); return 1; }
    bool live = opt.file == NULL;
    bool dash = !opt.list && !opt.quiet && (live || opt.dashboard);
    if (opt.adsb && !live) { fprintf(stderr, "--adsb needs the HackRF; ignored in --file mode\n"); opt.adsb = false; }
    if (opt.acfg.ac_callsign[0] && !opt.adsb && !opt.adsb_db && !opt.lol)
        fprintf(stderr, "warning: callsign filter without --adsb, --adsblol or --adsb-db can never resolve\n");
    if (opt.lol_radius < 1) opt.lol_radius = 1;
    if (opt.lol_radius > 250) opt.lol_radius = 250;            /* API maximum */

    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);

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
    analysis_set_sink(sink);
    struct demod *d = demod_new(opt.fs, opt.thr_db, opt.thr_abs, opt.dens_max, on_interrog, NULL);

    int fd = -1;
    hackrf_device *dev = NULL;
    struct run_info ri = { .fs = opt.fs, .freq = opt.freq, .live = live, .adsb = opt.adsb };
    {
        int o = 0;
        char rc[8];
        if (opt.acfg.radar_sel >= 0) {
            radar_code_str(opt.acfg.radar_sel, rc, sizeof(rc));
            o += snprintf(ri.filters + o, sizeof(ri.filters) - o, "radar %s  ", rc);
        }
        if (opt.acfg.modes_only) o += snprintf(ri.filters + o, sizeof(ri.filters) - o, "Mode S only  ");
        if (opt.acfg.ac_addr >= 0)
            o += snprintf(ri.filters + o, sizeof(ri.filters) - o, "aircraft %06llX  ", (unsigned long long)opt.acfg.ac_addr);
        if (opt.acfg.ac_callsign[0])
            o += snprintf(ri.filters + o, sizeof(ri.filters) - o, "aircraft %s  ", opt.acfg.ac_callsign);
    }
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
    double now = 0;
    bool adsb_on = false;
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
        }
        if (opt.duration > 0 && now >= opt.duration) break;

        /* ADS-B window scheduling, in capture time */
        if (opt.adsb) {
            double cap_t = captured() / 2.0 / opt.fs;
            double rot = analysis_rotation();
            ri.adsb_now = adsb_on;
            if (!adsb_on) {
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
                    if (go) { set_mark(dev, true); adsb_on = true; adsb_pos = captured(); }
                }
            } else if (captured() - adsb_pos >= (uint64_t)((opt.adsb_dwell + GUARD_S) * opt.fs) * 2) {
                set_mark(dev, false);
                adsb_on = false;
                adsb_last_end = captured() / 2.0 / opt.fs;
                if (opt.adsb_db) cs_save(opt.adsb_db);
            }
        }

        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        double el = (ts.tv_sec - last_draw.tv_sec) + (ts.tv_nsec - last_draw.tv_nsec) / 1e9;
        if (el > 0.5) {
            last_draw = ts;
            if (tio_raw) handle_keys(dev);
            ri.lna = opt.lna; ri.vga = opt.vga; ri.amp = opt.amp;
            if (dash) display_render(&ri, now, d, drops);
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
                (opt.amp ? 11 : 0) + opt.lna + opt.vga, gain_changes ? " - gain was changed during the run" : "");
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
