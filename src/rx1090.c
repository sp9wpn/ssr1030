/* 1090 MHz replies from a second receiver (any Beast source with its own 12 MHz
 * sample timestamps, e.g. dump1090 or readsb with an SDR), matched to the
 * interrogations heard on 1030 MHz.
 *
 * Time: a query has HackRF stream time t, a reply the dump1090 12 MHz
 * sample-clock timestamp r. A reply to that query arrives at
 *     r = t + off(t) + delay
 * off is the offset between the two receiver clocks (it drifts by ppm),
 * delay is the 128 us Mode S turnaround plus the extra path
 * radar -> aircraft -> receiver minus radar -> receiver (0 to ~3.5 ms).
 *  - coarse: both streams are stamped with the host clock on arrival; the
 *    smallest (wall - stream time) of each gives off to a few ms (seen: true
 *    offset 1-2 ms above it); locks search only +-LOCK_SPAN around it
 *  - lock: roll-calls and replies with the same address over 60 s; at the
 *    true offset they all fall within a few ms, so the 4 ms window holding
 *    the most pairs wins
 *  - track: an aircraft's delay barely changes between scans, so the
 *    residual of each re-matched aircraft follows the clock drift
 *  - drift: learned from matches more than 1 s apart and kept across lock
 *    losses (the host-clock minima jitter too much to measure it).
 *    Between rare matches the match window widens with the time since the
 *    last one, and the first match after a long gap resets the offset
 *  - re-check: while locked on positions, the geometry search runs again on
 *    the last minute of pairs; when nothing fits the tracked clock and
 *    another offset is clear, it re-locks there
 * Roll-calls (UF4/5/20/21) match by address, UF11 by the II/SI code in the
 * DF11 reply. Mode A/C replies carry neither and are not matched.
 * The second receiver usually hears only some of the replies (antenna
 * placement), so an unheard reply never counts against the lock; only
 * replies from the aircraft that do not fit the expected delay do.
 * Sync from 1090 windows (the best clock, used whenever it works): the HackRF
 * briefly listens on 1090 MHz (ADS-B windows) and its decoded frames carry
 * HackRF stream times. dump1090 hears the same frames; an identical frame
 * (56/112 bits, CRC good) within SYNC_SPAN of the coarse offset gives
 * off = r - t_preamble + SYNC_REF directly (dump1090 stamps preamble + 64 us;
 * a reply's preamble starts 128 us after the sync reversal, 4.75 us after P1).
 * The median of a window locks the clock to ~0.5 us, consecutive windows give
 * the drift, and network delay does not matter. Between windows the offset is
 * extrapolated; the error grows as the drift itself changes (warm-up), so the
 * window interval starts at SYNC_IV_MIN and doubles while windows agree with
 * the tracked clock. Matching, lock searches and
 * clock updates from replies stay off while these windows keep coming.
 * With the radar and receiver positions known (--radar-pos, --rx-pos), the
 * extra path comes from the aircraft's ADS-B position, so every heard reply
 * gives the offset itself (to ~10 us): the geometry lock. Without
 * positions the histogram lock above is used.
 * dump1090 drops DF4/5/20/21 from addresses it has not seen in a squitter;
 * every roll-called address is therefore announced to its Beast input
 * (port 30004) as a DF11 acquisition squitter, which puts it on its list. */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <pthread.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/time.h>
#include "ssr1030.h"

#define INQ         16384           /* reader thread -> main thread */
#define RP_N        65536           /* replies kept for matching */
#define QH_N        8192            /* roll-calls kept for locking */
#define AC_N        4096            /* per-aircraft delays */
#define WIN_LO      (-0.0003)       /* match window around t + off, s */
#define WIN_HI      0.0040
#define AC_TOL      0.00025         /* known aircraft: delay within this of last time, */
#define AC_TOL0     10e-6           /* ... with the clock from 1090 windows: this */
#define AC_MOVE     2e-6            /* ... plus this per s since (250 m/s there and back: 1.7 us/s) */
#define AC_AGE      60.0            /* s: a known delay is trusted this long */
#define LOCK_SPAN   0.01            /* lock search: coarse offset +- this, s (true offset seen 1-2 ms off) */
#define LOCK_BIN    0.00005
#define LOCK_NB     400             /* 2 * LOCK_SPAN / LOCK_BIN */
#define LOCK_W      80              /* 4 ms window */
#define LOCK_HIST   60.0            /* s of roll-calls used for a lock (few replies are heard) */
#define LOCK_MIN    4               /* pairs needed */
#define LOCK_RATIO  2               /* ... and this many times the best other window */
#define UNLOCK_N    5               /* roll-calls in a row whose aircraft was heard but did not fit */
#define NO_ADDR_RX  0xFFFFFFFFu
#define ANN_N       1024            /* addresses waiting to be announced */
#define ANN_EVERY   20.0            /* s: dump1090 forgets an address after ~60 s */
#define HOLD_UNLOCKED 7.0           /* s: queries wait this long while searching, so the passes that
                                       give the lock are still marked: it comes ~6 s after the first of its 2 scans */
#define BEAST_IN    "30004"
#define GEO_BIN     0.00001         /* geometry lock: 10 us bins */
#define GEO_NB      2000            /* 2 * LOCK_SPAN / GEO_BIN */
#define GEO_W       3               /* true pairs agree within 30 us */
#define GEO_TOL     0.00004         /* match: reply within this of the predicted time */
#define TOL_GROW    1e-6            /* ... plus this per s since the last match (drift error) */
#define TOL_MAX     0.0003
#define RECHECK     5.0             /* s between geometry re-checks while locked */
#define SCAN_GAP    0.5             /* s: pairs further apart than this come from different scans */
#define INDEP_MIN   2               /* a geometry lock needs agreement from this many scans or aircraft, */
#define INDEP_MOVE  3               /* a re-check moving an established lock this many */
#define PAIRS_ONE   6               /* ... or this many agreeing pairs from one source (false ones seen: 2) */
#define SYNC_N      8192            /* dump1090 frames kept for the window sync */
#define HQ_N        4096            /* HackRF window frames waiting for their dump1090 copy */
#define SYNC_SPAN   0.25            /* s: identical frame searched this far from the coarse offset */
#define SYNC_REF    4.75e-6         /* offset = r - t_preamble + this (see the top) */
#define SYNC_MIN    3               /* matched frames needed in a window */
#define SYNC_AGREE  3e-6            /* frames within this of the window median count */
#define SYNC_IV_MIN 30.0            /* s between sync windows while starting or unsure, */
#define SYNC_IV_MAX 300.0           /* doubling up to this while windows agree with the tracked clock */
#define SYNC_GOOD   5e-6            /* a window this close to the tracked clock is good, */
#define SYNC_BAD    15e-6           /* this far off: back to SYNC_IV_MIN */
#define SYNC_STALE  (SYNC_IV_MAX + 120) /* s without a window: back to locking from replies */
#define POS_AGE     10.0            /* s: a position is used this long (250 m/s -> ~8 us) */
#define TURNAROUND  128e-6          /* Mode S reply delay after the sync phase reversal */
#define C_LIGHT     299792458.0
#define RES_N       64
#define REC_N       16384           /* matched replies kept for display */

struct raw { uint64_t ts; double wall; uint8_t len; uint8_t msg[14]; };
struct reply { double r; uint32_t addr; int16_t code; uint8_t df, len; bool used; uint8_t msg[14]; };
struct qrec { double t; uint32_t addr; int kind, code; };    /* roll-calls and UF11 (code) */
struct acd {
    uint32_t addr;
    double delay, t, announced;     /* learned delay (no position), when, last announced to dump1090 */
    bool pos;                       /* ADS-B position: ECEF, receiver time of the fix */
    double ecef[3], pos_r;
    double h, h_t;                  /* height above the radar (reply or ADS-B altitude), when */
    bool used;
};
struct minage { double cur, prev, t; };     /* minimum over the last 30-60 s */

/* reader thread */
static char host[96], port[16];
static pthread_t thr;
static bool thr_on;
static volatile bool quit;
static pthread_mutex_t mx = PTHREAD_MUTEX_INITIALIZER;
static struct raw inq[INQ];
static uint64_t in_w, in_r;
static volatile bool connected, feed_failed;  /* feed_failed: the last try to announce to port 30004 failed */
static uint32_t ann[ANN_N];
static uint64_t ann_w, ann_r;

/* main thread */
static struct reply rp[RP_N];
static uint64_t rp_n;
static struct qrec qh[QH_N];
static uint64_t qh_n;
static struct acd act[AC_N];
static struct minage hw = { 1e300, 1e300, 0 }, rw = { 1e300, 1e300, 0 };
static double last_r = -1, last_rx_wall, t_stream;
static bool locked, lock_geo, lock_frames;
struct sframe { double r; uint8_t len, msg[14]; };      /* dump1090 frame for the window sync */
static struct sframe sf[SYNC_N];
static uint64_t sf_n;
struct hframe { double t; uint8_t len, msg[14]; };      /* HackRF window frame */
static struct hframe hq[HQ_N];
static uint64_t hq_r, hq_w;
static double wo[HQ_N], wt[HQ_N];                       /* current window: offsets, HackRF times */
static int wn;
static double sync_o, sync_t = -1e300;                  /* last window: offset at time */
static double sync_iv = SYNC_IV_MIN;                    /* wanted time between windows */
static int sync_good;                                   /* good windows in a row */
static bool geo_on;
static double g_radar[3], g_rx[3], g_d0, g_ref_lat, g_ref_lon;
static double g_ralt, g_up[3], g_xh, g_xz;  /* radar altitude; antenna in the radar's frame: horizontal, up */
static double res[RES_N];
static struct rx1090_rec recs[REC_N];
static uint64_t nrec;
static long nres;
static double c_off, c_k, c_t0;             /* off(t) = c_off + c_k (t - c_t0) */
static double last_try, last_check, last_match; /* last lock try, geometry re-check, clock update */
static int misfits;                         /* heard-but-unfitting roll-calls in a row */
static struct rx1090_stats S;
static FILE *rlog;                          /* RX1090_LOG: every reply, clock data */
static double rate_t0;
static long rate_n0;

static double wall_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void minage_put(struct minage *m, double v, double wall)
{
    if (wall - m->t > 30) { m->prev = m->cur; m->cur = 1e300; m->t = wall; }
    if (v < m->cur) m->cur = v;
}
static double minage_get(const struct minage *m) { return m->cur < m->prev ? m->cur : m->prev; }

static bool k_learned;                      /* c_k came from matches: keep it across lock losses */

/* drift a new lock starts from */
static double first_drift(void) { return k_learned ? c_k : 0; }

/* ---------------- reader thread: Beast binary frames ---------------- */

static int dial(const char *svc)
{
    struct addrinfo hints = { 0 }, *res, *a;
    hints.ai_socktype = SOCK_STREAM;
    if (getaddrinfo(host, svc, &hints, &res)) return -1;
    int fd = -1;
    for (a = res; a; a = a->ai_next) {
        fd = socket(a->ai_family, a->ai_socktype, a->ai_protocol);
        if (fd < 0) continue;
        if (connect(fd, a->ai_addr, a->ai_addrlen) == 0) break;
        close(fd);
        fd = -1;
    }
    freeaddrinfo(res);
    if (fd >= 0) {
        struct timeval tv = { 0, 300000 };
        setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    }
    return fd;
}

static void push(const uint8_t *f, int len, double wall)
{
    struct raw r;
    r.ts = 0;
    for (int k = 0; k < 6; k++) r.ts = r.ts << 8 | f[k];
    r.wall = wall;
    r.len = (uint8_t)len;
    memcpy(r.msg, f + 7, len);
    pthread_mutex_lock(&mx);
    if (in_w - in_r < INQ) inq[in_w++ % INQ] = r;
    pthread_mutex_unlock(&mx);
}

/* DF11 acquisition squitter (IID 0) for this address, as a Beast frame */
static int squitter_frame(uint32_t a, uint8_t *out)
{
    uint8_t m[7] = { 11 << 3 | 5, a >> 16, a >> 8, a };
    uint32_t pi = modes_crc(m, 56);
    m[4] = pi >> 16; m[5] = pi >> 8; m[6] = pi;
    uint8_t raw[1 + 6 + 1 + 7] = { '2', 0, 0, 0, 0, 0, 0, 0x80 };
    memcpy(raw + 8, m, 7);
    int n = 0;
    out[n++] = 0x1a;
    for (size_t k = 0; k < sizeof(raw); k++) {
        out[n++] = raw[k];
        if (k && raw[k] == 0x1a) out[n++] = 0x1a;
    }
    return n;
}

static void announce_pending(int *fd)
{
    for (;;) {
        pthread_mutex_lock(&mx);
        bool any = ann_r < ann_w;
        uint32_t a = any ? ann[ann_r++ % ANN_N] : 0;
        pthread_mutex_unlock(&mx);
        if (!any) return;
        if (*fd < 0) *fd = dial(BEAST_IN);
        if (*fd < 0) { feed_failed = true; continue; }
        uint8_t f[40];
        int n = squitter_frame(a, f);
        if (send(*fd, f, n, MSG_NOSIGNAL) != n) { close(*fd); *fd = -1; feed_failed = true; continue; }
        feed_failed = false;
    }
}

static void *reader(void *arg)
{
    (void)arg;
    uint8_t buf[4096], f[32];
    int out = -1;
    while (!quit) {
        int fd = dial(port);
        if (fd < 0) {
            for (int k = 0; k < 20 && !quit; k++) usleep(100000);
            continue;
        }
        connected = true;
        int state = 0, need = 0, fn = 0, type = 0;
        bool esc = false;
        while (!quit) {
            announce_pending(&out);
            ssize_t n = recv(fd, buf, sizeof(buf), 0);
            if (n < 0 && (errno == EAGAIN || errno == EWOULDBLOCK || errno == EINTR)) continue;
            if (n <= 0) break;
            double wall = wall_now();
            for (ssize_t i = 0; i < n; i++) {
                uint8_t c = buf[i];
                if (state == 2) {
                    if (esc) {
                        esc = false;
                        if (c != 0x1a) state = 1;           /* lone 0x1a: a new frame starts, c is its type */
                    } else if (c == 0x1a) { esc = true; continue; }
                    if (state == 2) {
                        f[fn++] = c;
                        if (fn == need) {
                            if (type != '1') push(f, need - 7, wall);
                            state = 0;
                        }
                        continue;
                    }
                }
                if (state == 0) { if (c == 0x1a) state = 1; continue; }
                /* state 1: type byte; payload = 6 timestamp + 1 signal + message */
                type = c;
                need = c == '1' ? 9 : c == '2' ? 14 : c == '3' ? 21 : 0;
                state = need ? 2 : (c == 0x1a ? 1 : 0);
                fn = 0;
            }
        }
        close(fd);
        connected = false;
    }
    if (out >= 0) close(out);
    return NULL;
}

int rx1090_start(const char *hostport)
{
    const char *c = strrchr(hostport, ':');
    if (c) {
        snprintf(host, sizeof(host), "%.*s", (int)(c - hostport), hostport);
        snprintf(port, sizeof(port), "%s", c + 1);
    } else {
        snprintf(host, sizeof(host), "%s", hostport);
        snprintf(port, sizeof(port), "30005");
    }
    if (!host[0]) snprintf(host, sizeof(host), "localhost");
    if (!port[0] || atoi(port) <= 0) return -1;
    snprintf(S.peer, sizeof(S.peer), "%s:%s", host, port);
    const char *lp = getenv("RX1090_LOG");
    if (lp && (rlog = fopen(lp, "w")))
        fprintf(rlog, "# R r_s df addr code hex | C t_stream hw_off rw_off (wall - stream time minima) drift_ppm "
                      "offset | L t_stream offset kind drift_ppm (r = t + offset + delay)\n");
    quit = false;
    if (pthread_create(&thr, NULL, reader, NULL)) return -1;
    thr_on = true;
    return 0;
}

void rx1090_stop(void)
{
    if (rlog) { fclose(rlog); rlog = NULL; }
    if (!thr_on) return;
    quit = true;
    pthread_join(thr, NULL);
    thr_on = false;
}

/* ---------------- main thread ---------------- */

void rx1090_set_geometry(double rlat, double rlon, double ralt, double xlat, double xlon, double xalt)
{
    geo_ecef(rlat, rlon, ralt, g_radar);
    geo_ecef(xlat, xlon, xalt, g_rx);
    g_d0 = geo_dist(g_radar, g_rx);
    double la = rlat * M_PI / 180, lo = rlon * M_PI / 180, v[3];
    g_up[0] = cos(la) * cos(lo); g_up[1] = cos(la) * sin(lo); g_up[2] = sin(la);
    for (int i = 0; i < 3; i++) v[i] = g_rx[i] - g_radar[i];
    g_xz = v[0] * g_up[0] + v[1] * g_up[1] + v[2] * g_up[2];
    g_xh = sqrt(g_d0 * g_d0 - g_xz * g_xz);
    g_ralt = ralt;
    g_ref_lat = xlat;
    g_ref_lon = xlon;
    geo_on = true;
}

static double off_at(double t) { return c_off + c_k * (t - c_t0); }

/* match window around the predicted reply: wider the longer no match corrected the clock */
static double geo_tol(double t)
{
    if (lock_frames) return GEO_TOL;         /* clock from 1090 windows: only the position age counts */
    double w = GEO_TOL + TOL_GROW * fabs(t - last_match);
    return w < TOL_MAX ? w : TOL_MAX;
}

/* known aircraft: how far its delay may be from the last one, age s later */
static double ac_tol(double age)
{
    if (!lock_frames) return AC_TOL;         /* the clock itself may have moved this much */
    double w = AC_TOL0 + AC_MOVE * age;
    return w < AC_TOL ? w : AC_TOL;
}

static struct acd *ac_get(uint32_t a, bool create)
{
    uint32_t h = (a * 2654435761u) >> 20;
    for (int k = 0; k < AC_N; k++) {
        struct acd *e = &act[(h + k) % AC_N];
        if (!e->used) {
            if (!create) return NULL;
            e->used = true;
            e->addr = a;
            e->t = -1e300;
            return e;
        }
        if (e->addr == a) return e;
    }
    return NULL;
}

/* predicted reply delay for this aircraft at receiver time r: turnaround
 * plus radar -> aircraft -> receiver minus radar -> receiver */
static bool predict(const struct acd *a, double r, double *d)
{
    if (!geo_on || !a || !a->pos || fabs(r - a->pos_r) > POS_AGE) return false;
    *d = TURNAROUND + (geo_dist(g_radar, a->ecef) + geo_dist(a->ecef, g_rx) - g_d0) / C_LIGHT;
    return true;
}

/* first reply index (absolute) with r >= x; replies are in time order */
static uint64_t first_at(double x)
{
    uint64_t lo = rp_n > RP_N ? rp_n - RP_N : 0, hi = rp_n;
    while (lo < hi) {
        uint64_t m = lo + (hi - lo) / 2;
        if (rp[m % RP_N].r < x) lo = m + 1; else hi = m;
    }
    return lo;
}

static bool rollcall(int kind) { return kind == K_UF4 || kind == K_UF5 || kind == K_UF20 || kind == K_UF21; }

/* the reply format a transponder may send to this uplink format */
static bool df_fits(int kind, int df)
{
    if (kind == K_UF4 || kind == K_UF20) return df == 4 || df == 20;
    return df == 5 || df == 21;
}

/* can reply y answer query q? roll-calls by address and DF, UF11 by II/SI code */
static bool pair_ok(const struct qrec *q, const struct reply *y)
{
    if (q->kind == K_UF11) return y->df == 11 && y->code == q->code;
    return y->df != 11 && y->addr == q->addr && df_fits(q->kind, y->df);
}

static void unlock(void)
{
    locked = lock_geo = lock_frames = false;
    for (int k = 0; k < AC_N; k++) act[k].t = -1e300;     /* forget delays, keep announce times */
}

static void try_lock(double t_now)
{
    static int hist[LOCK_NB + 1];
    double hwo = minage_get(&hw), rwo = minage_get(&rw);
    if (hwo > 1e299 || rwo > 1e299 || !qh_n) return;
    double oc = hwo - rwo;                  /* r - t, coarse */
    memset(hist, 0, sizeof(hist));
    uint64_t q0 = qh_n > QH_N ? qh_n - QH_N : 0;
    int npairs = 0, best = 0;
    double omin = 1e300;
    for (int pass = 0; pass < 2; pass++) {
        for (uint64_t i = q0; i < qh_n; i++) {
            const struct qrec *q = &qh[i % QH_N];
            if (q->t < t_now - LOCK_HIST || q->t > t_now - 0.5 || q->kind == K_UF11) continue;
            for (uint64_t j = first_at(q->t + oc - LOCK_SPAN); j < rp_n; j++) {
                const struct reply *y = &rp[j % RP_N];
                if (y->r > q->t + oc + LOCK_SPAN) break;
                if (!pair_ok(q, y)) continue;
                double o = y->r - q->t;
                int b = (int)((o - (oc - LOCK_SPAN)) / LOCK_BIN);
                if (b < 0 || b >= LOCK_NB) continue;
                if (pass == 0) { hist[b]++; npairs++; }
                else if (b >= best && b < best + LOCK_W && o < omin) omin = o;
            }
        }
        if (pass == 0) {
            if (npairs < LOCK_MIN) {
                if (getenv("RX1090_DEBUG")) fprintf(stderr, "rx1090 lock: t=%.1f roll-call/reply pairs %d\n", t_now, npairs);
                return;
            }
            int sum = 0, c1 = 0, c2 = 0;
            for (int b = 0; b < LOCK_W; b++) sum += hist[b];
            static int win[LOCK_NB];
            for (int b = 0; b + LOCK_W <= LOCK_NB; b++) {
                win[b] = sum;
                if (sum > c1) { c1 = sum; best = b; }
                sum += hist[b + LOCK_W] - hist[b];
            }
            for (int b = 0; b + LOCK_W <= LOCK_NB; b++)
                if (abs(b - best) > 2 * LOCK_W && win[b] > c2) c2 = win[b];
            if (getenv("RX1090_DEBUG"))
                fprintf(stderr, "rx1090 lock: t=%.1f pairs %d, best 4 ms window %d, best other %d\n", t_now, npairs, c1, c2);
            if (c1 < LOCK_MIN || c1 < LOCK_RATIO * c2) return;
        } else {
            unlock();
            locked = true;
            c_off = omin - 0.000128;        /* nearest-path aircraft: delay ~ turnaround only */
            c_k = first_drift();
            c_t0 = t_now;
            misfits = 0;
            S.locks++;
            if (rlog) fprintf(rlog, "L %.3f %.7f hist %.3f\n", t_now, c_off, 1e6 * c_k);
            S.rc = S.rc_ans = S.ac = S.ac_ans = 0;
            if (getenv("RX1090_DEBUG")) {
                fprintf(stderr, "rx1090 LOCKED at t=%.1f: offset %.6f s; pairs in the window (query t, addr, delay us):\n",
                        t_now, c_off);
                for (uint64_t i = q0; i < qh_n; i++) {
                    const struct qrec *q = &qh[i % QH_N];
                    if (q->t < t_now - LOCK_HIST || q->t > t_now - 0.5 || q->kind == K_UF11) continue;
                    for (uint64_t j = first_at(q->t + c_off - 0.001); j < rp_n; j++) {
                        const struct reply *y = &rp[j % RP_N];
                        if (y->r > q->t + c_off + WIN_HI) break;
                        if (pair_ok(q, y))
                            fprintf(stderr, "  %.4f %06x DF%d %.1f\n", q->t, q->addr, y->df, 1e6 * (y->r - q->t - c_off));
                    }
                }
            }
        }
    }
}

/* number of scans among query times: a radar pairs roll-calls (UF4 + UF5) at a
 * fixed spacing, so replies to another radar's pair to the same aircraft can
 * agree by chance within one scan; agreement across scans or aircraft cannot */
static int count_scans(double *ts, int n)
{
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && ts[j] < ts[j - 1]; j--) { double x = ts[j]; ts[j] = ts[j - 1]; ts[j - 1] = x; }
    int c = n > 0;
    for (int i = 1; i < n; i++) if (ts[i] - ts[i - 1] > SCAN_GAP) c++;
    return c;
}

/* geometry lock: off = r - t - predicted delay for every pair whose aircraft
 * has a position, over the last LOCK_HIST s, each carried to t_now with the
 * drift k0; true pairs agree within ~10 us, cross pairs (other queries to the
 * same aircraft) are hundreds of us away. The winning 30 us window needs
 * INDEP_MIN independent sources (scans, plus aircraft beyond the first); a
 * first lock also takes PAIRS_ONE pairs from one source. recheck (already
 * locked on positions): re-locks only when no pair fits the tracked clock and
 * another offset has INDEP_MOVE sources */
static bool try_lock_geo(double t_now, bool recheck)
{
    double k0 = recheck ? c_k : first_drift();
    static int hist[GEO_NB + GEO_W];
    double hwo = minage_get(&hw), rwo = minage_get(&rw);
    if (hwo > 1e299 || rwo > 1e299 || !qh_n) return false;
    double oc = hwo - rwo, base = oc - LOCK_SPAN;
    memset(hist, 0, sizeof(hist));
    uint64_t q0 = qh_n > QH_N ? qh_n - QH_N : 0;
    int npairs = 0, best = 0, nsum = 0, ncur = 0, nwin = 0;
    double sum = 0, twin[256];
    uint32_t awin[256];
    for (int pass = 0; pass < 2; pass++) {
        for (uint64_t i = q0; i < qh_n; i++) {
            const struct qrec *q = &qh[i % QH_N];
            if (q->t < t_now - LOCK_HIST || q->t > t_now - 0.5) continue;
            for (uint64_t j = first_at(q->t + oc - LOCK_SPAN); j < rp_n; j++) {
                const struct reply *y = &rp[j % RP_N];
                if (y->r > q->t + oc + LOCK_SPAN) break;
                double d;
                if (!pair_ok(q, y) || !predict(ac_get(y->addr, false), y->r, &d)) continue;
                double o = y->r - q->t - d + k0 * (t_now - q->t);
                int b = (int)((o - base) / GEO_BIN);
                if (b < 0 || b >= GEO_NB) continue;
                if (pass == 0) {
                    hist[b]++; npairs++;
                    if (recheck && fabs(y->r - q->t - d - off_at(q->t)) < geo_tol(q->t)) ncur++;   /* fits the tracked clock */
                } else if (b >= best && b < best + GEO_W) {
                    sum += o; nsum++;
                    if (nwin < 256) { awin[nwin] = y->addr; twin[nwin++] = q->t; }
                    if (getenv("RX1090_DEBUG"))
                        fprintf(stderr, "  %.4f %-5s %06x DF%d delay %.1f us (predicted)\n", q->t, kind_name(q->kind),
                                y->addr, y->df, 1e6 * d);
                }
            }
        }
        if (pass == 1) break;
        if (npairs < 2) return false;
        int c1 = 0, c2 = 0, w = 0;
        for (int b = 0; b < GEO_W; b++) w += hist[b];
        for (int b = 0; b + GEO_W <= GEO_NB; b++) {
            if (w > c1) { c1 = w; best = b; }
            w += hist[b + GEO_W] - hist[b];
        }
        w = 0;
        for (int b = 0; b < GEO_W; b++) w += hist[b];
        for (int b = 0; b + GEO_W <= GEO_NB; b++) {
            if (abs(b - best) > 10 && w > c2) c2 = w;
            w += hist[b + GEO_W] - hist[b];
        }
        if (getenv("RX1090_DEBUG"))
            fprintf(stderr, "rx1090 geo lock: t=%.1f pairs with position %d, best 30 us window %d, best other %d\n",
                    t_now, npairs, c1, c2);
        if (c1 < 2 || c1 < 2 * c2) return false;
        if (recheck && (ncur || fabs(base + (best + 0.5 * GEO_W) * GEO_BIN - off_at(t_now)) < GEO_TOL)) return false;
    }
    int scans = count_scans(twin, nwin), nac = 0;
    for (int i = 0; i < nwin; i++) {
        int j = 0;
        while (j < i && awin[j] != awin[i]) j++;
        nac += j == i;
    }
    bool indep = scans + nac - 1 >= (recheck ? INDEP_MOVE : INDEP_MIN) || (!recheck && nsum >= PAIRS_ONE);
    if (getenv("RX1090_DEBUG"))
        fprintf(stderr, "rx1090 geo window: %d pairs from %d scans, %d aircraft, %+.3f ms from coarse%s\n", nsum, scans, nac,
                1e3 * (sum / nsum - oc), indep ? "" : " (too few)");
    if (!indep) return false;
    if (getenv("RX1090_DEBUG")) fprintf(stderr, "rx1090 GEO %s at t=%.1f\n", recheck ? "RE-LOCKED" : "LOCKED", t_now);
    unlock();
    if (recheck) S.relocks++;
    locked = lock_geo = true;
    c_off = sum / nsum;
    c_k = k0;
    c_t0 = last_match = t_now;
    misfits = 0;
    nres = 0;
    if (!recheck) {
        S.locks++;
        S.rc = S.rc_ans = S.ac = S.ac_ans = 0;
    }
    if (rlog) fprintf(rlog, "L %.3f %.7f %s %.3f\n", t_now, c_off, recheck ? "regeo" : "geo", 1e6 * c_k);
    if (getenv("RX1090_DEBUG")) fprintf(stderr, "rx1090 offset %.7f s from %d pairs, drift %+.2f ppm\n", c_off, nsum, 1e6 * c_k);
    return true;
}

void rx1090_clock_reset(void)
{
    unlock();
    hw = (struct minage){ 1e300, 1e300, 0 };
}

void rx1090_hclock(double t, double wall)
{
    minage_put(&hw, wall - t, wall);
}

void rx1090_note_query(const struct interrog *it)
{
    if (it->tcas) return;
    if (it->kind == K_UF11) {
        int code = modes_own_code(it);
        if (code > 0) qh[qh_n++ % QH_N] = (struct qrec){ it->t, NO_ADDR_RX, it->kind, code };
        return;
    }
    if (!rollcall(it->kind)) return;
    qh[qh_n++ % QH_N] = (struct qrec){ it->t, it->addr, it->kind, -1 };
    struct acd *a = ac_get(it->addr, true);
    double w = wall_now();
    if (a && w - a->announced > ANN_EVERY) {
        a->announced = w;
        pthread_mutex_lock(&mx);
        if (ann_w - ann_r < ANN_N) ann[ann_w++ % ANN_N] = it->addr;
        pthread_mutex_unlock(&mx);
    }
}

/* ---------------- sync from 1090 windows ---------------- */

void rx1090_hframe(const uint8_t *msg, int nbits, double t)
{
    if (hq_w - hq_r >= HQ_N) hq_r++;
    struct hframe *h = &hq[hq_w++ % HQ_N];
    h->t = t; h->len = (uint8_t)(nbits / 8);
    memcpy(h->msg, msg, h->len);
}

static uint64_t sf_first_at(double r)
{
    uint64_t lo = sf_n > SYNC_N ? sf_n - SYNC_N : 0, hi = sf_n;
    while (lo < hi) {
        uint64_t m = lo + (hi - lo) / 2;
        if (sf[m % SYNC_N].r < r) lo = m + 1; else hi = m;
    }
    return lo;
}

static int cmp_d(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* a window's frames are in: lock or update the clock from their median */
static void sync_window(double t_now)
{
    int n = wn;
    wn = 0;
    if (n < SYNC_MIN) {
        if (getenv("RX1090_DEBUG")) fprintf(stderr, "rx1090 sync: window with %d matched frames, need %d\n", n, SYNC_MIN);
        return;
    }
    double so[HQ_N];
    memcpy(so, wo, n * sizeof(double));
    qsort(so, n, sizeof(double), cmp_d);
    double med = so[n / 2], sum = 0, st = 0, dev = 0;
    int k = 0;
    for (int i = 0; i < n; i++)
        if (fabs(wo[i] - med) < SYNC_AGREE) { sum += wo[i]; st += wt[i]; k++; }
    if (k < SYNC_MIN) return;
    double o = sum / k, t = st / k;
    for (int i = 0; i < n; i++) if (fabs(wo[i] - med) < SYNC_AGREE) dev += (wo[i] - o) * (wo[i] - o);
    double jump = locked ? o - off_at(t) : NAN;
    if (!lock_frames || isnan(jump) || fabs(jump) > SYNC_BAD) { sync_iv = SYNC_IV_MIN; sync_good = 0; }
    else if (fabs(jump) < SYNC_GOOD && ++sync_good >= 2 && sync_iv < SYNC_IV_MAX) {
        sync_iv = 2 * sync_iv < SYNC_IV_MAX ? 2 * sync_iv : SYNC_IV_MAX;
    }
    if (!lock_frames) {
        double k0 = first_drift();
        unlock();
        locked = lock_frames = true;
        lock_geo = geo_on;
        c_k = k0;
        misfits = 0;
        nres = 0;
        S.locks++;
        S.rc = S.rc_ans = S.ac = S.ac_ans = 0;
    } else if (t - sync_t > 5) {
        c_k = (o - sync_o) / (t - sync_t);    /* ~0.2 us over 30 s: 0.01 ppm, no smoothing needed */
        if (c_k > 1e-4) c_k = 1e-4;
        if (c_k < -1e-4) c_k = -1e-4;
        k_learned = true;
    }
    c_off = o;
    c_t0 = last_match = t;
    sync_o = o; sync_t = t;
    S.sync_windows++;
    S.sync_n = k;
    S.sync_rms_us = 1e6 * sqrt(dev / k);
    S.sync_jump_us = 1e6 * jump;
    S.sync_iv = sync_iv;
    if (rlog) fprintf(rlog, "S %.3f %.7f %d %.2f %.2f %.4f %.0f\n", t, o, k, S.sync_rms_us, S.sync_jump_us, 1e6 * c_k, sync_iv);
    if (getenv("RX1090_DEBUG"))
        fprintf(stderr, "rx1090 sync: t=%.1f %d of %d frames, offset %.7f s, rms %.2f us, vs tracked %+.2f us, drift %+.3f ppm, "
                "next in %.0f s\n", t, k, n, o, S.sync_rms_us, S.sync_jump_us, 1e6 * c_k, sync_iv);
    (void)t_now;
}

/* HackRF window frames whose dump1090 copy has had time to arrive */
static void sync_poll(double t_now)
{
    double hwo = minage_get(&hw), rwo = minage_get(&rw);
    if (hwo > 1e299 || rwo > 1e299) return;
    double oc = hwo - rwo;
    while (hq_r < hq_w) {
        struct hframe *h = &hq[hq_r % HQ_N];
        if (last_r - oc < h->t + SYNC_SPAN && t_now - h->t < 30) break;   /* dump1090 not past it yet */
        hq_r++;
        if (wn && h->t - wt[wn - 1] > 5) sync_window(t_now);           /* a new window began */
        int same = 0;
        double r = 0;
        for (uint64_t j = sf_first_at(h->t + oc - SYNC_SPAN); j < sf_n; j++) {
            const struct sframe *f = &sf[j % SYNC_N];
            if (f->r > h->t + oc + SYNC_SPAN) break;
            if (f->len == h->len && !memcmp(f->msg, h->msg, h->len)) { same++; r = f->r; }
        }
        if (same == 1 && wn < HQ_N) { wo[wn] = r - h->t + SYNC_REF; wt[wn++] = h->t; }
    }
    if (wn && t_now - wt[wn - 1] > 5 && hq_r == hq_w) sync_window(t_now);
}

/* time the next sync window should follow the last one; a window that is
 * overdue (it found too few frames) is retried soon */
double rx1090_sync_interval(double t_now)
{
    if (!lock_frames || t_now - sync_t > sync_iv + 45) return SYNC_IV_MIN;
    return sync_iv;
}

void rx1090_poll(double t_now)
{
    t_stream = t_now;
    static struct raw tmp[INQ];
    pthread_mutex_lock(&mx);
    uint64_t n = in_w - in_r;
    for (uint64_t k = 0; k < n; k++) tmp[k] = inq[(in_r + k) % INQ];
    in_r = in_w;
    pthread_mutex_unlock(&mx);

    for (uint64_t k = 0; k < n; k++) {
        const struct raw *w = &tmp[k];
        if (!w->ts) continue;                           /* no receiver timestamp */
        double r = w->ts / 12e6;
        if (last_r >= 0 && r < last_r - 1.0) {          /* dump1090 restarted: new clock */
            unlock();
            rw = (struct minage){ 1e300, 1e300, 0 };
            rp_n = 0;
        }
        last_r = r;
        last_rx_wall = w->wall;
        minage_put(&rw, w->wall - r, w->wall);
        S.frames++;
        int nbits = w->len * 8, df = w->msg[0] >> 3;
        if (nbits != 56 && nbits != 112) continue;
        adsb_frame(w->msg, nbits);                      /* callsigns */
        uint32_t syn = modes_crc(w->msg, nbits) ^ getbits(w->msg, nbits - 23, nbits);
        double la, lo, alt;
        if (geo_on && (df == 17 || df == 18) && nbits == 112 && !syn &&
            df17_position(w->msg, g_ref_lat, g_ref_lon, &la, &lo, &alt)) {
            double e[3];
            geo_ecef(la, lo, alt, e);
            struct acd *a = ac_get(getbits(w->msg, 9, 32), true);
            if (a && geo_dist(e, g_rx) < 400e3) {
                if (!a->pos) S.npos++;
                memcpy(a->ecef, e, sizeof(e)); a->pos = true; a->pos_r = r;
                a->h = alt - g_ralt; a->h_t = r;
            }
        }
        if (((df == 17 || df == 18) && nbits == 112 && !syn) || (df == 11 && nbits == 56 && !(syn & ~0x7Fu))) {
            struct sframe *f = &sf[sf_n++ % SYNC_N];
            f->r = r; f->len = w->len;
            memcpy(f->msg, w->msg, w->len);
        }
        struct reply y = { .r = r, .df = (uint8_t)df, .code = -1, .addr = NO_ADDR_RX, .len = w->len };
        memcpy(y.msg, w->msg, w->len);
        if (df == 11 && nbits == 56) {
            if (syn & ~0x7Fu) continue;
            int cl = syn >> 4 & 7, ic = syn & 15;
            if (syn) y.code = cl == 0 ? ic : cl <= 4 ? 16 + ic + 16 * (cl - 1) : -1;
            y.addr = getbits(w->msg, 9, 32);
        } else if ((df == 4 || df == 5) && nbits == 56) y.addr = syn;
        else if ((df == 20 || df == 21) && nbits == 112) y.addr = syn;
        else continue;
        rp[rp_n++ % RP_N] = y;
        if (rlog) {
            fprintf(rlog, "R %.7f %d %06x %d ", r, df, y.addr & 0xFFFFFF, y.code);
            for (int b = 0; b < w->len; b++) fprintf(rlog, "%02x", w->msg[b]);
            fputc('\n', rlog);
        }
    }
    static double last_c;
    if (rlog && t_now - last_c >= 10) {
        last_c = t_now;
        fprintf(rlog, "C %.3f %.6f %.6f %.3f %.7f\n", t_now, minage_get(&hw), minage_get(&rw),
                locked ? 1e6 * c_k : NAN, locked ? off_at(t_now) : NAN);
        fflush(rlog);
    }

    double wall = wall_now();
    if (wall - rate_t0 >= 5) {
        S.rate = rate_t0 > 0 ? (S.frames - rate_n0) / (wall - rate_t0) : 0;
        rate_t0 = wall;
        rate_n0 = S.frames;
    }
    sync_poll(t_now);
    if (lock_frames && t_now - sync_t > SYNC_STALE) {
        lock_frames = false;                /* windows stopped: lock from replies again */
        if (getenv("RX1090_DEBUG")) fprintf(stderr, "rx1090 sync: no 1090 window for %.0f s\n", SYNC_STALE);
    }
    if (lock_frames) return;
    if (locked && misfits >= UNLOCK_N) unlock();
    if (t_now - last_try > 1.0 && (!locked || (geo_on && !lock_geo))) {
        last_try = t_now;
        if (geo_on) try_lock_geo(t_now, false);   /* positions known: the histogram lock only adds false matches */
        else if (!locked) try_lock(t_now);
    }
    if (lock_geo && t_now - last_check > RECHECK) {
        last_check = t_now;
        try_lock_geo(t_now, true);
    }
}

bool rx1090_ready(double t)
{
    if (!connected) return true;
    if (!locked) return t < t_stream - HOLD_UNLOCKED;
    if (wall_now() - last_rx_wall > 1.0) return true;   /* feed stalled: do not hold the analysis */
    return last_r - off_at(t) > t + WIN_HI + 0.001;
}

static void clock_update(double t, double e)
{
    if (lock_frames) return;                 /* the 1090 windows are far more precise */
    double dt = t - c_t0;
    c_off = off_at(t) + (dt > 30 ? 1.0 : 0.3) * e;     /* after a long gap the drift error dominates */
    if (dt > 1.0) { c_k += (dt < 60 ? dt / 200 : 0.3) * e / dt; k_learned |= lock_geo; }   /* close matches: small steps, noise */
    last_match = t;
    if (c_k > 1e-4) c_k = 1e-4;
    if (c_k < -1e-4) c_k = -1e-4;
    c_t0 = t;
}

/* keep a matched reply for display; the query points at its replies */
static void keep(struct interrog *x, const struct reply *y, double d, float res_us, int req)
{
    struct rx1090_rec *k = &recs[nrec++ % REC_N];
    k->df = y->df;
    k->len = y->len;
    memcpy(k->msg, y->msg, y->len);
    k->delay_us = (float)(d * 1e6);
    k->res_us = res_us;
    k->req_bds = req;
    k->excess_m = lock_geo ? (float)((d - TURNAROUND) * C_LIGHT) : NAN;
    k->h_m = k->adsb_m = NAN;
    struct acd *a = ac_get(y->df == 11 ? y->addr : x->addr, false);
    int ft = reply_alt_ft(y->msg, y->df);
    if (ft != INT32_MIN) {
        k->h_m = (float)(ft * 0.3048 - g_ralt);
        if (a) { a->h = k->h_m; a->h_t = y->r; }
    } else if (a && a->h_t > 0 && fabs(y->r - a->h_t) < AC_AGE) k->h_m = (float)a->h;
    if (a && geo_on && a->pos && fabs(y->r - a->pos_r) <= POS_AGE) k->adsb_m = (float)geo_dist(g_radar, a->ecef);
    if (x->rep_n < 255) x->rep_n++;
}

/* Radar->aircraft range from the reply's path excess and the antenna angle az
 * (from the radar->antenna direction; its sign does not matter). Radar frame:
 * x towards the antenna, z up; the aircraft at ground range g sits at height
 * h - g^2/2R on the curved earth. Solve |A| + |A - X| = |X| + excess for g:
 * near the radar->antenna line there can be two solutions. */
int rx1090_range(const struct rx1090_rec *r, double az, double out_m[2])
{
    if (!geo_on || isnan(r->excess_m)) return 0;
    double h = isnan(r->h_m) ? 0 : r->h_m, k = g_d0 + r->excess_m, c = cos(az), sn = sin(az);
    int n = 0;
#define RANGE_F(g, f) do { double ax = (g) * c, ay = (g) * sn, z = h - (g) * (g) / (2 * 6371e3); \
        f = sqrt(ax * ax + ay * ay + z * z) + sqrt((ax - g_xh) * (ax - g_xh) + ay * ay + (z - g_xz) * (z - g_xz)) - k; } while (0)
    double g0 = 0, f0;
    RANGE_F(g0, f0);
    for (double g1 = 200; g1 <= 600e3 && n < 2; g1 += 200) {
        double f1;
        RANGE_F(g1, f1);
        if ((f0 < 0) != (f1 < 0)) {
            double lo = g0, hi = g1, flo = f0;
            for (int i = 0; i < 30; i++) {
                double m = 0.5 * (lo + hi), fm;
                RANGE_F(m, fm);
                if ((fm < 0) == (flo < 0)) { lo = m; flo = fm; } else hi = m;
            }
            double g = 0.5 * (lo + hi), z = h - g * g / (2 * 6371e3);
            out_m[n++] = sqrt(g * g + z * z);
        }
        g0 = g1; f0 = f1;
    }
#undef RANGE_F
    return n;
}

const struct rx1090_rec *rx1090_rec_get(uint64_t idx)
{
    return idx < nrec && nrec - idx <= REC_N ? &recs[idx % REC_N] : NULL;
}

void rx1090_match(struct interrog *x)
{
    x->answered = false;
    x->reply_df = 0;
    x->nreplies = 0;
    x->reply_us = NAN;
    x->reply_res_us = NAN;
    x->checked = false;
    x->rep_n = 0;
    if (!locked || x->tcas) return;
    bool rc = rollcall(x->kind);
    int code = x->kind == K_UF11 ? modes_own_code(x) : -1;
    if (!rc && code <= 0) return;           /* code 0 = II0: indistinguishable from squitters */
    x->checked = true;
    x->rep_idx = nrec;
    int req = -1, b1, b2;
    if (rc && modes_bds(x->msg, &b1, &b2)) req = b1 * 16 + b2;

    double off = off_at(x->t), lo = x->t + off + WIN_LO, hi = x->t + off + WIN_HI;
    uint64_t j = first_at(lo);

    if (rc) {
        S.rc++;
        struct acd *a = ac_get(x->addr, false);
        double pd = 0;
        bool geo = lock_geo && predict(a, x->t + off, &pd);
        bool known = a && x->t - a->t < AC_AGE;
        double tol = geo ? geo_tol(x->t) : known ? ac_tol(x->t - a->t) : 0;
        /* with a position: the reply at the predicted delay, nothing else.
         * Known aircraft: the reply at its usual delay, nothing else (other
         * radars' queries draw replies at other delays). Otherwise the
         * earliest one, and learn the delay */
        struct reply *best = NULL, *first = NULL;
        double bd = 0, fd = 0;
        for (; j < rp_n; j++) {
            struct reply *y = &rp[j % RP_N];
            if (y->r > hi) break;
            if (y->used || y->addr != x->addr || !df_fits(x->kind, y->df)) continue;
            double d = y->r - x->t - off;
            if (lock_frames && d < TURNAROUND - 5e-6) continue;     /* exact clock: nothing before the turnaround */
            if (!first) { first = y; fd = d; }
            if (geo) {
                if (fabs(d - pd) < tol && (!best || fabs(d - pd) < fabs(bd - pd))) { best = y; bd = d; }
            } else if (known && fabs(d - a->delay) < tol && (!best || fabs(d - a->delay) < fabs(bd - a->delay))) { best = y; bd = d; }
        }
        if (geo) {
            if (!best) { if (first) misfits++; return; }   /* replies there, none at the predicted time */
            clock_update(x->t, bd - pd);
            misfits = 0;
            x->reply_res_us = (float)((bd - pd) * 1e6);
            res[nres++ % RES_N] = bd - pd;
        } else if (best) { clock_update(x->t, bd - a->delay); misfits = 0; }
        else if (known) {                   /* replies there, none at this aircraft's delay */
            if (first) misfits++;
            return;
        } else if (first) { best = first; bd = fd; }
        else return;                        /* not heard: says nothing about the lock */
        best->used = true;
        x->answered = true;
        x->reply_df = best->df;
        x->reply_us = (float)(bd * 1e6);
        keep(x, best, bd, x->reply_res_us, req);
        S.rc_ans++;
        if (!a) a = ac_get(x->addr, true);
        if (a) { a->delay = best->r - x->t - off_at(x->t); a->t = x->t; }
    } else {
        S.ac++;
        for (; j < rp_n; j++) {
            struct reply *y = &rp[j % RP_N];
            if (y->r > hi) break;
            if (y->used || y->df != 11 || y->code != code) continue;
            double d = y->r - x->t - off, pd;
            struct acd *a = ac_get(y->addr, false);
            bool pos = lock_geo && predict(a, y->r, &pd);
            /* II/SI codes are reused by distant radars: with a position, only a
             * reply at the predicted time answers this radar's all-call; from a
             * known aircraft, only one at its usual delay */
            if (pos && fabs(d - pd) >= geo_tol(x->t)) continue;
            if (!pos && a && x->t - a->t < AC_AGE && fabs(d - a->delay) >= ac_tol(x->t - a->t)) continue;
            y->used = true;
            if (!x->nreplies || d * 1e6 < x->reply_us) x->reply_us = (float)(d * 1e6);
            if (pos) {
                if (isnan(x->reply_res_us)) x->reply_res_us = (float)((d - pd) * 1e6);
                clock_update(x->t, d - pd);
                res[nres++ % RES_N] = d - pd;
            }
            keep(x, y, d, pos ? (float)((d - pd) * 1e6) : NAN, -1);
            if (x->nreplies < 255) x->nreplies++;
        }
        if (x->nreplies) { x->answered = true; x->reply_df = 11; S.ac_ans++; }
    }
}

void rx1090_get_stats(struct rx1090_stats *s)
{
    S.connected = connected;
    S.feed_failed = feed_failed;
    S.locked = locked;
    S.geo_on = geo_on;
    S.frames_lock = lock_frames;
    S.off_ms = locked ? 1e3 * off_at(c_t0) : 0;
    S.drift_ppm = locked ? 1e6 * c_k : 0;
    S.geo = lock_geo;
    long n = nres < RES_N ? nres : RES_N;
    double ss = 0;
    for (long k = 0; k < n; k++) ss += res[k] * res[k];
    S.res_rms_us = n ? 1e6 * sqrt(ss / n) : NAN;
    S.nres = nres;
    *s = S;
}
