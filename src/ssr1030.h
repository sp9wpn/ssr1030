/* ssr1030 - SSR / Mode S 1030 MHz interrogation analyser for HackRF.
 * Detection approach derived from dump1030 (Fuziih/dump1030), extended with
 * pulse-timing classification and full Mode S uplink DPSK decoding. */
#ifndef SSR1030_H
#define SSR1030_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>
#include <complex.h>

#define MODES_MAX_BITS 112

/* Interrogation kinds. Order matters for the dashboard. */
enum kind {
    K_MODE1, K_MODE2, K_MODEA, K_MODEB, K_MODEC, K_MODED,
    K_ACONLY_A, K_ACONLY_C,           /* intermode, short P4 (0.8 us) */
    K_ACS_A, K_ACS_C,                 /* intermode, long P4 (1.6 us)  */
    K_INTERMODE_OTHER,                /* mode 1/2/B/D with P4          */
    K_UF0, K_UF4, K_UF5, K_UF11, K_UF16, K_UF19, K_UF20, K_UF21, K_UF24,
    K_UF_OTHER,
    K_S_NOSYNC,                       /* Mode S preamble+P6, sync masked by P5 (sidelobe) */
    K_S_GARBLE,                       /* Mode S preamble+P6, no sync reversal and no P5 */
    K_COUNT
};

struct pulse {
    double lead;        /* leading edge (half amplitude), global sample index */
    double width;       /* microseconds */
    float amp;          /* linear magnitude (full scale ~ 127) */
    int density;        /* strong pulses in the preceding 1 ms */
    bool used;
    bool strong;        /* from the main threshold (else only usable for Mode S) */
};

struct interrog {
    uint64_t seq;
    double t;           /* seconds since start of stream (P1 leading edge) */
    int kind;
    float amp_db;       /* query amplitude, dBFS: P6 for Mode S, P1 for A/C */
    float p1_db;        /* P1 amplitude, dBFS */
    float p6_db;        /* Mode S P6 plateau, dBFS; NAN for A/C */
    float p2_db;        /* P2 - P1 in dB; NAN when no P2 */
    float p5_db;        /* Mode S: P5 - P6 in dB (>= 0: sidelobe); NAN when no P5 detected */
    float p4_width;     /* us, 0 when none */
    bool s1;            /* whisper-shout S1 pulse 2 us before P1 */
    bool sidelobe;      /* P2 >= P1 (A/C) or no sync reversal (S) */
    bool tcas;          /* UF0/UF16/whisper-shout: counted, excluded from analysis */
    int burst;          /* burst id, -1 when not in a main-beam burst */
    int radar;          /* attributed radar code (see radar_code), -1 unknown */
    float gain_db;      /* nominal receiver gain when received (amp + LNA + VGA) */
    /* Mode S */
    int nbits;
    uint8_t msg[14];
    int uf;
    uint32_t addr;      /* address recovered from AP (FFFFFF for all-calls) */
    float sync;         /* sync reversal correlation, -1 = perfect */
    float quality;      /* weakest bit confidence, 0..1 */
    bool len_mismatch;  /* P6 width disagrees with UF */
    /* 1090 MHz reply (--rx1090) */
    bool checked;       /* a reply could be matched: clock locked, roll-call or UF11 */
    bool answered;      /* a reply was matched: roll-call by address, UF11 by II/SI */
    uint8_t reply_df;   /* DF of the (first) reply, 0 none */
    uint8_t nreplies;   /* UF11: DF11 replies matched */
    float reply_us;     /* reply delay after the query, us (128 us turnaround + extra path) */
    float reply_res_us; /* measured - predicted from the aircraft's position (geometry lock), NAN unknown */
    uint64_t rep_idx;   /* first matched reply in the rx1090 record store, rep_n of them */
    uint8_t rep_n;
};

/* ---- demod.c ---- */
typedef void (*interrog_cb)(struct interrog *it, void *ctx);

struct demod;
struct demod *demod_new(double fs, float thr_db, float thr_abs, float modes_db, int dens_max, interrog_cb cb, void *ctx);
void demod_process(struct demod *d, const int8_t *iq, size_t nsamples);
void demod_flush(struct demod *d);
void demod_skip(struct demod *d, size_t nsamples);   /* advance time without data (retune gaps) */
struct demod_stats {
    uint64_t samples;
    uint64_t clipped;
    uint64_t pulses;
    uint64_t stray_pulses;
    uint64_t garbled;      /* A/C candidates rejected in dense interference */
    float noise_dbfs;
    float thr_dbfs;
    float peak_dbfs;       /* peak since last read */
};
void demod_get_stats(struct demod *d, struct demod_stats *s, bool reset_peak);

/* ---- modes.c ---- */
uint32_t modes_crc(const uint8_t *msg, int nbits);       /* parity over nbits-24 data bits */
uint32_t modes_ap_to_addr(const uint8_t *msg, int nbits);
int modes_uf(const uint8_t *msg);
int modes_uf_bits(int uf);
unsigned getbits(const uint8_t *msg, int first, int last); /* 1-based, inclusive */
void modes_describe(const struct interrog *it, char *buf, size_t len);
int modes_bds(const uint8_t *msg, int *bds1, int *bds2);  /* returns 1 when Comm-B requested */
int modes_iisi(const uint8_t *msg, int *is_si);          /* UF11 II/SI code */
int modes_own_code(const struct interrog *it);           /* radar code carried by the message, -1 none */

/* Radar codes: II n -> n (0-15), SI n -> 16 + n (SI 1-63). */
int radar_code_parse(const char *s);                     /* "SI17", "II3"; -1 on error */
void radar_code_str(int code, char *buf, size_t len);    /* "SI17", "II3", "?" */

const char *kind_label(int k);
const char *kind_name(int k);    /* identifier for log files */
const char *kind_token(int k);
const char *bds_abbrev(int bds);                 /* shortcode for BDS1*16+BDS2, NULL if unnamed */
uint16_t seq_token(const struct interrog *it);   /* kind + requested BDS, for burst sequences */
const char *seq_token_str(uint16_t t, char *buf, size_t len);
bool kind_is_modes(int k);
bool kind_is_allcall(int k);

/* ---- analysis.c ---- */
#define AC_MAX 16
struct ac_sel { int64_t addr; char cs[9]; };   /* hex address (addr >= 0) or callsign */

struct analysis_cfg {
    double burst_gap;      /* s: silence that closes a burst */
    int burst_min;         /* minimum interrogations for a burst */
    float member_db;       /* in-beam: at least this far above background */
    float core_db;         /* a burst needs one interrogation this far above background */
    /* filters */
    int radar_sel;         /* radar code to keep, -1 = all */
    bool modes_only;       /* drop Mode A/C and intermode interrogations */
    bool ac_only;          /* keep only plain Mode A/C (no P4), as sent by older SSR radars */
    bool rollcall_only;    /* keep only interrogations addressed to an aircraft */
    int nac;               /* aircraft to keep (0 = all) */
    struct ac_sel ac[AC_MAX];
};
void analysis_init(const struct analysis_cfg *cfg);
void analysis_set_floor(double thr_dbfs);
void analysis_set_sink(void (*sink)(const struct interrog *it));
void analysis_set_rx(bool (*ready)(double t), void (*match)(struct interrog *x));   /* --rx1090 */
void analysis_add(struct interrog *it);
void analysis_tick(double now);
void analysis_print_summary(FILE *f);
void analysis_gap(double t0, double t1);    /* no 1030 data between t0 and t1 (ADS-B window) */
double analysis_rotation(void);             /* rotation of the selected / busiest radar, 0 unknown */
double analysis_rotation_of(const char *key);  /* rotation of radar "SI17", "?" ..., 0 unknown */
double analysis_last_burst_end(void);       /* stream time of last committed burst end, -1 none */
int analysis_main_code(void);               /* selected radar, else the one with most bursts; -2 none */
void analysis_gain_shift(float delta_db);   /* receiver gain changed: shift background levels */
void analysis_set_gain(float gain_db);      /* nominal gain, stamped on committed bursts */
bool analysis_addr_known(uint32_t addr);    /* address was roll-called */
int ac_list_parse(const char *s, struct ac_sel *out, char *bad, size_t badlen);  /* count, -1 on error */
void ac_list_str(const struct analysis_cfg *c, bool resolved, char *buf, size_t n); /* "48AF01,LOT3KM=4A2182" */

/* ---- adsb.c ---- */
struct adsb_stats {
    long windows, frames, df17, df11, idents, commb_idents;
    double seconds;
};
void adsb_init(double fs);
void adsb_window_start(void);
void adsb_process(const int8_t *iq, size_t nsamples);
void adsb_get_stats(struct adsb_stats *s);
const char *cs_lookup(uint32_t addr);       /* callsign or NULL */
const char *ac_name(uint32_t addr);         /* callsign when known, else hex */
int64_t cs_find(const char *callsign);      /* address or -1 */
int cs_count(void);
unsigned cs_version(void);
int cs_load(const char *path);
int cs_save(const char *path);
int cs_fetch_adsblol(double lat, double lon, int radius_nm);   /* startup bootstrap, -1 on error */
void adsb_frame(const uint8_t *msg, int nbits);   /* decoded 1090 frame from elsewhere (--rx1090) */

/* ---- geo.c ---- */
bool df17_position(const uint8_t *msg, double ref_lat, double ref_lon, double *lat, double *lon, double *alt_m);
void geo_ecef(double lat, double lon, double alt_m, double xyz[3]);
double geo_dist(const double a[3], const double b[3]);
int geo_parse(const char *s, double *lat, double *lon, double *alt_m);   /* "LAT,LON[,ALT_M]", -1 on error */

/* ---- rx1090.c ---- */
struct rx1090_rec {         /* a reply matched to a query, kept for display */
    uint8_t df, len;        /* len: bytes */
    uint8_t msg[14];
    float delay_us, res_us; /* res_us NAN without a position */
    int req_bds;            /* register the query asked for (BDS1*16+BDS2), -1 none */
    float excess_m;         /* radar->aircraft->antenna minus radar->antenna, NAN without a position lock */
    float h_m;              /* aircraft height above the radar (reply or recent altitude), NAN unknown */
    float adsb_m;           /* radar->aircraft from a recent ADS-B position, NAN none */
};
const struct rx1090_rec *rx1090_rec_get(uint64_t idx);   /* NULL once overwritten */
void reply_describe(const struct rx1090_rec *r, char *buf, size_t n);   /* modes.c: decoded content */
int reply_alt_ft(const uint8_t *msg, int df);            /* modes.c: DF4/20 altitude, INT32_MIN none */
int rx1090_range(const struct rx1090_rec *r, double az, double out_m[2]);   /* radar->aircraft for antenna
                            angle az (rad) from the antenna direction: 0-2 solutions */
struct rx1090_stats {
    bool connected, locked;
    bool feeding;           /* roll-called addresses are announced to dump1090 */
    bool geo;               /* locked from aircraft positions */
    bool geo_on;            /* radar and receiver positions given */
    double res_rms_us;      /* rms of measured - predicted delay, last 64 matches */
    long nres;
    long npos;              /* aircraft with an ADS-B position */
    long frames;            /* Mode S frames received */
    double rate;            /* frames/s over the last few seconds */
    double off_ms, drift_ppm;
    long locks;             /* times the clock lock was found */
    long relocks;           /* times the geometry re-check moved the lock */
    long rc, rc_ans;        /* roll-calls checked / answered (since lock) */
    long ac, ac_ans;        /* UF11 all-calls checked / answered */
    char peer[128];
};
int rx1090_start(const char *hostport);              /* starts the reader thread, -1 on bad address */
void rx1090_stop(void);
void rx1090_hclock(double t, double wall);           /* HackRF stream time t arrived at host time wall */
void rx1090_clock_reset(void);                       /* HackRF samples were lost: lock again */
void rx1090_note_query(const struct interrog *it);   /* every demodulated interrogation, in time order */
void rx1090_poll(double t_now);                      /* main loop: take new replies, try to lock */
bool rx1090_ready(double t);                         /* replies for a query at t are in */
void rx1090_match(struct interrog *x);               /* set answered/reply fields, queries in time order */
void rx1090_get_stats(struct rx1090_stats *s);
void rx1090_set_geometry(double rlat, double rlon, double ralt, double xlat, double xlon, double xalt);

/* ---- display.c ---- */
struct run_info {
    const char *source;
    double fs;
    double freq;
    int lna, vga, amp;
    double bbw;            /* baseband filter, Hz */
    char hw[160];          /* board, firmware, serial */
    bool live;
    bool adsb;             /* ADS-B windows enabled */
    bool adsb_now;         /* currently listening on 1090 */
    bool rx1090;           /* replies from dump1090 (--rx1090) */
    double adsb_next;      /* seconds until the next window (approx.), <0 waiting for burst */
    char filters[256];
    bool agc;              /* automatic gain */
    bool graph;            /* burst graph screen instead of the dashboard */
    bool frozen;           /* space: screen frozen, data collection continues */
    bool autofreeze;       /* F: freeze on the first burst with matched replies */
    bool cmd_mode;         /* filter input line open */
    char cmd_kind;         /* 'r' radar or 'a' aircraft */
    char cmd[160];
    char cmd_msg[160];     /* result of the last command */
    double gain_age;       /* seconds since the last gain change, <0 never */
    int gain_delta;        /* last change, dB */
    char gain_why[64];
};
void display_init(void);
void display_restore(void);
void display_render(const struct run_info *ri, double now, struct demod *d, uint64_t drops);
void display_graph_key(int key);
void display_graph_pin(int burst_id);   /* show this burst in the graph (autofreeze) */
void display_status_only(const struct run_info *ri);   /* frozen dashboard: redraw the bottom row only */   /* 'U' 'D' 'P' 'N' 'H' 'E': scroll the burst graph */

/* ---- output.c ---- */
int output_open(const char *path, char sep);
void output_write(const struct interrog *it, double start_epoch);
void output_close(void);
void print_interrog_line(FILE *f, const struct interrog *it);

#endif
