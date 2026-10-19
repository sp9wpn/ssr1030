/* Analysis state shared between analysis.c and display.c. */
#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "ssr1030.h"

#define BURST_HIST   64
#define SEQ_MAX      2048
#define MAX_RADARS   8
#define ADDR_HASH    8192
#define RECENT_N     8
#define NCODES       80          /* radar codes: II 0-15, SI 16+1..63 */
#define NO_TARGET    255
#define TR_N         1024        /* distinct token transitions tracked */

/* query-type groups used in per-burst tables; every interrogation is in exactly one */
enum group { G_AC, G_ACONLY, G_ACS, G_SALL, G_SALT, G_SSQK, G_CAL, G_SETALT, G_TRK, G_HDGSPD, G_OBDS,
             G_MSG, G_SOTHER, G_COUNT };

/* a committed burst, built from the interrogations that passed the filters */
struct burst {
    int id;
    int code;               /* radar code, -1 unknown */
    char key[12];           /* "SI17", "II3", "?" */
    bool partial;           /* touched an ADS-B gap: excluded from per-burst stats */
    double t0, t1, tpeak;
    float peak;
    int n;
    int cnt[K_COUNT];
    int gcnt[G_COUNT];      /* by query type (tok_group) */
    uint16_t ks[SEQ_MAX];   /* sequence tokens (seq_token) in arrival order */
    uint8_t ka[SEQ_MAX];    /* target index into addrs[], NO_TARGET when none */
    int nks;
    int naddr;              /* distinct roll-call addresses */
    uint32_t addrs[64];
    double period;          /* time since previous burst of the same radar */
};

/* burst detector state (all interrogations, before filtering) */
struct open_burst {
    bool open;
    int id;
    double t0, t1;
    int n_raw;
    int votes[NCODES];
    bool partial;
};

struct radar {
    int code;
    char key[12];
    int nbursts;             /* complete bursts */
    int npartial;
    double last_peak;
    bool gap_since;          /* an ADS-B gap happened since last_peak */
    double periods[64]; int np;
    double dwell_sum;
    long cnt_sum[G_COUNT];
    int cnt_min[G_COUNT], cnt_max[G_COUNT];
    long total_sum;
};

struct trans {
    uint16_t a, b;          /* sequence tokens */
    long n;
    double sum, min, max;
};

struct addr_entry {
    uint32_t addr;
    int n;
    double last;
    bool tcas;
};

struct analysis_state {
    struct analysis_cfg cfg;
    double now;
    double floor_db;           /* detection threshold dBFS (empty-bin value) */
    float bg_db;               /* background (sidelobe) level */

    long total[K_COUNT], inbeam[K_COUNT], side[K_COUNT];
    long tcas_total, tcas_uf0, tcas_uf16, tcas_ws;
    long n_all, n_pass;
    long filt_radar, filt_modes, filt_rollcall, filt_aircraft;

    struct open_burst ob;
    struct burst hist[BURST_HIST];
    int nhist;                  /* total committed bursts */
    double last_burst_end;
    struct radar radars[MAX_RADARS];
    int nradars;

    struct trans tr[TR_N];      /* in-beam token -> token transitions */
    int ntr;
    int prev_tok, prev_burst; double prev_t;

    long uf11_code[2][64][16];  /* [is_si][code][pr] */
    long bds[256];
    long pc[8];
    long di[8];
    struct addr_entry addrs[ADDR_HASH];
    int naddr_radar, naddr_tcas;

    int64_t target;             /* resolved --aircraft address */
    double gap_end;

    struct interrog recent[RECENT_N];
    int nrecent;
};

extern struct analysis_state AS;
int tok_group(uint16_t tok);     /* seq_token -> group */
const char *group_label(int g);
double median_d(const double *v, int n);
int analysis_top_transitions(struct trans *out, int max);
void burst_ranks(const struct burst *b, int *rank);          /* alphabetical target numbers, 1-based */
int burst_item(const struct burst *b, int s, const int *rank, char *out, size_t n);  /* run at s -> next */
void burst_targets(const struct burst *b, const int *rank, char *out, size_t n);   /* "1/DLH123 2/..." */

#endif
