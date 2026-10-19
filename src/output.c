/* Per-interrogation text output: scrolling list and CSV / space-delimited dump. */
#include <math.h>
#include <time.h>
#include "ssr1030.h"

static FILE *out;
static char osep = ',';

int output_open(const char *path, char sep)
{
    out = fopen(path, "w");
    if (!out) return -1;
    osep = sep;
    const char *cols[] = { "t_s", "utc", "kind", "uf", "amp_dbfs", "p2_rel_db", "p4_us",
                           "s1", "sidelobe", "tcas", "burst", "radar", "addr", "callsign", "ii_si", "pr", "rr",
                           "di", "bds", "sync", "quality", "hex" };
    for (size_t k = 0; k < sizeof(cols) / sizeof(cols[0]); k++)
        fprintf(out, "%s%s", k ? (char[]){ osep, 0 } : "", cols[k]);
    fputc('\n', out);
    return 0;
}

void output_write(const struct interrog *it, double start_epoch)
{
    if (!out) return;
    char s = osep;
    double ut = start_epoch + it->t;
    time_t sec = (time_t)ut;
    struct tm tm;
    gmtime_r(&sec, &tm);
    char ts[40];
    strftime(ts, sizeof(ts), "%Y-%m-%dT%H:%M:%S", &tm);

    char iisi[8] = "-", pr[4] = "-", rr[4] = "-", di[4] = "-", bds[8] = "-", addr[8] = "-", hex[32] = "-";
    char radar[8], cs[12] = "-";
    radar_code_str(it->radar, radar, sizeof(radar));
    if (radar[0] == '?') snprintf(radar, sizeof(radar), "-");
    if (kind_is_modes(it->kind) && it->kind != K_S_NOSYNC) {
        snprintf(addr, sizeof(addr), "%06X", it->addr);
        const char *c = cs_lookup(it->addr);
        if (c) snprintf(cs, sizeof(cs), "%s", c);
        int o = 0;
        for (int k = 0; k < it->nbits / 8; k++) o += snprintf(hex + o, sizeof(hex) - o, "%02X", it->msg[k]);
        if (it->uf == 11) {
            int si, c = modes_iisi(it->msg, &si);
            snprintf(iisi, sizeof(iisi), "%s%d", si ? "SI" : "II", c);
            snprintf(pr, sizeof(pr), "%u", getbits(it->msg, 6, 9));
        }
        if (it->uf == 4 || it->uf == 5 || it->uf == 20 || it->uf == 21) {
            snprintf(rr, sizeof(rr), "%u", getbits(it->msg, 9, 13));
            snprintf(di, sizeof(di), "%u", getbits(it->msg, 14, 16));
            int b1, b2;
            if (modes_bds(it->msg, &b1, &b2)) snprintf(bds, sizeof(bds), "%X%X", b1, b2);
        }
    }
    fprintf(out, "%.7f%c%s.%06dZ%c%s%c%d%c%.1f%c", it->t, s, ts, (int)((ut - sec) * 1e6), s,
            kind_name(it->kind), s, kind_is_modes(it->kind) ? it->uf : -1, s, it->amp_db, s);
    if (isnan(it->p2_db)) fprintf(out, "-%c", s); else fprintf(out, "%.1f%c", it->p2_db, s);
    fprintf(out, "%.2f%c%d%c%d%c%d%c%d%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%s%c%.2f%c%.2f%c%s\n",
            it->p4_width, s, it->s1, s, it->sidelobe, s, it->tcas, s, it->burst, s, radar, s, addr, s,
            cs, s, iisi, s, pr, s, rr, s, di, s, bds, s, it->sync, s, it->quality, s, hex);
}

void output_close(void)
{
    if (out) fclose(out);
    out = NULL;
}

void print_interrog_line(FILE *f, const struct interrog *it)
{
    char desc[160] = "";
    if (kind_is_modes(it->kind) && it->kind != K_S_NOSYNC) modes_describe(it, desc, sizeof(desc));
    char p2[16] = "";
    if (!isnan(it->p2_db)) snprintf(p2, sizeof(p2), "P2%+.0fdB", it->p2_db);
    char q[12] = "";
    if (kind_is_modes(it->kind) && it->kind != K_S_NOSYNC) snprintf(q, sizeof(q), "q%.2f ", it->quality);
    char rc[8];
    radar_code_str(it->radar, rc, sizeof(rc));
    fprintf(f, "%12.6f %-24s %6.1fdBFS %-8s %-4s %s%s%s%s%s %s%s\n", it->t, kind_label(it->kind),
            it->amp_db, p2, rc, it->s1 ? "S1 " : "", it->sidelobe ? "SLS " : "",
            it->tcas ? "TCAS " : "", it->len_mismatch ? "LEN? " : "",
            it->burst >= 0 ? "B" : "-", q, desc);
}
