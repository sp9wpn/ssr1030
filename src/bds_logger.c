/* bds_logger: reads a Beast stream (dump1090, readsb, ...) and prints the
 * Comm-B replies (DF20/21) of ADS-B aircraft that are BDS 1,0 (data link
 * capability), 2,0 (identification), 3,0 (ACAS RA), 4,0 (selected altitude),
 * 4,4 / 4,5 (meteo routine / hazard), 5,0 (track and turn), 5,1 / 5,2
 * (position report coarse / fine) or 6,0 (heading and speed), one line each
 * (--bds chooses which; default all but 5,0).
 *
 * A Comm-B reply does not say which register it holds; that is in the
 * interrogation, which a 1090 receiver does not hear. So the register is
 * inferred from the content, as in pyModeS:
 *  - every field with its status bit clear must be all zeros, and the
 *    values must be in range (roll, speeds, rates, Mach);
 *  - BDS 1,0, 2,0 and 3,0 carry their own code in MB bits 1-8;
 *  - BDS 4,4 and 4,5 only from data parity replies (below);
 *  - the aircraft must be known from its own ADS-B (DF17) in the last minute,
 *    and the content must agree with its ADS-B velocity: BDS 5,0 ground
 *    speed and track, BDS 6,0 heading (wind correction) and Mach / IAS
 *    against its altitude and ground speed; BDS 4,0 is checked on its own
 *    (reserved bits, selected altitude, barometric setting), BDS 5,1 and 5,2
 *    against its pressure altitude and ADS-B position (a 5,1 of an aircraft
 *    without one: against the area of the ADS-B positions heard);
 *  - a reply that fits more than one register is printed with each decode
 *    marked '?'.
 * The address of a DF20/21 reply is its parity syndrome (AP field), or, when
 * the interrogator asked for data parity (DP), the syndrome XOR (BDS1,BDS2
 * << 16). A DP reply thus names its own register: it is decoded as that
 * register only, with the range checks but without the ADS-B comparison
 * (which fails on the runway, where there is no airborne ADS-B velocity),
 * and marked 'D'. A reply with bit errors in the top of its parity field
 * looks like data parity for some register too (bit 92 flipped: 1,0), so
 * the content checks apply also there.
 *
 * Callsign: from ADS-B identification and BDS 2,0 replies, the latest.
 *
 * Altitude: the last pressure altitude (DF20 reply, DF0/4/16, ADS-B), made
 * AMSL with --qnh, or else with a QNH read from BDS 4,0 replies (barometric
 * setting) of aircraft below 4000 ft that are not set to standard.
 */
#define _GNU_SOURCE
#include <errno.h>
#include <math.h>
#include <netdb.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#define AC_N        4096            /* aircraft table */
#define KNOWN_AGE   60.0            /* s: an aircraft is known this long after its last DF17 */
#define VEL_AGE     30.0            /* s: an ADS-B velocity is used this long */
#define ALT_AGE     60.0            /* s: a pressure altitude is shown this long */
#define QNH_N       16              /* barometric settings kept for the derived QNH */
#define QNH_AGE     1800.0          /* s */
#define QNH_MAX_FT  4000            /* settings only from aircraft below this (above the transition: standard) */
#define POS_AGE     10.0            /* s: an ADS-B position checks BDS 5,1 / 5,2 this long (250 m/s: 0.02 deg) */
#define POS_TOL     0.05            /* deg: BDS 5,1 / 5,2 position within this of the ADS-B one */
#define POS_ALT_FT  100             /* BDS 5,1/5,2 altitude within this of the known pressure altitude */
#define B51_AGE     300.0           /* s: a BDS 5,1 position resolves 5,2 (+-0.7 deg) this long */
#define CPR_AGE     10.0            /* s: even and odd CPR halves this close decode globally */
#define AREA_DEG    4.0             /* deg: a BDS 5,1 of an aircraft without an ADS-B position must be
                                       this close to the centre of the ADS-B positions heard */
#define AREA_MIN    20              /* positions decoded before the centre is used */
#define AREA_AGE    600.0           /* s: the centre is used this long after the last position */

struct ac {
    uint32_t addr;
    bool used;
    double seen;                    /* last DF17 */
    char cs[9];
    double vel_t;                   /* ADS-B velocity: ground speed kt, track deg, vertical rate ft/min */
    double gs, trk, vr;
    bool vr_ok;
    double alt_t;
    int alt;                        /* pressure altitude, ft */
    double cpr_t[2];                /* ADS-B CPR position halves: even, odd */
    unsigned cpr_lat[2], cpr_lon[2];
    double pos_t, lat, lon;         /* ADS-B position (global CPR decode) */
    double b51_t, b51_lat, b51_lon; /* last BDS 5,1 position: resolves 5,2 without ADS-B */
};

static struct ac act[AC_N];
static double area_lat, area_lon, area_t = -1e300;  /* running centre of the ADS-B positions */
static long area_n;
static double qnh_cmd = NAN;        /* --qnh */
static int max_alt = INT32_MAX;     /* --max-alt: AMSL with a QNH, else pressure altitude, ft */
#define SEL_N 32
static char sel[SEL_N][9];          /* -A: hex addresses or callsigns */
static int nsel;
static double qset[QNH_N], qset_t[QNH_N];
static int qset_n;
static bool verbose;
/* registers decoded */
enum { R10, R20, R30, R40, R44, R45, R50, R51, R52, R60, NREG };
static const uint8_t reg_code[NREG] = { 0x10, 0x20, 0x30, 0x40, 0x44, 0x45, 0x50, 0x51, 0x52, 0x60 };
static bool show[NREG] = { true, true, true, true, true, true, false, true, true, true };   /* --bds */
static long n_reg[NREG];            /* replies that fit only this register */
static long n_commb, n_unknown, n_both, n_none, n_high, n_other, n_dp, n_dp_other;

/* registers that may be requested with data parity (the syndrome is then
 * address ^ (BDS << 16)) */
static const uint8_t dp_registers[] = {
    0x10, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x30,
    0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x48, 0x50, 0x51, 0x52, 0x53, 0x5f, 0x60
};

static double now_s(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static unsigned getbits(const uint8_t *m, int first, int last)
{
    unsigned v = 0;
    for (int b = first; b <= last; b++) v = v << 1 | ((m[(b - 1) >> 3] >> (7 - ((b - 1) & 7))) & 1);
    return v;
}

/* Mode S CRC-24 remainder over the whole message (0 for a good DF17, the address for DF20/21) */
static uint32_t crc(const uint8_t *m, int nbits)
{
    uint32_t r = 0;
    for (int b = 0; b < nbits - 24; b++) {
        uint32_t bit = (m[b >> 3] >> (7 - (b & 7))) & 1;
        uint32_t top = (r >> 23) & 1;
        r = (r << 1 | bit) & 0xFFFFFF;
        if (top) r ^= 0xFFF409;
    }
    for (int b = 0; b < 24; b++) {
        uint32_t top = (r >> 23) & 1;
        r = (r << 1) & 0xFFFFFF;
        if (top) r ^= 0xFFF409;
    }
    return r ^ getbits(m, nbits - 23, nbits);
}

static struct ac *ac_get(uint32_t a, bool create)
{
    uint32_t h = (a * 2654435761u) >> 20;
    for (int k = 0; k < AC_N; k++) {
        struct ac *e = &act[(h + k) % AC_N];
        if (!e->used) {
            if (!create) return NULL;
            memset(e, 0, sizeof(*e));
            e->used = true;
            e->addr = a;
            e->seen = e->vel_t = e->alt_t = e->pos_t = e->b51_t = e->cpr_t[0] = e->cpr_t[1] = -1e300;
            return e;
        }
        if (e->addr == a) return e;
    }
    return NULL;
}

/* ---------------- altitude ---------------- */

/* 13-bit field (C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4) -> 0xABCD, one octal digit per nibble */
static unsigned abcd13(unsigned v)
{
#define B13(n) ((v >> (13 - (n))) & 1)
    return B13(6) << 14 | B13(4) << 13 | B13(2) << 12 | B13(12) << 10 | B13(10) << 9 | B13(8) << 8 |
           B13(5) << 6 | B13(3) << 5 | B13(1) << 4 | B13(13) << 2 | B13(11) << 1 | B13(9);
#undef B13
}

/* Gillham code -> 100 ft units, INT32_MIN if invalid */
static int gillham(unsigned a)
{
    if ((a & 0xFFFF8889) || !(a & 0xF0)) return INT32_MIN;
    unsigned h = 0, f = 0;
    if (a & 0x10) h ^= 7;
    if (a & 0x20) h ^= 3;
    if (a & 0x40) h ^= 1;
    if ((h & 5) == 5) h ^= 2;
    if (h > 5) return INT32_MIN;
    if (a & 0x0002) f ^= 0xFF;
    if (a & 0x0004) f ^= 0x7F;
    if (a & 0x1000) f ^= 0x3F;
    if (a & 0x2000) f ^= 0x1F;
    if (a & 0x4000) f ^= 0x0F;
    if (a & 0x0100) f ^= 0x07;
    if (a & 0x0200) f ^= 0x03;
    if (a & 0x0400) f ^= 0x01;
    if (f & 1) h = 6 - h;
    return (int)(f * 5 + h) - 13;
}

/* 13-bit altitude code (DF0/4/16/20) -> feet, INT32_MIN if absent, metric or invalid */
static int ac13_ft(unsigned v)
{
    if (!v || (v & 0x40)) return INT32_MIN;
    if (v & 0x10) return (int)(((v & 0x1F80) >> 2 | (v & 0x20) >> 1 | (v & 0x0F)) * 25) - 1000;
    int h = gillham(abcd13(v));
    return h == INT32_MIN ? INT32_MIN : h * 100;
}

/* ADS-B 12-bit altitude: the 13-bit code without its M bit */
static int ac12_ft(unsigned v) { return ac13_ft((v & 0xFC0) << 1 | (v & 0x3F)); }

static void set_alt(struct ac *a, int ft, double t)
{
    if (ft == INT32_MIN) return;
    a->alt = ft;
    a->alt_t = t;
}

/* QNH: --qnh, else the median of recent BDS 4,0 settings; NAN when none */
static double qnh_now(double t)
{
    if (!isnan(qnh_cmd)) return qnh_cmd;
    double v[QNH_N];
    int n = 0;
    for (int k = 0; k < qset_n && k < QNH_N; k++)
        if (t - qset_t[k] < QNH_AGE) v[n++] = qset[k];
    if (!n) return NAN;
    for (int i = 1; i < n; i++)
        for (int j = i; j > 0 && v[j] < v[j - 1]; j--) { double x = v[j]; v[j] = v[j - 1]; v[j - 1] = x; }
    return n & 1 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

/* pressure altitude -> AMSL with the ISA pressure-height relation */
static double amsl_ft(int pa, double qnh)
{
    double h_qnh = 145366.45 * (1 - pow(qnh / 1013.25, 0.190284));   /* pressure altitude of the QNH level */
    return pa - h_qnh;
}

/* ---------------- ADS-B (DF17) ---------------- */

static const char ais[] = "#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######";

/* CPR number of longitude zones at a latitude */
static int cpr_nl(double lat)
{
    lat = fabs(lat);
    if (lat < 1e-9) return 59;
    if (lat >= 87) return lat > 87 ? 1 : 2;
    double c = cos(M_PI / 180 * lat);
    return (int)floor(2 * M_PI / acos(1 - (1 - cos(M_PI / 30)) / (c * c)));
}

static double cpr_mod(double a, double b) { double r = fmod(a, b); return r < 0 ? r + b : r; }

/* airborne position: global decode from an even and an odd half, the newer one f */
static void cpr_decode(struct ac *a, int f, double t)
{
    if (fabs(a->cpr_t[0] - a->cpr_t[1]) > CPR_AGE) return;
    double lat0 = a->cpr_lat[0] / 131072.0, lat1 = a->cpr_lat[1] / 131072.0;
    double lon0 = a->cpr_lon[0] / 131072.0, lon1 = a->cpr_lon[1] / 131072.0;
    double j = floor(59 * lat0 - 60 * lat1 + 0.5);
    double rl0 = 360.0 / 60 * (cpr_mod(j, 60) + lat0), rl1 = 360.0 / 59 * (cpr_mod(j, 59) + lat1);
    if (rl0 >= 270) rl0 -= 360;
    if (rl1 >= 270) rl1 -= 360;
    if (cpr_nl(rl0) != cpr_nl(rl1)) return;                    /* the halves straddle a zone boundary */
    double lat = f ? rl1 : rl0;
    int ni = cpr_nl(lat) - f;
    if (ni < 1) ni = 1;
    double m = floor(lon0 * (cpr_nl(lat) - 1) - lon1 * cpr_nl(lat) + 0.5);
    double lon = 360.0 / ni * (cpr_mod(m, ni) + (f ? lon1 : lon0));
    if (lon >= 180) lon -= 360;
    a->lat = lat; a->lon = lon; a->pos_t = t;
    /* reception area: plain mean of the first positions, then exponential */
    double k = area_n < 50 ? 1.0 / (area_n + 1) : 0.02;
    area_lat += (lat - area_lat) * k;
    area_lon = remainder(area_lon + remainder(lon - area_lon, 360) * k, 360);
    area_n++;
    area_t = t;
}

static void df17(const uint8_t *m, double t)
{
    struct ac *a = ac_get(getbits(m, 9, 32), true);
    if (!a) return;
    a->seen = t;
    int tc = getbits(m, 33, 37);
    if (tc >= 1 && tc <= 4) {
        for (int k = 0; k < 8; k++) a->cs[k] = ais[getbits(m, 41 + 6 * k, 46 + 6 * k)];
        a->cs[8] = 0;
        for (int k = 7; k >= 0 && (a->cs[k] == ' ' || a->cs[k] == '#'); k--) a->cs[k] = 0;
    } else if ((tc >= 9 && tc <= 18) || (tc >= 20 && tc <= 22)) {
        if (tc <= 18) set_alt(a, ac12_ft(getbits(m, 41, 52)), t);
        int f = getbits(m, 54, 54);
        a->cpr_lat[f] = getbits(m, 55, 71); a->cpr_lon[f] = getbits(m, 72, 88); a->cpr_t[f] = t;
        cpr_decode(a, f, t);
    } else if (tc == 19) {
        int st = getbits(m, 38, 40);
        if (st == 1 || st == 2) {
            int ew = getbits(m, 47, 56), ns = getbits(m, 58, 67);
            if (!ew || !ns) return;
            double k = st == 2 ? 4 : 1;
            double vx = (getbits(m, 46, 46) ? -1 : 1) * (ew - 1) * k;
            double vy = (getbits(m, 57, 57) ? -1 : 1) * (ns - 1) * k;
            a->gs = hypot(vx, vy);
            a->trk = fmod(atan2(vx, vy) * 180 / M_PI + 360, 360);
            int vr = getbits(m, 70, 78);
            a->vr_ok = vr != 0;
            a->vr = (getbits(m, 69, 69) ? -1 : 1) * (vr - 1) * 64.0;
            a->vel_t = t;
        }
    }
}

/* ---------------- Comm-B inference ---------------- */

/* signed field with its sign bit just before it */
static double sgn(const uint8_t *m, int sbit, int first, int last)
{
    int v = (int)getbits(m, first, last);
    return getbits(m, sbit, sbit) ? v - (1 << (last - first + 1)) : v;
}

#include "commb_fmt.h"

/* a field whose status bit is clear must be all zeros */
static bool field_ok(const uint8_t *m, int status, int last)
{
    return getbits(m, status, status) || !getbits(m, status + 1, last);
}

static double angdiff(double a, double b)
{
    double d = fmod(fabs(a - b), 360);
    return d > 180 ? 360 - d : d;
}

struct b50 { bool roll_ok, trk_ok, gs_ok, rate_ok, tas_ok; double roll, trk, gs, rate, tas; };
struct b60 { bool hdg_ok, ias_ok, mach_ok, baro_ok, ins_ok; double hdg, ias, mach, baro, ins; };

/* -A lists this address as hex: its replies are decoded even without ADS-B */
static bool hex_given(uint32_t addr)
{
    for (int k = 0; k < nsel; k++) {
        char *e;
        unsigned long v = strtoul(sel[k], &e, 16);
        if (strlen(sel[k]) == 6 && !*e && v == addr) return true;
    }
    return false;
}

/* ref: an ADS-B velocity to check against. Without one (an aircraft given as
 * hex with -A and no recent ADS-B), every field must be present instead */
static bool is50(const uint8_t *m, const struct ac *a, bool ref, bool dp, struct b50 *o)
{
    if (!field_ok(m, 33, 43) || !field_ok(m, 44, 55) || !field_ok(m, 56, 66) || !field_ok(m, 67, 77) || !field_ok(m, 78, 88))
        return false;
    o->roll_ok = getbits(m, 33, 33); o->roll = sgn(m, 34, 35, 43) * 45.0 / 256;
    o->trk_ok = getbits(m, 44, 44);  o->trk = fmod(sgn(m, 45, 46, 55) * 90.0 / 512 + 360, 360);
    o->gs_ok = getbits(m, 56, 56);   o->gs = getbits(m, 57, 66) * 2.0;
    o->rate_ok = getbits(m, 67, 67); o->rate = sgn(m, 68, 69, 77) * 8.0 / 256;
    o->tas_ok = getbits(m, 78, 78);  o->tas = getbits(m, 79, 88) * 2.0;
    if (!o->trk_ok || !o->gs_ok) return false;                  /* the fields checked against ADS-B */
    if (o->roll_ok && fabs(o->roll) > 50) return false;
    if (o->gs < 40 || o->gs > 700) return false;
    if (o->tas_ok && (o->tas < 40 || o->tas > 600 || fabs(o->gs - o->tas) > 200)) return false;
    if (o->rate_ok && fabs(o->rate) > 8) return false;
    if (dp) return true;                                        /* register named by the reply itself */
    if (!ref) return o->roll_ok && o->tas_ok;
    return fabs(o->gs - a->gs) <= 25 && angdiff(o->trk, a->trk) <= 12;
}

/* Mach from calibrated airspeed (kt) at pressure altitude (ft), ISA */
static double cas_to_mach(double cas, double alt_ft)
{
    double p = 101325 * pow(1 - 6.8756e-6 * alt_ft, 5.2559);
    if (alt_ft > 36089) p = 22632 * exp(-4.8063e-5 * (alt_ft - 36089));
    double qc = 101325 * (pow(1 + 0.2 * pow(cas / 661.47, 2), 3.5) - 1);
    return sqrt(5 * (pow(qc / p + 1, 2.0 / 7) - 1));
}

/* ISA speed of sound at pressure altitude, kt */
static double sound_kt(double alt_ft)
{
    double tk = alt_ft < 36089 ? 288.15 - 1.98e-3 * alt_ft : 216.65;
    return 38.967854 * sqrt(tk);
}

static bool is60(const uint8_t *m, const struct ac *a, double t, bool ref, bool dp, struct b60 *o)
{
    if (!field_ok(m, 33, 44) || !field_ok(m, 45, 55) || !field_ok(m, 56, 66) || !field_ok(m, 67, 77) || !field_ok(m, 78, 88))
        return false;
    o->hdg_ok = getbits(m, 33, 33);  o->hdg = fmod(sgn(m, 34, 35, 44) * 90.0 / 512 + 360, 360);
    o->ias_ok = getbits(m, 45, 45);  o->ias = getbits(m, 46, 55);
    o->mach_ok = getbits(m, 56, 56); o->mach = getbits(m, 57, 66) * 2.048 / 512;
    o->baro_ok = getbits(m, 67, 67); o->baro = sgn(m, 68, 69, 77) * 32;
    o->ins_ok = getbits(m, 78, 78);  o->ins = sgn(m, 79, 80, 88) * 32;
    if (!o->hdg_ok || (!o->ias_ok && !o->mach_ok)) return false;
    if (o->ias_ok && (o->ias < 40 || o->ias > 500)) return false;
    if (o->mach_ok && (o->mach < 0.1 || o->mach > 0.95)) return false;
    if (o->baro_ok && fabs(o->baro) > 6000) return false;
    if (o->ins_ok && fabs(o->ins) > 6000) return false;
    if (o->baro_ok && o->ins_ok && fabs(o->baro - o->ins) > 2000) return false;
    bool alt = t - a->alt_t <= ALT_AGE;
    if (dp)                                                     /* register named by the reply itself */
        return !(o->ias_ok && o->mach_ok && alt && fabs(cas_to_mach(o->ias, a->alt) - o->mach) > 0.04);
    if (!ref)
        return o->ias_ok && o->mach_ok && (o->baro_ok || o->ins_ok) &&
               (!alt || fabs(cas_to_mach(o->ias, a->alt) - o->mach) <= 0.04);
    if (angdiff(o->hdg, a->trk) > 40) return false;             /* wind correction angle */
    if (a->vr_ok && o->ins_ok && fabs(o->ins - a->vr) > 1500) return false;
    if (a->vr_ok && !o->ins_ok && o->baro_ok && fabs(o->baro - a->vr) > 2000) return false;
    if (o->ias_ok && o->mach_ok && alt && fabs(cas_to_mach(o->ias, a->alt) - o->mach) > 0.04) return false;
    double m_ = o->mach_ok ? o->mach : alt ? cas_to_mach(o->ias, a->alt) : NAN;
    if (alt && !isnan(m_) && fabs(m_ * sound_kt(a->alt) - a->gs) > 200) return false;   /* TAS vs ground speed: wind */
    return true;
}

struct b40 { bool mcp_ok, fms_ok, baro_ok, mode_ok, src_ok; double mcp, fms, baro; bool vnav, alth, app; int src; };

/* BDS 4,0 (selected vertical intention): reserved bits 72-79 and 84-85 zero,
 * a selected altitude near a multiple of 100 ft (16 ft resolution), a
 * plausible barometric setting. It needs no ADS-B to check against */
static bool is40(const uint8_t *m, struct b40 *o)
{
    if (getbits(m, 72, 79) || getbits(m, 84, 85)) return false;
    if (!field_ok(m, 33, 45) || !field_ok(m, 46, 58) || !field_ok(m, 59, 71) || !field_ok(m, 80, 83) || !field_ok(m, 86, 88))
        return false;
    o->mcp_ok = getbits(m, 33, 33);  o->mcp = getbits(m, 34, 45) * 16.0;
    o->fms_ok = getbits(m, 46, 46);  o->fms = getbits(m, 47, 58) * 16.0;
    o->baro_ok = getbits(m, 59, 59); o->baro = getbits(m, 60, 71) * 0.1 + 800;
    o->mode_ok = getbits(m, 80, 80); o->vnav = getbits(m, 81, 81); o->alth = getbits(m, 82, 82); o->app = getbits(m, 83, 83);
    o->src_ok = getbits(m, 86, 86);  o->src = getbits(m, 87, 88);
    if (!o->mcp_ok && !o->fms_ok) return false;
    if (o->mcp_ok && (o->mcp > 50000 || fabs(o->mcp - 100 * round(o->mcp / 100)) > 16)) return false;
    if (o->fms_ok && (o->fms > 50000 || fabs(o->fms - 100 * round(o->fms / 100)) > 16)) return false;
    if (o->baro_ok && (o->baro < 950 || o->baro > 1060)) return false;
    return true;
}

struct b51 { double lat, lon, alt; };
struct b52 { int fom; double fla, flo, alt, lat, lon; bool gnss, resolved; };
#define B51_LSB  (360.0 / 1048576)
#define B52_LSB  (90.0 / 16777216)
#define B52_SPAN (180.0 / 128)          /* the fine position is the position modulo this */

/* position error, deg of latitude (longitude scaled by cos lat) */
static double pos_err(double la, double lo, double rla, double rlo)
{
    return fmax(fabs(la - rla), angdiff(lo, rlo) * cos(rla * M_PI / 180));
}

/* BDS 5,1 (position report coarse): status, latitude and longitude (sign + 19
 * bits, LSB 360/2^20 deg), pressure altitude (sign + 14 bits, LSB 8 ft).
 * Inferred (no dp): the altitude must match the aircraft's pressure altitude,
 * and the position its ADS-B one, or without one (*area set) lie in the
 * reception area. err: the error to the ADS-B position, deg (NAN unchecked) */
static bool is51(const uint8_t *m, const struct ac *a, double t, bool dp, struct b51 *o, double *err, bool *area)
{
    *err = NAN;
    *area = false;
    if (!getbits(m, 33, 33)) return false;
    o->lat = sgn(m, 34, 35, 53) * B51_LSB;
    o->lon = sgn(m, 54, 55, 73) * B51_LSB;
    o->alt = sgn(m, 74, 75, 88) * 8;
    if (fabs(o->lat) > 90 || o->alt < -1500 || o->alt > 60000) return false;
    if (t - a->pos_t <= POS_AGE) *err = pos_err(o->lat, o->lon, a->lat, a->lon);
    if (dp) return true;
    if (t - a->alt_t > ALT_AGE || fabs(o->alt - a->alt) > POS_ALT_FT) return false;
    if (!isnan(*err)) return *err <= POS_TOL;
    if (area_n < AREA_MIN || t - area_t > AREA_AGE) return false;
    *area = true;
    return pos_err(o->lat, o->lon, area_lat, area_lon) <= AREA_DEG;
}

/* BDS 5,2 (position report fine): status, FOM/source (4 bits), latitude and
 * longitude fine (18-bit two's complement, LSB 90/2^24 deg: the position
 * modulo 180/128 deg), altitude as in 5,1 (GNSS height with FOM 11-13).
 * Resolved with the ADS-B position, else the last 5,1. Inferred (no dp): FOM
 * 1-5, the fine position must match the ADS-B one and the altitude its
 * pressure altitude. 5,1 and 5,2 share the status and altitude fields, so the
 * position is what tells them apart */
static bool is52(const uint8_t *m, const struct ac *a, double t, bool dp, struct b52 *o, double *err)
{
    *err = NAN;
    if (!getbits(m, 33, 33)) return false;
    o->fom = (int)getbits(m, 34, 37);
    o->fla = sgn(m, 38, 39, 55) * B52_LSB;
    o->flo = sgn(m, 56, 57, 73) * B52_LSB;
    o->alt = sgn(m, 74, 75, 88) * 8;
    o->gnss = o->fom >= 11 && o->fom <= 13;
    double rla = NAN, rlo = NAN;
    if (t - a->pos_t <= POS_AGE) { rla = a->lat; rlo = a->lon; }
    else if (t - a->b51_t <= B51_AGE) { rla = a->b51_lat; rlo = a->b51_lon; }
    o->resolved = !isnan(rla);
    if (o->resolved) {
        o->lat = rla + remainder(o->fla - rla, B52_SPAN);
        o->lon = remainder(rlo + remainder(o->flo - rlo, B52_SPAN), 360);
        if (t - a->pos_t <= POS_AGE) *err = pos_err(o->lat, o->lon, rla, rlo);
    }
    if (o->alt < -1500 || o->alt > 60000) return false;
    if (dp) return o->fom <= 5 || o->gnss;
    if (o->fom < 1 || o->fom > 5 || isnan(*err) || *err > POS_TOL) return false;
    return t - a->alt_t <= ALT_AGE && fabs(o->alt - a->alt) <= POS_ALT_FT;
}

/* BDS 1,0, 2,0 and 3,0 carry their own code in MB bits 1-8, so they need no
 * ADS-B check; the other fields must still be consistent, since a reply with
 * a bit error in its parity field can look like data parity for them */
static bool is10(const uint8_t *m)
{
    return getbits(m, 33, 40) == 0x10 && !getbits(m, 42, 46) && getbits(m, 49, 55) <= 7;   /* reserved, subnet version */
}

static bool is20(const uint8_t *m, char cs[9])
{
    if (getbits(m, 33, 40) != 0x20) return false;
    for (int k = 0; k < 8; k++)
        if ((cs[k] = ais[getbits(m, 41 + 6 * k, 46 + 6 * k)]) == '#') return false;
    cs[8] = 0;
    for (int k = 7; k >= 0 && cs[k] == ' '; k--) cs[k] = 0;
    return cs[0] && cs[0] != ' ';
}

static bool is30(const uint8_t *m)
{
    if (getbits(m, 33, 40) != 0x30 || getbits(m, 48, 54) >= 48) return false;   /* ARA bits reserved for ACAS III */
    unsigned tti = getbits(m, 61, 62);
    return tti != 3 && (tti || !getbits(m, 63, 88)) && (tti != 1 || !getbits(m, 87, 88));
}

/* BDS 4,4 and 4,5 (meteorological) say nothing about themselves and are
 * rare, so they are taken only from data parity replies, with the checks
 * pyModeS uses */
static bool is44(const uint8_t *m)
{
    if (getbits(m, 33, 36) > 4 || !getbits(m, 37, 37) || getbits(m, 38, 46) > 250) return false;   /* FOM, wind */
    if (!field_ok(m, 67, 78) || !field_ok(m, 79, 81) || !field_ok(m, 82, 88)) return false;
    double sat = sat44(m);
    return sat >= -80 && sat <= 60 && (getbits(m, 38, 46) || getbits(m, 47, 55) || getbits(m, 57, 66));
}

static bool is45(const uint8_t *m)
{
    if (getbits(m, 84, 88)) return false;                       /* reserved */
    if (getbits(m, 39, 39) && !getbits(m, 57, 88)) return false;  /* the shape of BDS 1,7 */
    for (int k = 0; k < 5; k++) if (!field_ok(m, 33 + 3 * k, 35 + 3 * k)) return false;
    if (!field_ok(m, 48, 58) || !field_ok(m, 59, 70) || !field_ok(m, 71, 83)) return false;
    double sat = sgn(m, 49, 50, 58) * 0.25;
    if (getbits(m, 48, 48) && (sat < -80 || sat > 60)) return false;
    return getbits(m, 33, 88) != 0;
}

static const char *fom_str(int f)
{
    static const char *s[] = { "no nav", "RNP 20", "RNP 5", "RNP 1", "DGNSS 10 m", "DGNSS 2.5 m" };
    return f <= 5 ? s[f] : f <= 13 && f >= 11 ? s[f - 8] : "FOM ?";
}

/* QNH from a BDS 4,0 barometric setting: aircraft low enough to fly on QNH,
 * setting not standard */
static void qnh_note(const struct b40 *f, const struct ac *a, double t)
{
    if (!isnan(qnh_cmd) || !f->baro_ok || t - a->alt_t > ALT_AGE || a->alt > QNH_MAX_FT) return;
    if (fabs(f->baro - 1013.2) < 0.4) return;
    double before = qnh_now(t);
    qset[qset_n % QNH_N] = f->baro;
    qset_t[qset_n % QNH_N] = t;
    qset_n++;
    double after = qnh_now(t);
    if (isnan(before) || fabs(after - before) >= 0.5)
        fprintf(stderr, "# QNH %.1f hPa (from BDS 4,0 of %06x %s, %d ft)\n", after, a->addr, a->cs, a->alt);
}

static void print_head(const struct ac *a, double t, double rssi, char *buf, size_t n)
{
    time_t s = (time_t)t;
    struct tm tm;
    localtime_r(&s, &tm);
    int o = snprintf(buf, n, "%02d:%02d:%02d.%03d  %06x  %-8s  ", tm.tm_hour, tm.tm_min, tm.tm_sec,
                     (int)((t - floor(t)) * 1000), a->addr, a->cs[0] ? a->cs : "-");
    double q = qnh_now(t);
    if (t - a->alt_t > ALT_AGE) o += snprintf(buf + o, n - o, "%14s", "alt ?");
    else if (isnan(q)) o += snprintf(buf + o, n - o, "%8d ft STD", a->alt);
    else o += snprintf(buf + o, n - o, "%9.0f ft   ", amsl_ft(a->alt, q));
    if (isnan(rssi)) snprintf(buf + o, n - o, "  %9s", "-");
    else snprintf(buf + o, n - o, "  %5.1f dBFS", rssi);
}

/* -A: the aircraft's hex address or callsign is in the list */
static bool selected(const struct ac *a)
{
    if (!nsel) return true;
    char hex[8];
    snprintf(hex, sizeof(hex), "%06x", a->addr);
    for (int k = 0; k < nsel; k++)
        if (!strcasecmp(sel[k], hex) || (a->cs[0] && !strcasecmp(sel[k], a->cs))) return true;
    return false;
}

static void commb(const uint8_t *m, int df, double t, double rssi)
{
    n_commb++;
    uint32_t syn = crc(m, 112), addr = syn;
    struct ac *a = ac_get(addr, false);
    bool given = hex_given(addr);
    unsigned dp = 0;                                            /* register named by data parity */
    if ((!a || t - a->seen > KNOWN_AGE) && !given) {
        /* data parity: syndrome = address ^ (BDS << 16) */
        for (unsigned k = 0; k < sizeof(dp_registers) && !dp; k++) {
            uint32_t x = syn ^ ((uint32_t)dp_registers[k] << 16);
            struct ac *b = ac_get(x, false);
            if ((b && t - b->seen <= KNOWN_AGE) || hex_given(x)) { dp = dp_registers[k]; addr = x; a = b; given = hex_given(x); }
        }
        if (!dp) { n_unknown++; return; }
        n_dp++;
    }
    if (!a && !(a = ac_get(addr, true))) return;
    if (df == 20) set_alt(a, ac13_ft(getbits(m, 20, 32)), t);
    int dreg = -1;
    for (int k = 0; k < NREG && dp; k++) if (reg_code[k] == dp) dreg = k;
    if (dp && dreg < 0) { n_dp_other++; return; }
    struct b40 f4;
    struct b50 f5;
    struct b51 f51;
    struct b52 f52;
    struct b60 f6;
    char cs[9];
    bool r[NREG] = { false }, ref = t - a->vel_t <= VEL_AGE, area51 = false;
    double e51 = NAN, e52 = NAN;
    if (dp) {
        /* only the named register, no ADS-B comparison */
        switch (dp) {
        case 0x10: r[R10] = is10(m); break;
        case 0x20: r[R20] = is20(m, cs); break;
        case 0x30: r[R30] = is30(m); break;
        case 0x40: r[R40] = is40(m, &f4); break;
        case 0x44: r[R44] = is44(m); break;
        case 0x45: r[R45] = is45(m); break;
        case 0x50: r[R50] = is50(m, a, ref, true, &f5); break;
        case 0x51: r[R51] = is51(m, a, t, true, &f51, &e51, &area51); break;
        case 0x52: r[R52] = is52(m, a, t, true, &f52, &e52); break;
        case 0x60: r[R60] = is60(m, a, t, ref, true, &f6); break;
        }
    } else {
        r[R10] = is10(m);
        r[R20] = is20(m, cs);
        r[R30] = is30(m);
        r[R40] = is40(m, &f4);
        /* 5,0 and 6,0 need an ADS-B velocity to check against, unless given as hex */
        bool chk = ref || given;
        r[R50] = chk && is50(m, a, ref, false, &f5);
        r[R51] = is51(m, a, t, false, &f51, &e51, &area51);
        r[R52] = is52(m, a, t, false, &f52, &e52);
        r[R60] = chk && is60(m, a, t, ref, false, &f6);
        /* both: a genuine one is close to the ADS-B position, the other reading lands anywhere */
        if (r[R51] && r[R52]) { if (e51 > 3 * e52) r[R51] = false; else if (e52 > 3 * e51) r[R52] = false; }
    }
    int nfit = 0, only = -1;
    for (int k = 0; k < NREG; k++) if (r[k]) { nfit++; only = k; }
    if (r[R40]) qnh_note(&f4, a, t);
    if (r[R51] && nfit == 1) { a->b51_lat = f51.lat; a->b51_lon = f51.lon; a->b51_t = t; }   /* resolves a later 5,2 */
    /* the callsign: from BDS 2,0 as well as ADS-B identification, the latest wins */
    char cs_was[9] = "";
    if (r[R20] && nfit == 1 && strcmp(a->cs, cs)) { memcpy(cs_was, a->cs, sizeof(cs_was)); memcpy(a->cs, cs, sizeof(a->cs)); }
    if (!nfit) { n_none++; return; }
    /* altitude filter: above it, or not known, not shown */
    if (max_alt != INT32_MAX) {
        double q = qnh_now(t), shown = isnan(q) ? a->alt : amsl_ft(a->alt, q);   /* the altitude the line shows */
        if (t - a->alt_t > ALT_AGE || shown > max_alt) { n_high++; return; }
    }
    if (!selected(a)) { n_other++; return; }
    /* fits more than one register: each decode marked '?' (also when the
     * other one is not shown) */
    const char *mark = dp ? "D" : nfit > 1 ? "?" : " ";
    const char *nochk = ref || dp ? "" : "  (no ADS-B check)";
    if (nfit > 1) n_both++;
    else n_reg[only]++;
    char head[96], body[256];
    print_head(a, t, rssi, head, sizeof(head));
    for (int k = 0; k < NREG; k++) {
        if (!r[k] || !show[k]) continue;
        int o = snprintf(body, sizeof(body), "BDS%X,%X%s", reg_code[k] >> 4, reg_code[k] & 15, mark);
        const char *tail = "";
        size_t n = sizeof(body) - o;
        char *b = body + o;
        switch (reg_code[k]) {
        case 0x10: fmt10(m, b, n); break;
        case 0x20:
            o = snprintf(b, n, " id %s", cs);
            if (cs_was[0]) snprintf(b + o, n - o, "  (was %s)", cs_was);
            else if (strcmp(a->cs, cs)) snprintf(b + o, n - o, "  (also fits another register; callsign %s kept)", a->cs[0] ? a->cs : "-");
            break;
        case 0x30: fmt30(m, b, n); break;
        case 0x40: {
            o = 0;
            if (f4.mcp_ok) o += snprintf(b + o, n - o, " MCP %5.0f ft", f4.mcp);
            if (f4.fms_ok) o += snprintf(b + o, n - o, "  FMS %5.0f ft", f4.fms);
            if (f4.baro_ok) o += snprintf(b + o, n - o, "  baro %6.1f hPa", f4.baro);
            if (f4.mode_ok) o += snprintf(b + o, n - o, "  mode%s%s%s%s", f4.vnav ? " VNAV" : "", f4.alth ? " ALT" : "",
                                          f4.app ? " APP" : "", f4.vnav || f4.alth || f4.app ? "" : " -");
            static const char *srcs[] = { "unknown", "aircraft alt", "MCP/FCU", "FMS" };
            if (f4.src_ok) snprintf(b + o, n - o, "  target %s", srcs[f4.src]);
            break;
        }
        case 0x44: fmt44(m, b, n); break;
        case 0x45: fmt45(m, b, n); break;
        case 0x50:
            o = 0;
            if (f5.roll_ok) o += snprintf(b + o, n - o, " roll %+5.1f°", f5.roll);
            o += snprintf(b + o, n - o, "  track %5.1f°  GS %3.0f kt", f5.trk, f5.gs);
            if (f5.tas_ok) o += snprintf(b + o, n - o, "  TAS %3.0f kt", f5.tas);
            if (f5.rate_ok) snprintf(b + o, n - o, "  turn %+.2f°/s", f5.rate);
            tail = nochk;
            break;
        case 0x51:
            snprintf(b, n, " %9.4f %9.4f  alt %6.0f ft", f51.lat, f51.lon, f51.alt);
            if (area51) tail = "  (no ADS-B position: area check)";
            break;
        case 0x52:
            if (f52.resolved) o = snprintf(b, n, " %11.6f %11.6f", f52.lat, f52.lon);
            else o = snprintf(b, n, " fine %+.6f %+.6f (mod %.4f°, no position known)", f52.fla, f52.flo, B52_SPAN);
            snprintf(b + o, n - o, "  %s %6.0f ft  %s", f52.gnss ? "GNSS ht" : "alt", f52.alt, fom_str(f52.fom));
            break;
        case 0x60:
            o = snprintf(b, n, " hdg %5.1f°", f6.hdg);
            if (f6.ias_ok) o += snprintf(b + o, n - o, "  IAS %3.0f kt", f6.ias);
            if (f6.mach_ok) o += snprintf(b + o, n - o, "  M%.3f", f6.mach);
            if (f6.baro_ok) o += snprintf(b + o, n - o, "  baro %+5.0f ft/min", f6.baro);
            if (f6.ins_ok) snprintf(b + o, n - o, "  inertial %+5.0f ft/min", f6.ins);
            tail = nochk;
            break;
        }
        printf("%s  %s%s\n", head, body, tail);
    }
    fflush(stdout);
}

/* rssi: received level from the Beast signal byte, dBFS (NAN when the receiver gives none) */
static void frame(const uint8_t *m, int len, double t, double rssi)
{
    int df = m[0] >> 3;
    if (len == 14 && (df == 17) && crc(m, 112) == 0) df17(m, t);
    else if (len == 14 && (df == 20 || df == 21)) commb(m, df, t, rssi);
    else if (len == 7 && (df == 0 || df == 4)) {
        uint32_t addr = crc(m, 56);
        struct ac *a = ac_get(addr, false);
        if (!a && hex_given(addr)) a = ac_get(addr, true);
        if (a && (t - a->seen <= KNOWN_AGE || hex_given(addr))) set_alt(a, ac13_ft(getbits(m, 20, 32)), t);
    }
}

/* ---------------- Beast stream ---------------- */

static int dial(const char *host, const char *port)
{
    struct addrinfo hints = { .ai_socktype = SOCK_STREAM }, *res;
    if (getaddrinfo(host, port, &hints, &res)) return -1;
    int fd = -1;
    for (struct addrinfo *r = res; r; r = r->ai_next) {
        fd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);
        if (fd < 0) continue;
        if (!connect(fd, r->ai_addr, r->ai_addrlen)) break;
        close(fd);
        fd = -1;
    }
    freeaddrinfo(res);
    return fd;
}

/* frame time: the receiver's 12 MHz timestamp, placed on the wall clock by the
 * smallest (wall - timestamp) seen, since frames reach us in batches */
static double frame_time(const uint8_t *ts)
{
    static double off = NAN, last;
    uint64_t v = 0;
    for (int k = 0; k < 6; k++) v = v << 8 | ts[k];
    double w = now_s(), r = v / 12e6;
    if (!v) return w;                                   /* no timestamp */
    if (isnan(off) || r < last - 1 || w - r < off) off = w - r;      /* first, receiver restarted, or less delay */
    last = r;
    return r + off;
}

/* Beast: 0x1a, type ('1' A/C 2 bytes, '2' short 7, '3' long 14), 6-byte timestamp,
 * 1 signal byte, the message; 0x1a inside is doubled */
static void beast(int fd)
{
    uint8_t buf[4096], f[32];
    int st = 0, need = 0, got = 0;
    bool esc = false;
    double last_stat = now_s();
    for (;;) {
        ssize_t n = read(fd, buf, sizeof(buf));
        if (n <= 0) return;
        for (ssize_t i = 0; i < n; i++) {
            uint8_t c = buf[i];
            if (st == 0) { if (c == 0x1a) st = 1; continue; }
            if (st == 1) {
                need = c == '1' ? 2 : c == '2' ? 7 : c == '3' ? 14 : 0;
                if (!need) { st = c == 0x1a ? 1 : 0; continue; }
                need += 7; got = 0; esc = false; st = 2;
                continue;
            }
            if (c == 0x1a && !esc) { esc = true; continue; }
            if (esc && c != 0x1a) { st = 1; i--; esc = false; continue; }   /* unescaped 0x1a: a new frame */
            esc = false;
            f[got++] = c;
            if (got == need) {
                /* signal byte: sqrt of the signal power (full scale 1) x 255, as dump1090 writes it */
                if (need > 9) frame(f + 7, need - 7, frame_time(f), f[6] ? 20 * log10(f[6] / 255.0) : NAN);
                st = 0;
            }
        }
        double t = now_s();
        if (verbose && t - last_stat > 60) {
            last_stat = t;
            char regs[256];
            int o = 0;
            for (int k = 0; k < NREG; k++)
                o += snprintf(regs + o, sizeof(regs) - o, "BDS%X,%X %ld, ", reg_code[k] >> 4, reg_code[k] & 15, n_reg[k]);
            fprintf(stderr, "# Comm-B %ld: unknown aircraft %ld, data parity %ld (other registers %ld), %sabove --max-alt %ld, other aircraft (-A) %ld, ambiguous %ld, neither %ld\n",
                    n_commb, n_unknown, n_dp, n_dp_other, regs, n_high, n_other, n_both, n_none);
        }
    }
}

static void usage(void)
{
    fprintf(stderr,
            "usage: bds_logger [HOST:PORT] [--qnh HPA] [--max-alt FT] [-A LIST] [--bds LIST] [-v]\n"
            "  HOST:PORT   Beast output (default localhost:30005)\n"
            "  --qnh HPA   QNH for AMSL altitudes; default: from BDS 4,0 settings of aircraft below %d ft,\n"
            "              or the pressure altitude (STD) until one is seen\n"
            "  --max-alt FT  hide aircraft above this altitude: AMSL when a QNH is known (--qnh or\n"
            "              derived), else pressure altitude; also those with no altitude known\n"
            "  -A LIST     only these aircraft: hex addresses or callsigns, separated by commas\n"
            "              (e.g. -A 48af10,LOT3TD; a callsign matches once its ADS-B identification or a\n"
            "              BDS 2,0 reply with it is heard).\n"
            "              An aircraft given as hex is decoded also without ADS-B: then all fields of the\n"
            "              register must be present, and the line ends '(no ADS-B check)'\n"
            "  -v          statistics every minute on stderr\n"
            "  --bds LIST  registers to print, of 10,20,30,40,44,45,50,51,52,60 (default: all but 50);\n"
            "              the others are still checked, so a reply that also fits one of them is marked '?'.\n"
            "              44 and 45 come only from data parity replies\n"
            "Prints Comm-B replies of ADS-B aircraft, one per line: BDS 1,0 capability, 2,0 identification,\n"
            "3,0 ACAS RA, 4,0 selected altitude, 4,4 / 4,5 meteo, 5,0 track and turn, 5,1 / 5,2 position,\n"
            "6,0 heading and speed.\n"
            "Replies with data parity (address XOR register) name their register: marked 'D'.\n", QNH_MAX_FT);
}

int main(int argc, char **argv)
{
    const char *hp = "localhost:30005";
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--qnh") && i + 1 < argc) {
            qnh_cmd = atof(argv[++i]);
            if (qnh_cmd < 900 || qnh_cmd > 1100) { fprintf(stderr, "--qnh: %s is not a QNH in hPa\n", argv[i]); return 1; }
        } else if (!strcmp(argv[i], "--max-alt") && i + 1 < argc) {
            char *e;
            long v = strtol(argv[++i], &e, 10);
            if (*e || v < -1000 || v > 60000) { fprintf(stderr, "--max-alt: %s is not an altitude in ft\n", argv[i]); return 1; }
            max_alt = (int)v;
        } else if (!strcmp(argv[i], "-A") && i + 1 < argc) {
            char list[512];
            snprintf(list, sizeof(list), "%s", argv[++i]);
            for (char *sv, *tok = strtok_r(list, ", ", &sv); tok; tok = strtok_r(NULL, ", ", &sv)) {
                if (nsel == SEL_N || strlen(tok) > 8) { fprintf(stderr, "-A: too many entries or '%s' too long\n", tok); return 1; }
                snprintf(sel[nsel++], sizeof(sel[0]), "%s", tok);
            }
        } else if (!strcmp(argv[i], "--bds") && i + 1 < argc) {
            char list[128];
            snprintf(list, sizeof(list), "%s", argv[++i]);
            memset(show, 0, sizeof(show));
            for (char *sv, *tok = strtok_r(list, ", ", &sv); tok; tok = strtok_r(NULL, ", ", &sv)) {
                char *e;
                unsigned long v = strtoul(tok, &e, 16);
                int k = 0;
                while (k < NREG && (*e || reg_code[k] != v)) k++;
                if (k == NREG) { fprintf(stderr, "--bds: '%s' is not one of 10,20,30,40,44,45,50,51,52,60\n", tok); return 1; }
                show[k] = true;
            }
        } else if (!strcmp(argv[i], "-v")) verbose = true;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(); return 0; }
        else if (argv[i][0] != '-') hp = argv[i];
        else { usage(); return 1; }
    }
    char host[256], port[16] = "30005";
    snprintf(host, sizeof(host), "%s", hp);
    char *c = strrchr(host, ':');
    if (c) { snprintf(port, sizeof(port), "%s", c + 1); *c = 0; }
    for (;;) {
        int fd = dial(host, port);
        if (fd < 0) { fprintf(stderr, "# cannot connect to %s:%s, retrying in 5 s\n", host, port); sleep(5); continue; }
        fprintf(stderr, "# connected to %s:%s\n", host, port);
        beast(fd);
        close(fd);
        fprintf(stderr, "# connection to %s:%s lost, reconnecting\n", host, port);
        sleep(1);
    }
}
