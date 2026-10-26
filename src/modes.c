/* Mode S uplink format helpers: CRC, address/parity recovery, field decode. */
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ssr1030.h"

#define GEN 0x1FFF409u   /* x^24 + ... + 1 */

unsigned getbits(const uint8_t *msg, int first, int last)
{
    unsigned v = 0;
    for (int b = first; b <= last; b++)
        v = (v << 1) | ((msg[(b - 1) / 8] >> (7 - (b - 1) % 8)) & 1);
    return v;
}

uint32_t modes_crc(const uint8_t *msg, int nbits)
{
    uint32_t r = 0;
    int nd = nbits - 24;
    for (int b = 0; b < nd + 24; b++) {
        int bit = b < nd ? (msg[b / 8] >> (7 - b % 8)) & 1 : 0;
        r = (r << 1) | bit;
        if (r & (1u << 24)) r ^= GEN;
    }
    return r & 0xFFFFFF;
}

/* Uplink AP = parity XOR high-order 24 bits of A(x)G(x) (Gertz, MIT LL
 * ATC-117, eq. 5-1). Since the low-order part has degree < 24 it does not
 * change the quotient, so A(x) = quotient of (parity XOR AP)(x) x^24 / G(x).
 * Verified on air: UF11 and UF16 broadcast decode to FFFFFF. */
uint32_t modes_ap_to_addr(const uint8_t *msg, int nbits)
{
    uint32_t ap = getbits(msg, nbits - 23, nbits);
    uint64_t num = (uint64_t)(modes_crc(msg, nbits) ^ ap) << 24;
    uint32_t q = 0;
    for (int deg = 47; deg >= 24; deg--) {
        if ((num >> deg) & 1) {
            num ^= (uint64_t)GEN << (deg - 24);
            q |= 1u << (deg - 24);
        }
    }
    return q;
}

int modes_uf(const uint8_t *msg)
{
    int uf = msg[0] >> 3;
    if ((uf & 0x18) == 0x18) uf = 24;  /* ELM: only the first two bits are UF */
    return uf;
}

int modes_uf_bits(int uf)
{
    return (uf == 16 || uf == 20 || uf == 21 || uf == 24 || uf == 19) ? 112 : 56;
}

int modes_iisi(const uint8_t *msg, int *is_si)
{
    int ic = getbits(msg, 10, 13), cl = getbits(msg, 14, 16);
    if (cl == 0) { *is_si = 0; return ic; }
    *is_si = 1;
    return (cl >= 1 && cl <= 4) ? ic + 16 * (cl - 1) : -1;
}

int modes_own_code(const struct interrog *it)
{
    const uint8_t *m = it->msg;
    switch (it->kind) {
    case K_UF11: {
        int si, c = modes_iisi(m, &si);
        if (c < 0) return -1;
        return si ? (c >= 1 ? 16 + c : -1) : c;
    }
    case K_UF4: case K_UF5: case K_UF20: case K_UF21: {
        int di = getbits(m, 14, 16);
        if (di == 3) { int sis = getbits(m, 17, 22); return sis ? 16 + sis : -1; }
        if (di == 0 || di == 1 || di == 7) return getbits(m, 17, 20);
        return -1;
    }
    default:
        return -1;
    }
}

int radar_code_parse(const char *s)
{
    if ((s[0] == 'S' || s[0] == 's') && (s[1] == 'I' || s[1] == 'i')) {
        int v = atoi(s + 2);
        return v >= 1 && v <= 63 ? 16 + v : -1;
    }
    if ((s[0] == 'I' || s[0] == 'i') && (s[1] == 'I' || s[1] == 'i')) {
        int v = atoi(s + 2);
        return v >= 0 && v <= 15 ? v : -1;
    }
    return -1;
}

void radar_code_str(int code, char *buf, size_t len)
{
    if (code < 0) snprintf(buf, len, "?");
    else if (code < 16) snprintf(buf, len, "II%d", code);
    else snprintf(buf, len, "SI%d", (code - 16) & 63);
}

/* UF4/5/20/21: RR >= 16 requests Comm-B; BDS2 comes from RRS when DI = 7 or 3. */
int modes_bds(const uint8_t *msg, int *bds1, int *bds2)
{
    int rr = getbits(msg, 9, 13), di = getbits(msg, 14, 16);
    if (rr < 16) return 0;
    *bds1 = rr - 16;
    *bds2 = 0;
    if (di == 7) *bds2 = getbits(msg, 21, 24);
    else if (di == 3) *bds2 = getbits(msg, 24, 27);
    return 1;
}

static const char *pc_txt[8] = { "", "PC=1 AC-lockout", "PC=2", "PC=3", "PC=4 cancelB",
                                 "PC=5 cancelC", "PC=6 cancelD", "PC=7" };

static int hexmsg(const struct interrog *it, char *b, size_t n)
{
    int o = 0;
    for (int k = 0; k < it->nbits / 8 && (size_t)o + 3 < n; k++)
        o += snprintf(b + o, n - o, "%02x", it->msg[k]);
    return o;
}

void modes_describe(const struct interrog *it, char *buf, size_t len)
{
    const uint8_t *m = it->msg;
    char hx[32];
    hexmsg(it, hx, sizeof(hx));
    int o = 0;
    switch (it->uf) {
    case 4: case 5: case 20: case 21: {
        int pc = getbits(m, 6, 8), rr = getbits(m, 9, 13), di = getbits(m, 14, 16);
        o = snprintf(buf, len, "%-8s RR=%d DI=%d", ac_name(it->addr), rr, di);
        if (di == 0 || di == 1 || di == 7)
            o += snprintf(buf + o, len - o, " IIS=%u", getbits(m, 17, 20));
        if (di == 3)
            o += snprintf(buf + o, len - o, " SIS=%u", getbits(m, 17, 22));
        if (di == 1) {
            unsigned mbs = getbits(m, 21, 22), mes = getbits(m, 23, 25);
            if (mbs) o += snprintf(buf + o, len - o, " MBS=%u", mbs);
            if (mes) o += snprintf(buf + o, len - o, " MES=%u", mes);
        }
        if ((di == 1 || di == 7) && getbits(m, 26, 26)) o += snprintf(buf + o, len - o, " LOS");
        int b1, b2;
        if (modes_bds(m, &b1, &b2)) o += snprintf(buf + o, len - o, " BDS%X,%X", b1, b2);
        if (pc) o += snprintf(buf + o, len - o, " %s", pc_txt[pc]);
        break;
    }
    case 11: {
        int si, code = modes_iisi(m, &si), pr = getbits(m, 6, 9);
        o = snprintf(buf, len, "%s=%d PR=%d%s", si ? "SI" : "II", code, pr,
                     it->addr == 0xFFFFFF ? "" : " AP?");
        break;
    }
    case 0: case 16:
        o = snprintf(buf, len, "%-8s RL=%u AQ=%u", it->addr == 0xFFFFFF ? "ffffff" : ac_name(it->addr),
                     getbits(m, 9, 9), getbits(m, 14, 14));
        break;
    case 24:
        o = snprintf(buf, len, "%-8s RC=%u NC=%u", ac_name(it->addr),
                     getbits(m, 3, 4), getbits(m, 5, 8));
        break;
    default:
        o = snprintf(buf, len, "%s", ac_name(it->addr));
    }
    if (o < 0) o = 0;
    snprintf(buf + o, len > (size_t)o ? len - o : 0, " [%s]", hx);
}

static const char *labels[K_COUNT] = {
    "Mode 1", "Mode 2", "Mode A", "Mode B", "Mode C", "Mode D",
    "A/C-only all-call (A)", "A/C-only all-call (C)",
    "A/C/S all-call (A)", "A/C/S all-call (C)", "Intermode other",
    "UF0  ACAS short", "UF4  surv. altitude", "UF5  surv. identity",
    "UF11 Mode S all-call", "UF16 ACAS long", "UF19 military",
    "UF20 Comm-A altitude", "UF21 Comm-A identity", "UF24 Comm-C ELM",
    "UF other", "Mode S no sync (P5/SLS)", "Mode S no sync (garbled)"
};
static const char *tokens[K_COUNT] = {
    "1", "2", "AA", "B", "CC", "D", "A", "C", "As", "Cs", "I?",
    "S-acas", "alt", "sqk", "S", "S-acasL", "S-mil", "S-msg", "S-msg", "S-elm", "S-?", "S-sls", "S-grb"
};

static const char *names[K_COUNT] = {
    "M1", "M2", "MA", "MB", "MC", "MD", "ACONLY_A", "ACONLY_C", "ACS_A", "ACS_C", "INTERMODE",
    "UF0", "UF4", "UF5", "UF11", "UF16", "UF19", "UF20", "UF21", "UF24", "UF_OTHER", "S_NOSYNC", "S_GARBLE"
};

const char *kind_name(int k) { return (k >= 0 && k < K_COUNT) ? names[k] : "?"; }
const char *kind_label(int k) { return (k >= 0 && k < K_COUNT) ? labels[k] : "?"; }
const char *kind_token(int k) { return (k >= 0 && k < K_COUNT) ? tokens[k] : "?"; }
bool kind_is_modes(int k) { return k >= K_UF0; }
bool kind_is_allcall(int k)
{
    return k == K_ACONLY_A || k == K_ACONLY_C || k == K_ACS_A || k == K_ACS_C || k == K_UF11;
}

/* Comm-B register (BDS) abbreviations for burst sequences */
const char *bds_abbrev(int b)
{
    switch (b) {
    case 0x00: return "aicb";      /* air-initiated Comm-B */
    case 0x05: return "esap";      /* ES airborne position */
    case 0x06: return "essp";      /* ES surface position */
    case 0x07: return "esst";      /* ES status */
    case 0x08: return "esid";      /* ES identification */
    case 0x09: return "esvel";     /* ES velocity */
    case 0x10: return "cap";       /* data link capability */
    case 0x17: return "gicb";      /* common-usage GICB capability */
    case 0x18: case 0x19: case 0x1A: case 0x1B: case 0x1C: return "mscap";
    case 0x1D: case 0x1E: case 0x1F: return "mscap";
    case 0x20: return "acid";     /* aircraft identification (callsign) */
    case 0x21: return "reg";       /* registration */
    case 0x25: return "type";      /* aircraft type */
    case 0x30: return "ra";        /* ACAS resolution advisory */
    case 0x40: return "selAlt";      /* selected vertical intention (selected altitude) */
    case 0x41: case 0x42: case 0x43: return "wpt";
    case 0x44: return "met";       /* meteorological routine report */
    case 0x45: return "haz";       /* meteorological hazard report */
    case 0x48: return "vhf";
    case 0x50: return "trk";       /* track and turn */
    case 0x51: case 0x52: return "pos";
    case 0x53: return "asv";       /* air-referenced state vector */
    case 0x54: case 0x55: case 0x56: return "wpt";
    case 0x5F: return "qspm";      /* quasi-static parameter monitoring */
    case 0x60: return "hdgSpd";      /* heading and speed */
    case 0x65: return "opst";      /* aircraft operational status */
    default: return NULL;
    }
}

/* sequence token: kind in bits 0-4; for UF4/5/20/21 requesting Comm-B,
 * 1 + BDS (BDS1 * 16 + BDS2) in bits 5-13 */
uint16_t seq_token(const struct interrog *it)
{
    int k = it->kind;
    uint16_t t = (uint16_t)k;
    if (k == K_UF4 || k == K_UF5 || k == K_UF20 || k == K_UF21) {
        int b1, b2;
        if (modes_bds(it->msg, &b1, &b2)) t |= (uint16_t)(1 + b1 * 16 + b2) << 5;
    }
    return t;
}

bool bds_human;

const char *seq_token_str(uint16_t t, char *buf, size_t len)
{
    int k = t & 31, x = t >> 5;
    if (!x) {
        snprintf(buf, len, "%s", kind_token(k));
        return buf;
    }
    if (!bds_human) {                       /* uplink format and register: UF4b50 */
        int uf = k == K_UF4 ? 4 : k == K_UF5 ? 5 : k == K_UF20 ? 20 : 21;
        snprintf(buf, len, "UF%db%X%X", uf, (x - 1) >> 4, (x - 1) & 15);
        return buf;
    }
    /* -h: named registers are shown bare (acid, selAlt, hdgSpd...), others as S-bXY */
    char hx[8];
    const char *a = bds_abbrev(x - 1);
    bool named = a != NULL;
    if (!a) { snprintf(hx, sizeof(hx), "b%X%X", (x - 1) >> 4, (x - 1) & 15); a = hx; }
    if (k == K_UF20 || k == K_UF21) snprintf(buf, len, "S-msg/%s", a);
    else snprintf(buf, len, named ? "%s" : "S-%s", a);
    return buf;
}

/* ---------------- downlink (1090 MHz replies) ---------------- */

/* 13-bit identity / altitude field (bits 20-32: C1 A1 C2 A2 C4 A4 X B1 D1 B2 D2 B4 D4)
 * in 0xABCD layout, one octal digit per nibble */
static unsigned abcd13(unsigned v)
{
#define B13(n) ((v >> (13 - (n))) & 1)
    return B13(6) << 14 | B13(4) << 13 | B13(2) << 12 | B13(12) << 10 | B13(10) << 9 | B13(8) << 8 |
           B13(5) << 6 | B13(3) << 5 | B13(1) << 4 | B13(13) << 2 | B13(11) << 1 | B13(9);
#undef B13
}

/* Gillham (Mode C) code in 0xABCD layout -> altitude in 100 ft, INT32_MIN if invalid */
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

/* 13-bit altitude code (AC13) -> feet, INT32_MIN if absent, metric or invalid */
static int ac13_ft(unsigned v)
{
    if (!v || (v & 0x40)) return INT32_MIN;
    if (v & 0x10) return (int)(((v & 0x1F80) >> 2 | (v & 0x20) >> 1 | (v & 0x0F)) * 25) - 1000;
    int h = gillham(abcd13(v));
    return h == INT32_MIN ? INT32_MIN : h * 100;
}

static int alt13(unsigned v, char *buf, size_t n)
{
    if (!v) return snprintf(buf, n, "alt -");
    if (v & 0x40) return snprintf(buf, n, "alt (metric)");
    int ft = ac13_ft(v);
    return ft == INT32_MIN ? snprintf(buf, n, "alt ?") : snprintf(buf, n, "alt %d ft", ft);
}

int reply_alt_ft(const uint8_t *msg, int df)
{
    return df == 4 || df == 20 ? ac13_ft(getbits(msg, 20, 32)) : INT32_MIN;
}

/* signed field with a separate sign bit before it */
static double sgn(const uint8_t *m, int sbit, int first, int last)
{
    int v = (int)getbits(m, first, last);
    return getbits(m, sbit, sbit) ? v - (1 << (last - first + 1)) : v;
}

static const char ais[] = "#ABCDEFGHIJKLMNOPQRSTUVWXYZ##### ###############0123456789######";

/* Comm-B content of the register the query asked for (MB = bits 33-88) */
static void commb(const uint8_t *m, int bds, char *buf, size_t n)
{
    int o = snprintf(buf, n, "  BDS %X,%X:", bds >> 4, bds & 15);
    switch (bds) {
    case 0x20: {
        char cs[9];
        for (int k = 0; k < 8; k++) cs[k] = ais[getbits(m, 41 + 6 * k, 46 + 6 * k)];
        cs[8] = 0;
        for (int k = 7; k >= 0 && cs[k] == ' '; k--) cs[k] = 0;
        snprintf(buf + o, n - o, getbits(m, 33, 40) == 0x20 ? " id %s" : " (not BDS 2,0)", cs);
        return;
    }
    case 0x40:
        if (getbits(m, 33, 33)) o += snprintf(buf + o, n - o, " MCP alt %u ft", getbits(m, 34, 45) * 16);
        if (getbits(m, 46, 46)) o += snprintf(buf + o, n - o, " FMS alt %u ft", getbits(m, 47, 58) * 16);
        if (getbits(m, 59, 59)) snprintf(buf + o, n - o, " QNH %.1f hPa", getbits(m, 60, 71) * 0.1 + 800);
        return;
    case 0x50:
        if (getbits(m, 33, 33)) o += snprintf(buf + o, n - o, " roll %+.1f°", sgn(m, 34, 35, 43) * 45.0 / 256);
        if (getbits(m, 44, 44)) o += snprintf(buf + o, n - o, " track %.1f°", fmod(sgn(m, 45, 46, 55) * 90.0 / 512 + 360, 360));
        if (getbits(m, 56, 56)) o += snprintf(buf + o, n - o, " GS %u kt", getbits(m, 57, 66) * 2);
        if (getbits(m, 67, 67)) o += snprintf(buf + o, n - o, " turn %+.2f°/s", sgn(m, 68, 69, 77) * 8.0 / 256);
        if (getbits(m, 78, 78)) snprintf(buf + o, n - o, " TAS %u kt", getbits(m, 79, 88) * 2);
        return;
    case 0x60:
        if (getbits(m, 33, 33)) o += snprintf(buf + o, n - o, " hdg %.1f°", fmod(sgn(m, 34, 35, 44) * 90.0 / 512 + 360, 360));
        if (getbits(m, 45, 45)) o += snprintf(buf + o, n - o, " IAS %u kt", getbits(m, 46, 55));
        if (getbits(m, 56, 56)) o += snprintf(buf + o, n - o, " M%.3f", getbits(m, 57, 66) * 2.048 / 512);
        if (getbits(m, 67, 67)) o += snprintf(buf + o, n - o, " baro %+.0f ft/min", sgn(m, 68, 69, 77) * 32);
        if (getbits(m, 78, 78)) snprintf(buf + o, n - o, " inertial %+.0f ft/min", sgn(m, 79, 80, 88) * 32);
        return;
    default:
        o += snprintf(buf + o, n - o, " MB ");
        for (int k = 4; k < 11 && o + 3 < (int)n; k++) o += snprintf(buf + o, n - o, "%02x", m[k]);
    }
}

void reply_describe(const struct rx1090_rec *r, char *buf, size_t n)
{
    const uint8_t *m = r->msg;
    int o = 0;
    switch (r->df) {
    case 11: o = snprintf(buf, n, "%06x", getbits(m, 9, 32)); break;   /* only what the reply carries */
    case 4: case 20: o = alt13(getbits(m, 20, 32), buf, n); break;
    case 5: case 21: o = snprintf(buf, n, "squawk %04x", abcd13(getbits(m, 20, 32))); break;
    default: o = snprintf(buf, n, "DF%d", r->df);
    }
    if ((r->df == 20 || r->df == 21) && r->len == 14 && o < (int)n) {
        if (r->req_bds >= 0) commb(m, r->req_bds, buf + o, n - o);
        else snprintf(buf + o, n - o, "  Comm-B (no register requested)");
    }
}
