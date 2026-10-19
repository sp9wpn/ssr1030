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
        o += snprintf(b + o, n - o, "%02X", it->msg[k]);
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
        o = snprintf(buf, len, "%-8s RL=%u AQ=%u", it->addr == 0xFFFFFF ? "FFFFFF" : ac_name(it->addr),
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

const char *seq_token_str(uint16_t t, char *buf, size_t len)
{
    int k = t & 31, x = t >> 5;
    if (!x) {
        snprintf(buf, len, "%s", kind_token(k));
        return buf;
    }
    /* named registers are shown bare (call, selAlt, hdgSpd...), others as S-bXY */
    char hx[8];
    const char *a = bds_abbrev(x - 1);
    bool named = a != NULL;
    if (!a) { snprintf(hx, sizeof(hx), "b%X%X", (x - 1) >> 4, (x - 1) & 15); a = hx; }
    if (k == K_UF20 || k == K_UF21) snprintf(buf, len, "S-msg/%s", a);
    else snprintf(buf, len, named ? "%s" : "S-%s", a);
    return buf;
}
