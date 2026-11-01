/* Comm-B register formatters shared by ssr1030 (modes.c) and bds_logger.
 * The includer defines getbits(m, first, last) (bits numbered from 1),
 * sgn(m, signbit, first, last) and ac13_ft(13-bit altitude code). MB bit k
 * is message bit k + 32. Each returns the characters written. */
#ifndef COMMB_FMT_H
#define COMMB_FMT_H

static const char *const hazard[] = { "nil", "light", "moderate", "severe" };

/* BDS 1,0 data link capability, in groups (see docs/bds10.md):
 * transponder (OVC, subnetwork version, level 5, specific services), ELM
 * throughput, surveillance (ident, SCS, SI codes, 1,7 changed toggle),
 * ACAS (operating, hybrid surveillance, RA, version), DTE subaddresses */
static int fmt10(const uint8_t *m, char *buf, size_t n)
{
    static const char *const rtca[] = { "DO-185", "DO-185A", "DO-185B", "DO-185?" };
#define FL(bit, name) (getbits(m, bit, bit) ? " " name : "")
    int o = snprintf(buf, n, "%s subnet v%u%s%s", FL(47, "OVC"), getbits(m, 49, 55), FL(56, "L5"), FL(57, "SSS"));
    o += snprintf(buf + o, n - o, "  ELM up %u down %u", getbits(m, 58, 60), getbits(m, 61, 64));
    o += snprintf(buf + o, n - o, " %s%s%s%s", FL(65, "ident"), FL(66, "SCS"), FL(67, "SI"), FL(68, "1,7-chg"));
    if (getbits(m, 48, 48))
        o += snprintf(buf + o, n - o, "  ACAS%s%s %s", FL(69, "HS"), FL(70, "RA"), rtca[getbits(m, 71, 72)]);
    else o += snprintf(buf + o, n - o, "  no ACAS");
    o += snprintf(buf + o, n - o, "  DTE %04x%s", getbits(m, 73, 88), FL(41, " +cont"));
#undef FL
    return o;
}

/* BDS 3,0 ACAS active resolution advisory: ARA (MB 9-22), RAC (23-26), RA
 * terminated, multiple threat, threat type and identity (address, or
 * altitude, range and bearing) */
static int fmt30(const uint8_t *m, char *buf, size_t n)
{
    int o;
    if (getbits(m, 41, 41))
        o = snprintf(buf, n, " RA %s %s%s%s%s%s", getbits(m, 42, 42) ? "corrective" : "preventive",
                     getbits(m, 43, 43) ? "down" : "up", getbits(m, 44, 44) ? " increase" : "",
                     getbits(m, 45, 45) ? " reversal" : "", getbits(m, 46, 46) ? " crossing" : "",
                     getbits(m, 47, 47) ? " positive" : "");
    else o = snprintf(buf, n, " no RA");
    static const char *const rac[] = { " no-below", " no-above", " no-left", " no-right" };
    for (int k = 0; k < 4; k++) if (getbits(m, 55 + k, 55 + k)) o += snprintf(buf + o, n - o, "%s", rac[k]);
    if (getbits(m, 59, 59)) o += snprintf(buf + o, n - o, " terminated");
    if (getbits(m, 60, 60)) o += snprintf(buf + o, n - o, " multi-threat");
    unsigned tti = getbits(m, 61, 62);
    if (tti == 1) o += snprintf(buf + o, n - o, "  threat %06x", getbits(m, 63, 86));
    else if (tti == 2) {
        int ft = ac13_ft(getbits(m, 63, 75));
        unsigned rg = getbits(m, 76, 82), br = getbits(m, 83, 88);
        o += snprintf(buf + o, n - o, "  threat");
        if (ft != INT32_MIN) o += snprintf(buf + o, n - o, " %d ft", ft);
        if (rg) o += snprintf(buf + o, n - o, " %.1f NM", (rg - 1) / 10.0);
        if (br) o += snprintf(buf + o, n - o, " %u°", 6 * (br - 1) + 3);
    }
    return o;
}

/* BDS 4,4 meteorological routine air report (Doc 9871 layout as pyModeS
 * reads it): FOM/source, wind, static air temperature (sign + 10 bits,
 * 0.25 C, no status bit), average static pressure, turbulence, humidity */
static double sat44(const uint8_t *m) { return sgn(m, 56, 57, 66) * 0.25; }
static int fmt44(const uint8_t *m, char *buf, size_t n)
{
    int o = snprintf(buf, n, " FOM %u", getbits(m, 33, 36));
    if (getbits(m, 37, 37)) o += snprintf(buf + o, n - o, "  wind %u kt %.0f°", getbits(m, 38, 46), getbits(m, 47, 55) * 180.0 / 256);
    o += snprintf(buf + o, n - o, "  SAT %+.2f°C", sat44(m));
    if (getbits(m, 67, 67)) o += snprintf(buf + o, n - o, "  p %u hPa", getbits(m, 68, 78));
    if (getbits(m, 79, 79)) o += snprintf(buf + o, n - o, "  turb %s", hazard[getbits(m, 80, 81)]);
    if (getbits(m, 82, 82)) o += snprintf(buf + o, n - o, "  RH %.0f%%", getbits(m, 83, 88) * 100.0 / 64);
    return o;
}

/* BDS 4,5 meteorological hazard report: turbulence, wind shear, microburst,
 * icing, wake vortex (nil..severe), temperature, pressure, radio height */
static int fmt45(const uint8_t *m, char *buf, size_t n)
{
    static const char *const hz[] = { "turb", "shear", "microburst", "icing", "wake" };
    int o = 0;
    buf[0] = 0;
    for (int k = 0; k < 5; k++)
        if (getbits(m, 33 + 3 * k, 33 + 3 * k)) o += snprintf(buf + o, n - o, " %s %s", hz[k], hazard[getbits(m, 34 + 3 * k, 35 + 3 * k)]);
    if (getbits(m, 48, 48)) o += snprintf(buf + o, n - o, "  SAT %+.2f°C", sgn(m, 49, 50, 58) * 0.25);
    if (getbits(m, 59, 59)) o += snprintf(buf + o, n - o, "  p %u hPa", getbits(m, 60, 70));
    if (getbits(m, 71, 71)) o += snprintf(buf + o, n - o, "  radio ht %u ft", getbits(m, 72, 83) * 16);
    if (!o) o = snprintf(buf, n, " (no fields)");
    return o;
}

#endif
