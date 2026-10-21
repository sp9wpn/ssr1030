/* Aircraft positions from ADS-B (DF17 airborne position, CPR local decoding
 * against a reference point) and the radar -> aircraft -> receiver path
 * geometry used to predict reply delays. */
#include <math.h>
#include "ssr1030.h"

#define NZ      15
#define CPR_MAX 131072.0                 /* 2^17 */
#define A_WGS   6378137.0
#define E2_WGS  6.69437999014e-3

/* number of longitude zones at this latitude (1090-WP-9-14, 1.7.2) */
static int cpr_nl(double lat)
{
    lat = fabs(lat);
    if (lat < 1e-9) return 59;
    if (lat > 87.0) return lat > 87.0 + 1e-9 ? 1 : 2;
    double a = 1 - cos(M_PI / (2 * NZ)), c = cos(M_PI / 180 * lat);
    return (int)floor(2 * M_PI / acos(1 - a / (c * c)));
}

static double modpos(double x, double y) { double r = fmod(x, y); return r < 0 ? r + y : r; }

bool df17_position(const uint8_t *msg, double ref_lat, double ref_lon, double *lat, double *lon, double *alt_m)
{
    if ((msg[0] >> 3) != 17 && (msg[0] >> 3) != 18) return false;
    int tc = getbits(msg, 33, 37);
    bool baro = tc >= 9 && tc <= 18, gnss = tc >= 20 && tc <= 22;
    if (!baro && !gnss) return false;
    unsigned alt = getbits(msg, 41, 52);
    if (baro) {
        if (!(alt & 0x10)) return false;            /* Q = 0: Gillham code, not decoded */
        unsigned n = (alt >> 5) << 4 | (alt & 0xF);
        *alt_m = (n * 25.0 - 1000.0) * 0.3048;
    } else *alt_m = alt;                            /* GNSS height, m */
    if (!alt) return false;

    int odd = getbits(msg, 54, 54);
    double yz = getbits(msg, 55, 71) / CPR_MAX, xz = getbits(msg, 72, 88) / CPR_MAX;
    double dlat = 360.0 / (4 * NZ - odd);
    double j = floor(ref_lat / dlat) + floor(0.5 + modpos(ref_lat, dlat) / dlat - yz);
    double la = dlat * (j + yz);
    if (la < -90 || la > 90) return false;
    int nl = cpr_nl(la) - odd;
    double dlon = 360.0 / (nl > 0 ? nl : 1);
    double m = floor(ref_lon / dlon) + floor(0.5 + modpos(ref_lon, dlon) / dlon - xz);
    double lo = dlon * (m + xz);
    if (lo > 180) lo -= 360;
    *lat = la;
    *lon = lo;
    return true;
}

void geo_ecef(double lat, double lon, double alt_m, double xyz[3])
{
    double la = lat * M_PI / 180, lo = lon * M_PI / 180;
    double n = A_WGS / sqrt(1 - E2_WGS * sin(la) * sin(la));
    xyz[0] = (n + alt_m) * cos(la) * cos(lo);
    xyz[1] = (n + alt_m) * cos(la) * sin(lo);
    xyz[2] = (n * (1 - E2_WGS) + alt_m) * sin(la);
}

double geo_dist(const double a[3], const double b[3])
{
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

int geo_parse(const char *s, double *lat, double *lon, double *alt_m)
{
    *alt_m = 0;
    int n = sscanf(s, "%lf,%lf,%lf", lat, lon, alt_m);
    if (n < 2 || fabs(*lat) > 90 || fabs(*lon) > 180) return -1;
    return 0;
}
