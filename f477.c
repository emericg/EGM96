/*
 * Copyright (c) 2006 D.Ineiev <ineiev@yahoo.co.uk>
 * Copyright (c) 2020 Emeric Grange <emeric.grange@gmail.com>
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from
 * the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

/*
 this program is designed for the calculation of a geoid undulation
 at a point whose latitude and longitude is specified. the program
 is designed to use the potential coefficient model egm96 and a
 set of spherical harmonic coefficients of a correction term.
 the correction term is composed of several different components
 the primary one being the conversion of a height anomaly to a geoid
 undulation. the principles of this procedure were initially
 described in the paper: use of potential coefficient models for geoid
 undulation determination using a spherical harmonic representation
 of the height anomaly/geoid undulation difference by R.H. Rapp,
 Journal of Geodesy, 1996.
 this program is designed to be used with the constants of egm96
 and those of the wgs84(g873) system. the undulation will refer to
 the wgs84 ellipsoid.
 specific details on the undulation computation will be found in the
 joint project report describing the development of egm96.
 this program is a modification of the program described in the
 following report:
 a fortran program for the computation of gravimetric quantities from
 high degree spherical harmonic expansions, Richard H. Rapp,
 report 334, Department of Geodetic Science and Surveying, the Ohio
 State University, Columbus, 1982
 this program was put in this form in Dec 1996.
 rhrapp.f477.nonly

dimensions of p,q,hc,hs must be at least ((maxn+1)*(maxn+2))/2,
dimensions of sinml,cosml must be at least maxn, where maxn is maximum order of computation
*/

#include <stdio.h>
#include <math.h>

/* ************************************************************************** */

#define _coeffs (65341) //!< Size of correction and harmonic coefficients arrays (361*181)
#define _nmax   (360)   //!< Maximum degree and orders of harmonic coefficients.
#define _361    (361)

static double cc[_coeffs+1], cs[_coeffs+1], hc[_coeffs+1], hs[_coeffs+1];
static double p[_coeffs+1], sinml[_361+1], cosml[_361+1], rleg[_361+1];

/* ************************************************************************** */

double hundu(double p[_coeffs+1],
             double hc[_coeffs+1], double hs[_coeffs+1],
             double sinml[_361+1], double cosml[_361+1],
             double gr, double re,
             double cc[_coeffs+1], double cs[_coeffs+1])
{
    // WGS 84 gravitational constant in m³/s² (mass of Earth’s atmosphere included)
    const double GM = 0.3986004418e15;
    // WGS 84 datum surface equatorial radius
    const double ae = 6378137.0;

    double ar = ae/re;
    double arn = ar;
    double ac = 0;
    double a = 0;

    unsigned k = 3;
    for (unsigned n = 2; n <= _nmax; n++)
    {
        arn *= ar;
        k++;
        double sum = p[k]*hc[k];
        double sumc = p[k]*cc[k];

        for (unsigned m = 1; m <= n; m++)
        {
            k++;
            double tempc = cc[k]*cosml[m] + cs[k]*sinml[m];
            double temp  = hc[k]*cosml[m] + hs[k]*sinml[m];
            sumc += p[k]*tempc;
            sum  += p[k]*temp;
        }
        ac += sumc;
        a += sum*arn;
    }
    ac += cc[1] + (p[2]*cc[2]) + (p[3] * (cc[3]*cosml[1] + cs[3]*sinml[1]));

    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
    // Add -0.53m to make undulation refer to the WGS84 ellipsoid

    return ((a * GM) / (gr * re)) + (ac / 100.0) - 0.53;
}

void dscml(double rlon, double sinml[_361+1], double cosml[_361+1])
{
    double a = sin(rlon);
    double b = cos(rlon);

    sinml[1] = a;
    cosml[1] = b;
    sinml[2] = 2*b*a;
    cosml[2] = 2*b*b - 1;

    for (unsigned m = 3; m <= _nmax; m++)
    {
        sinml[m] = 2*b*sinml[m-1] - sinml[m-2];
        cosml[m] = 2*b*cosml[m-1] - cosml[m-2];
    }
}

void dhcsin(FILE *f_12, double hc[_coeffs+1], double hs[_coeffs+1])
{
    double c, s, ec, es;

    // The even degree zonal coefficients given below were computed for the WGS84(g873)
    // system of constants and are identical to those values used in the NIMA gridding procedure.
    // Computed using subroutine grs written by N.K. PAVLIS

    const double j2 = 0.108262982131e-2;
    const double j4 = -.237091120053e-05;
    const double j6 = 0.608346498882e-8;
    const double j8 = -0.142681087920e-10;
    const double j10 = 0.121439275882e-13;

    unsigned n;
    unsigned m = (((_nmax+1) * (_nmax+2)) / 2);
    for (n = 1; n <= m; n++)
    {
        hc[n] = hs[n] = 0;
    }

    while (6 == fscanf(f_12, "%i %i %lf %lf %lf %lf", &n, &m, &c, &s, &ec, &es))
    {
        if (n > _nmax) continue;
        n = ((n * (n + 1)) / 2) + m + 1;
        hc[n] = c;
        hs[n] = s;
    }
    hc[4]  += j2 / sqrt(5);
    hc[11] += j4 / 3.0;
    hc[22] += j6 / sqrt(13);
    hc[37] += j8 / sqrt(17);
    hc[56] += j10 / sqrt(21);
}

/*!
 * \param m: order.
 * \param theta: Colatitude (radians).
 * \param rleg: Normalized legendre function.
 *
 * This subroutine computes all normalized legendre function in 'rleg'.
 * All calculations are in double precision.
 *
 * 'ir' must be set to zero before the first call to this sub.
 * the dimensions of arrays 'rleg' must be at least equal to nmax+1.
 *
 * Original programmer: Oscar L. Colombo, Dept. of Geodetic Science the Ohio State University, August 1980.
 * ineiev: I removed the derivatives, for they are never computed here.
 */
void legfdn(unsigned m, double theta, double rleg[_361+1])
{
    static double drts[1301], dirt[1301], cothet, sithet, rlnn[_361+1];
    static int ir; // TODO 'ir' must be set to zero before the first call to this sub.

    unsigned nmax1 = _nmax + 1;
    unsigned nmax2p = (2 * _nmax) + 1;
    unsigned m1 = m + 1;
    unsigned m2 = m + 2;
    unsigned m3 = m + 3;
    unsigned n, n1, n2;

    if (!ir)
    {
        ir = 1;
        for (n = 1; n <= nmax2p; n++)
        {
            drts[n] = sqrt(n);
            dirt[n] = 1 / drts[n];
        }
    }

    cothet = cos(theta);
    sithet = sin(theta);

    // compute the legendre functions
    rlnn[1] = 1;
    rlnn[2] = sithet * drts[3];
    for (n1 = 3; n1 <= m1; n1++)
    {
        n = n1 - 1;
        n2 = 2 * n;
        rlnn[n1] = drts[n2 + 1] * dirt[n2] * sithet * rlnn[n];
    }

    switch (m)
    {
        case 1:
            rleg[2] = rlnn[2];
            rleg[3] = drts[5] * cothet * rleg[2];
            break;
        case 0:
            rleg[1] = 1;
            rleg[2] = cothet * drts[3];
            break;
    }
    rleg[m1] = rlnn[m1];

    if (m2 <= nmax1)
    {
        rleg[m2] = drts[m1*2 + 1] * cothet * rleg[m1];
        if (m3 <= nmax1)
        {
            for (n1 = m3; n1 <= nmax1; n1++)
            {
                n = n1 - 1;
                if ((!m && n < 2) || (m == 1 && n < 3)) continue;
                n2 = 2 * n;
                rleg[n1] = drts[n2+1] * dirt[n+m] * dirt[n-m] * (drts[n2-1] * cothet * rleg[n1-1] - drts[n+m-1] * drts[n-m-1] * dirt[n2-3] * rleg[n1-2]);
            }
        }
    }
}

/*!
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \param re: Geocentric radius.
 * \param rlat: Geocentric latitude.
 * \param gr: Normal gravity (m/sec²).
 *
 * This subroutine computes geocentric distance to the point, the geocentric
 * latitude, and an approximate value of normal gravity at the point based the
 * constants of the WGS84(g873) system are used.
 */
void radgra(double lat, double lon, double *rlat, double *gr, double *re)
{
    const double a = 6378137.0;
    const double e2 = 0.00669437999013;
    const double geqt = 9.7803253359;
    const double k = 0.00193185265246;
    double t1 = sin(lat) * sin(lat);
    double n = a / sqrt(1.0 - (e2 * t1));
    double t2 = n * cos(lat);
    double x = t2 * cos(lon);
    double y = t2 * sin(lon);
    double z = (n * (1 - e2)) * sin(lat);

    *re = sqrt((x * x) + (y * y) + (z * z));            // compute the geocentric radius
    *rlat = atan(z / sqrt((x * x) + (y * y)));          // compute the geocentric latitude
    *gr = geqt * (1 + (k * t1)) / sqrt(1 - (e2 * t1));  // compute normal gravity (m/sec²)
}

/*!
 * \param lat: Latitude in radians.
 * \param lon: Longitude in radians.
 * \return The geoid undulation (in meters) from the EGM96 potential coefficient model, for given lat/lon.
 */
double undulation(double lat, double lon)
{
    double rlat, gr, re;
    unsigned nmax1 = _nmax + 1;

    radgra(lat, lon, &rlat, &gr, &re);
    rlat = (M_PI / 2) - rlat;

    for (unsigned j = 1; j <= nmax1; j++)
    {
        unsigned m = j - 1;
        legfdn(m, rlat, rleg);
        for (unsigned i = j ; i <= nmax1; i++)
        {
            p[(((i - 1) * i) / 2) + m + 1] = rleg[i];
        }
     }
     dscml(lon, sinml, cosml);

     return hundu(p, hc, hs, sinml, cosml, gr, re, cc, cs);
}

/* ************************************************************************** */

void init_arrays()
{
    int ig, n, m;
    double t1, t2;

    FILE *f_1 = fopen("CORCOEF", "rb");   // correction coefficient file: modified with 'sed -e"s/D/e/g"' to be read with fscanf
    FILE *f_12 = fopen("EGM96", "rb");    // potential coefficient file

    if (f_1 && f_12)
    {
        for (unsigned i = 1; i <= _coeffs; i++)
        {
            cc[i] = cs[i] = 0;
        }

        while (4 == fscanf(f_1, "%i %i %lg %lg", &n, &m, &t1, &t2))
        {
            ig = ((n * (n+1)) / 2) + m + 1;
            cc[ig] = t1;
            cs[ig] = t2;
        }

        // the correction coefficients are now read in
        // the potential coefficients are now read in,
        // and the reference even degree zonal harmonic coefficients removed to degree 6
        dhcsin(f_12, hc, hs);

        fclose(f_1);
        fclose(f_12);
    }
}

/* ************************************************************************** */

/*!
 * \brief Main function.
 * \return 0 if success.
 *
 * The input files consist of:
 * - correction coefficient set ("CORRCOEF") => unit = 1
 * - potential coefficient set ("EGM96") => unit = 12
 * - points at which to compute ("INPUT.dat") => unit = 14
 * The output file is:
 * - computed geoid heights ("OUTF477") => unit = 20
 */
int main(void)
{
    const double rad = (180.0 / M_PI);

    init_arrays();

    FILE *f_14 = fopen("INPUT.DAT", "rb");
    FILE *f_20 = fopen("OUTF477.DAT", "wb");

    if (f_14 && f_20)
    {
        // read geodetic latitude,longitude at point undulation is wanted
        double flat, flon;
        while (2 == fscanf(f_14, "%lg %lg", &flat, &flon))
        {
            // compute the geocentric latitude, geocentric radius, normal gravity
            double u = undulation(flat/rad, flon/rad);

            // u is the geoid undulation from the EGM96 potential coefficient model
            // including the height anomaly to geoid undulation correction term
            // and a correction term to have the undulations refer to the
            // WGS84 ellipsoid. the geoid undulation unit is meters.

            fprintf(f_20, "%14.7f %14.7f %10.7f\n", flat, flon, u);
        }

        fclose(f_14);
        fclose(f_20);
    }

    return 0;
}

/* ************************************************************************** */
