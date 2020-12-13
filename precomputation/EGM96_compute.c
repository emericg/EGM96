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

#include "../EGM96.h"

#include <stdio.h>
#include <math.h>

/* ************************************************************************** */

#define _coeffs (65341) //!< Size of correction and harmonic coefficients arrays (361*181)
#define _nmax   (360)   //!< Maximum degree and orders of harmonic coefficients.
#define _361    (361)

//! EGM96 correction and harmonic coefficients
static double egm96_data[_coeffs+1][4];

/* ************************************************************************** */

/*!
 * \param f_12: EGM96 coefficients file.
 *
 * The even degree zonal coefficients given below were computed for the WGS84(g873)
 * system of constants and are identical to those values used in the NIMA gridding procedure.
 * Computed using subroutine grs written by N.K. PAVLIS.
 */
void dhcsin(FILE *f_12)
{
    double c, s, ec, es;

    const double j2 = 0.108262982131e-2;
    const double j4 = -.237091120053e-05;
    const double j6 = 0.608346498882e-8;
    const double j8 = -0.142681087920e-10;
    const double j10 = 0.121439275882e-13;

    unsigned n;
    unsigned m = (((_nmax+1) * (_nmax+2)) / 2);
    for (n = 1; n <= m; n++)
    {
        egm96_data[n][2] = egm96_data[n][3] = 0;
    }

    while (6 == fscanf(f_12, "%i %i %lf %lf %lf %lf", &n, &m, &c, &s, &ec, &es))
    {
        if (n > _nmax) continue;
        n = ((n * (n + 1)) / 2) + m + 1;
        egm96_data[n][2] = c;
        egm96_data[n][3] = s;
    }
    egm96_data[4][2] += j2 / sqrt(5);
    egm96_data[11][2] += j4 / 3.0;
    egm96_data[22][2] += j6 / sqrt(13);
    egm96_data[37][2] += j8 / sqrt(17);
    egm96_data[56][2] += j10 / sqrt(21);
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
            egm96_data[i][0] = egm96_data[i][1] = 0;
        }

        while (4 == fscanf(f_1, "%i %i %lg %lg", &n, &m, &t1, &t2))
        {
            ig = ((n * (n+1)) / 2) + m + 1;
            egm96_data[ig][0] = t1;
            egm96_data[ig][1] = t2;
        }

        // the correction coefficients are now read in
        // the potential coefficients are now read in,
        // and the reference even degree zonal harmonic coefficients removed to degree 6
        dhcsin(f_12);

        fclose(f_1);
        fclose(f_12);
    }
}

/*!
 * \brief Write precomputed EGM96 correction and harmonic coefficients to egm96_data.h
 */
void write_arrays()
{
    FILE *precomp_out = fopen("EGM96_data.h", "wb");

    if (precomp_out)
    {
        fprintf(precomp_out, "#ifndef EGM96_DATA_H\n");
        fprintf(precomp_out, "#define EGM96_DATA_H\n\n");
        fprintf(precomp_out, "//! Precomputed EGM96 correction and harmonic coefficients\n");
        fprintf(precomp_out, "static const double egm96_data[65342][4] = {\n");

        for (unsigned i = 0; i <= _coeffs; i++)
        {
            fprintf(precomp_out, "{%g,%g,%g,%g},\n", egm96_data[i][0], egm96_data[i][1], egm96_data[i][2], egm96_data[i][3]);
        }

        fprintf(precomp_out, "};\n\n");
        fprintf(precomp_out, "#endif // EGM96_DATA_H\n");

        fclose(precomp_out);
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
 * - computed geoid heights ("OUTPUT.dat") => unit = 20
 * - precomputed egm96_data.h (to use with the library)
 */
int main(void)
{
    init_arrays();
    write_arrays();

    FILE *f_14 = fopen("INPUT.DAT", "rb");
    FILE *f_20 = fopen("OUTPUT.DAT", "wb");

    if (f_14 && f_20)
    {
        // read geodetic latitude,longitude at point undulation is wanted
        double flat, flon;
        while (2 == fscanf(f_14, "%lg %lg", &flat, &flon))
        {
            // compute the geocentric latitude, geocentric radius, normal gravity
            double u = egm96_compute_altitude_offset(flat, flon);

            // u is the geoid undulation from the EGM96 potential coefficient model
            // including the height anomaly to geoid undulation correction term
            // and a correction term to have the undulations refer to the
            // WGS84 ellipsoid. the geoid undulation unit is meters.

            fprintf(f_20, "%14.7f %14.7f %16.7f\n", flat, flon, u);
        }

        fclose(f_14);
        fclose(f_20);
    }

    return 0;
}

/* ************************************************************************** */
