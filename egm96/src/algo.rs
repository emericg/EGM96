#![allow(clippy::all)]
#![allow(clippy::needless_return)]

/*
 * Copyright (c) 2006 D.Ineiev <ineiev@yahoo.co.uk>
 * Copyright (c) 2020 Emeric Grange <emeric.grange@gmail.com>
 * Copyright (c) 2025 Micah Chambers <micahc.vt@gmail.com>
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
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgment in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 **/

/*
 * This program is designed for the calculation of a geoid undulation at a point
 * whose latitude and longitude is specified.
 *
 * This program is designed to be used with the constants of EGM96 and those of
 * the WGS84(g873) system. The undulation will refer to the WGS84 ellipsoid.
 *
 * It's designed to use the potential coefficient model EGM96 and a set of
 * spherical harmonic coefficients of a correction term.
 * The correction term is composed of several different components, the primary
 * one being the conversion of a height anomaly to a geoid undulation.
 * The principles of this procedure were initially described in the paper:
 * - use of potential coefficient models for geoid undulation determination using
 * a spherical harmonic representation of the height anomaly/geoid undulation
 * difference by R.H. Rapp, Journal of Geodesy, 1996.
 *
 * This program is a modification of the program described in the following report:
 * - a fortran program for the computation of gravimetric quantities from high
 * degree spherical harmonic expansions, Richard H. Rapp, report 334, Department
 * of Geodetic Science and Surveying, the Ohio State University, Columbus, 1982
 **/

use std::f64::consts::PI;
use std::sync::OnceLock;

use crate::egm96_data::EGM96_DATA;

/****************************************************************************/

// Maximum degree and orders of harmonic coefficients
const NMAX: usize = 360;
const N361: usize = 361;
// Size of correction and harmonic coefficients arrays (361*181)
const COEFFS: usize = 65341;

/***************************************************************************/

/// Compute sine and cosine values for the given longitude
fn dscml(rlon: f64, sinml: &mut [f64; N361 + 1], cosml: &mut [f64; N361 + 1]) {
    let a = rlon.sin();
    let b = rlon.cos();

    sinml[1] = a;
    cosml[1] = b;
    sinml[2] = 2.0 * b * a;
    cosml[2] = 2.0 * b * b - 1.0;

    for m in 3..=NMAX {
        sinml[m] = 2.0 * b * sinml[m - 1] - sinml[m - 2];
        cosml[m] = 2.0 * b * cosml[m - 1] - cosml[m - 2];
    }
}

/// Compute height undulation based on coefficients
fn hundu(
    p: &[f64; COEFFS + 1],
    sinml: &[f64; N361 + 1],
    cosml: &[f64; N361 + 1],
    gr: f64,
    re: f64,
) -> f64 {
    // WGS 84 gravitational constant in m^3/s^2 (mass of Earth's atmosphere included)
    const GM: f64 = 0.3986004418e15;
    // WGS 84 datum surface equatorial radius
    const AE: f64 = 6378137.0;

    let ar = AE / re;
    let mut arn = ar;
    let mut ac = 0.0;
    let mut a = 0.0;

    let mut k = 3;
    for n in 2..=NMAX {
        arn *= ar;
        k += 1;
        let mut sum = p[k] * EGM96_DATA[k][2] as f64;
        let mut sumc = p[k] * EGM96_DATA[k][0] as f64;

        for m in 1..=n {
            k += 1;
            let tempc = EGM96_DATA[k][0] as f64 * cosml[m] + EGM96_DATA[k][1] as f64 * sinml[m];
            let temp = EGM96_DATA[k][2] as f64 * cosml[m] + EGM96_DATA[k][3] as f64 * sinml[m];
            sumc += p[k] * tempc;
            sum += p[k] * temp;
        }
        ac += sumc;
        a += sum * arn;
    }
    ac += EGM96_DATA[1][0] as f64
        + (p[2] * EGM96_DATA[2][0] as f64)
        + (p[3] * (EGM96_DATA[3][0] as f64 * cosml[1] + EGM96_DATA[3][1] as f64 * sinml[1]));

    // Add haco = ac/100 to convert height anomaly on the ellipsoid to the undulation
    // Add -0.53m to make undulation refer to the WGS84 ellipsoid

    ((a * GM) / (gr * re)) + (ac / 100.0) - 0.53
}

fn legfdn(m: usize, theta: f64, rleg: &mut [f64; N361 + 1]) {
    static DRTS_DIRT: OnceLock<(Vec<f64>, Vec<f64>)> = OnceLock::new();
    let (drts, dirt) = DRTS_DIRT.get_or_init(|| {
        let nmax2p = (2 * NMAX) + 1;
        let mut drts = vec![0.0; 1301];
        let mut dirt = vec![0.0; 1301];
        for n in 1..=nmax2p {
            drts[n] = (n as f64).sqrt();
            dirt[n] = 1.0 / drts[n];
        }

        return (drts, dirt);
    });

    let mut rlnn = [0.0; N361 + 1];

    let nmax1 = NMAX + 1;
    let m1 = m + 1;
    let m2 = m + 2;
    let m3 = m + 3;

    let cothet = theta.cos();
    let sithet = theta.sin();

    // compute the legendre functions
    rlnn[1] = 1.0;
    rlnn[2] = sithet * drts[3];
    for n1 in 3..=m1 {
        let n = n1 - 1;
        let n2 = 2 * n;
        rlnn[n1] = drts[n2 + 1] * dirt[n2] * sithet * rlnn[n];
    }

    match m {
        1 => {
            rleg[2] = rlnn[2];
            rleg[3] = drts[5] * cothet * rleg[2];
        }
        0 => {
            rleg[1] = 1.0;
            rleg[2] = cothet * drts[3];
        }
        _ => {}
    }
    rleg[m1] = rlnn[m1];

    if m2 <= nmax1 {
        rleg[m2] = drts[m1 * 2 + 1] * cothet * rleg[m1];
        if m3 <= nmax1 {
            for n1 in m3..=nmax1 {
                let n = n1 - 1;
                if (!m == 0 && n < 2) || (m == 1 && n < 3) {
                    continue;
                }
                let n2 = 2 * n;
                rleg[n1] = drts[n2 + 1]
                    * dirt[n + m]
                    * dirt[n - m]
                    * (drts[n2 - 1] * cothet * rleg[n1 - 1]
                        - drts[n + m - 1] * drts[n - m - 1] * dirt[n2 - 3] * rleg[n1 - 2]);
            }
        }
    }
}

/// Computes geocentric distance, geocentric latitude, and approximate normal gravity
fn radgra(lat: f64, lon: f64, rlat: &mut f64, gr: &mut f64, re: &mut f64) {
    const A: f64 = 6378137.0;
    const E2: f64 = 0.00669437999013;
    const GEQT: f64 = 9.7803253359;
    const K: f64 = 0.00193185265246;
    let t1 = lat.sin().powi(2);
    let n = A / (1.0 - (E2 * t1)).sqrt();
    let t2 = n * lat.cos();
    let x = t2 * lon.cos();
    let y = t2 * lon.sin();
    let z = (n * (1.0 - E2)) * lat.sin();

    *re = (x * x + y * y + z * z).sqrt(); // compute the geocentric radius
    *rlat = (z / (x * x + y * y).sqrt()).atan(); // compute the geocentric latitude
    *gr = GEQT * (1.0 + (K * t1)) / (1.0 - (E2 * t1)).sqrt(); // compute normal gravity (m/sec²)
}

/// Compute the geoid undulation from the EGM96 model
fn undulation(lat: f64, lon: f64) -> f64 {
    let mut p = [0.0; COEFFS + 1];
    let mut sinml = [0.0; N361 + 1];
    let mut cosml = [0.0; N361 + 1];
    let mut rleg = [0.0; N361 + 1];

    let mut rlat = 0.0;
    let mut gr = 0.0;
    let mut re = 0.0;
    let nmax1 = NMAX + 1;

    // Compute the geocentric latitude, geocentric radius, normal gravity
    radgra(lat, lon, &mut rlat, &mut gr, &mut re);
    rlat = (PI / 2.0) - rlat;

    for j in 1..=nmax1 {
        let m = j - 1;
        legfdn(m, rlat, &mut rleg);
        for i in j..=nmax1 {
            p[((i - 1) * i) / 2 + m + 1] = rleg[i];
        }
    }
    dscml(lon, &mut sinml, &mut cosml);

    hundu(&p, &sinml, &cosml, gr, re)
}

fn wrap_degrees(mut degrees: f64) -> f64 {
    degrees += 180.0;
    degrees = degrees.rem_euclid(360.0);
    degrees - 180.0
}

/// Public function to compute altitude offset using EGM96 model
pub fn egm96_compute_altitude_offset(lat: f64, lon: f64) -> f64 {
    let lon = wrap_degrees(lon);
    let lat = lat.clamp(-90.0, 90.0);
    undulation(lat.to_radians(), lon.to_radians())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_radgra_at_equator() {
        let lat = 0.0;
        let lon = 0.0;
        let mut rlat = 0.0;
        let mut gr = 0.0;
        let mut re = 0.0;
        radgra(lat, lon, &mut rlat, &mut gr, &mut re);
        assert!((rlat).abs() < 1e-9); // Geocentric latitude should be close to 0 at equator
        assert!((re - 6378137.0).abs() < 1e-9); // Geocentric radius should be close to semi-major axis
        assert!((gr - 9.7803253359).abs() < 1e-9); // Normal gravity at equator
    }

    #[test]
    fn test_radgra_at_pole() {
        let lat = FRAC_PI_2;
        let lon = 0.0;
        let mut rlat = 0.0;
        let mut gr = 0.0;
        let mut re = 0.0;
        radgra(lat, lon, &mut rlat, &mut gr, &mut re);
        assert!((rlat - FRAC_PI_2).abs() < 1e-9); // Geocentric latitude should be close to pi/2 at pole
        assert!((re - 6356752.314245).abs() < 1e-6); // Geocentric radius at pole (semi-minor axis)
        assert!((gr - 9.8321863685).abs() < 1e-5); // Normal gravity at pole
    }

    #[test]
    fn test_undulation_at_locations() {
        struct Check {
            lat: f64,
            lon: f64,
            geoid: f64,
        }

        let checks = [
            //Houston       :
            Check {
                lat: 29.7604,
                lon: -95.3698,
                geoid: -28.41,
            },
            //San Antonio   :
            Check {
                lat: 29.4241,
                lon: -98.4936,
                geoid: -26.52,
            },
            //San Diego     :
            Check {
                lat: 32.7157,
                lon: -117.1611,
                geoid: -35.22,
            },
            //Dallas        :
            Check {
                lat: 32.7767,
                lon: -96.797,
                geoid: -27.34,
            },
            //San Jose      :
            Check {
                lat: 37.3382,
                lon: -121.8863,
                geoid: -32.37,
            },
            //Los Angeles   :
            Check {
                lat: 34.0522,
                lon: -118.2437,
                geoid: -35.17,
            },
            //New York      :
            Check {
                lat: 40.7128,
                lon: -74.006,
                geoid: -32.73,
            },
            //San Francisco :
            Check {
                lat: 37.7749,
                lon: -122.4194,
                geoid: -32.17,
            },
            //Chicago       :
            Check {
                lat: 41.8781,
                lon: -87.6298,
                geoid: -33.93,
            },
            //London        :
            Check {
                lat: 51.5074,
                lon: 0.1278,
                geoid: 45.78,
            },
            //Paris         :
            Check {
                lat: 48.8566,
                lon: 2.3522,
                geoid: 44.61,
            },
            //Toky          :
            Check {
                lat: 35.6895,
                lon: 139.6917,
                geoid: 36.71,
            },
            //Philadelphia  :
            Check {
                lat: 40.05,
                lon: -75.45,
                geoid: -34.32,
            },
            //Phoenix       :
            Check {
                lat: 33.4484,
                lon: -112.074,
                geoid: -30.25,
            },
            //null island
            Check {
                lat: 0.0,
                lon: 0.0,
                geoid: 17.22,
            },
        ];

        for check in checks {
            let computed = egm96_compute_altitude_offset(check.lat, check.lon);
            let expected = check.geoid;
            let err = (computed - expected).abs();
            if err.is_nan() || err.is_infinite() || err > 0.5 {
                panic!(
                    "Lat: {}, Lon: {}, Expected: {expected}, Computed: {computed}",
                    check.lat, check.lon
                );
            }
        }
    }

    #[test]
    fn test_wrap_degrees() {
        assert_eq!(wrap_degrees(0.0), 0.0);
        assert_eq!(wrap_degrees(179.8), 179.8);
        assert_eq!(wrap_degrees(-179.0), -179.0);
        assert_eq!(wrap_degrees(-181.0), 179.0);
        assert_eq!(wrap_degrees(190.0), -170.0);
        assert_eq!(wrap_degrees(-190.0), 170.0);
        assert_eq!(wrap_degrees(-190.0 - 360.0), 170.0);
        assert_eq!(wrap_degrees(360.0), 0.0);
        assert_eq!(wrap_degrees(540.0), -180.0);
        assert_eq!(wrap_degrees(-540.0), -180.0);
        assert_eq!(wrap_degrees(1000.0), -80.0);
        assert_eq!(wrap_degrees(-1000.0), 80.0);
    }
}
