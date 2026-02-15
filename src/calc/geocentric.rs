use std::marker::PhantomData;

use crate::types::angle::{Angle, Degrees};
use crate::types::julian_day::{JDay, JulianEphemerisCentury, JulianEphemerisMiliennium};

/// Marker type for Geocentric Longitude
#[derive(Debug, Clone, Copy)]
pub struct GeocentricLongitude;

/// Marker type for Geocentric Latitude
#[derive(Debug, Clone, Copy)]
pub struct GeocentricLatitude;

/// Marker type for Nutation in Longitude (Δψ)
#[derive(Debug, Clone, Copy)]
pub struct NutationLongitude;

/// Marker type for Nutation in Obliquity (Δε)
#[derive(Debug, Clone, Copy)]
pub struct NutationObliquity;

/// Marker type for Mean Ecliptic Obliquity (ε₀)
#[derive(Debug, Clone, Copy)]
pub struct MeanEclipticObliquity;

/// Type-safe geocentric coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geocentric<T>(f64, PhantomData<T>);

impl<T> Geocentric<T> {
    /// Create a new geocentric value
    pub fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    /// Get the underlying value
    pub fn value(&self) -> f64 {
        self.0
    }
}

/// Type-safe nutation value (in degrees)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Nutation<T>(f64, PhantomData<T>);

impl<T> Nutation<T> {
    /// Create a new nutation value in degrees
    pub fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    /// Get the value in degrees
    pub fn degrees(&self) -> f64 {
        self.0
    }

    /// Get the value in arcseconds
    pub fn arcseconds(&self) -> f64 {
        self.0 * 3600.0
    }
}

/// Nutation coefficients: [a, b, c, d] for each term
/// Δψ = Σ[(a + b*JCE) * sin(arg)]
/// Δε = Σ[(c + d*JCE) * cos(arg)]
const NUTATION_ABCD: [[f64; 4]; 63] = [
    [-171996.0, -174.2, 92025.0, 8.9],
    [-13187.0, -1.6, 5736.0, -3.1],
    [-2274.0, -0.2, 977.0, -0.5],
    [2062.0, 0.2, -895.0, 0.5],
    [1426.0, -3.4, 54.0, -0.1],
    [712.0, 0.1, -7.0, 0.0],
    [-517.0, 1.2, 224.0, -0.6],
    [-386.0, -0.4, 200.0, 0.0],
    [-301.0, 0.0, 129.0, -0.1],
    [217.0, -0.5, -95.0, 0.3],
    [-158.0, 0.0, 0.0, 0.0],
    [129.0, 0.1, -70.0, 0.0],
    [123.0, 0.0, -53.0, 0.0],
    [63.0, 0.0, 0.0, 0.0],
    [63.0, 0.1, -33.0, 0.0],
    [-59.0, 0.0, 26.0, 0.0],
    [-58.0, -0.1, 32.0, 0.0],
    [-51.0, 0.0, 27.0, 0.0],
    [48.0, 0.0, 0.0, 0.0],
    [46.0, 0.0, -24.0, 0.0],
    [-38.0, 0.0, 16.0, 0.0],
    [-31.0, 0.0, 13.0, 0.0],
    [29.0, 0.0, 0.0, 0.0],
    [29.0, 0.0, -12.0, 0.0],
    [26.0, 0.0, 0.0, 0.0],
    [-22.0, 0.0, 0.0, 0.0],
    [21.0, 0.0, -10.0, 0.0],
    [17.0, -0.1, 0.0, 0.0],
    [16.0, 0.0, -8.0, 0.0],
    [-16.0, 0.1, 7.0, 0.0],
    [-15.0, 0.0, 9.0, 0.0],
    [-13.0, 0.0, 7.0, 0.0],
    [-12.0, 0.0, 6.0, 0.0],
    [11.0, 0.0, 0.0, 0.0],
    [-10.0, 0.0, 5.0, 0.0],
    [-8.0, 0.0, 3.0, 0.0],
    [7.0, 0.0, -3.0, 0.0],
    [-7.0, 0.0, 0.0, 0.0],
    [-7.0, 0.0, 3.0, 0.0],
    [-7.0, 0.0, 3.0, 0.0],
    [6.0, 0.0, 0.0, 0.0],
    [6.0, 0.0, -3.0, 0.0],
    [6.0, 0.0, -3.0, 0.0],
    [-6.0, 0.0, 3.0, 0.0],
    [-6.0, 0.0, 3.0, 0.0],
    [5.0, 0.0, 0.0, 0.0],
    [-5.0, 0.0, 3.0, 0.0],
    [-5.0, 0.0, 3.0, 0.0],
    [-5.0, 0.0, 3.0, 0.0],
    [4.0, 0.0, 0.0, 0.0],
    [4.0, 0.0, 0.0, 0.0],
    [4.0, 0.0, 0.0, 0.0],
    [-4.0, 0.0, 0.0, 0.0],
    [-4.0, 0.0, 0.0, 0.0],
    [-4.0, 0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
    [-3.0, 0.0, 0.0, 0.0],
];

/// Y-terms for nutation calculation: [Y0, Y1, Y2, Y3, Y4]
/// arg = Y0*X0 + Y1*X1 + Y2*X2 + Y3*X3 + Y4*X4
const NUTATION_YTERMS: [[i32; 5]; 63] = [
    [0, 0, 0, 0, 1],
    [-2, 0, 0, 2, 2],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [-2, 1, 0, 2, 2],
    [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2],
    [-2, -1, 0, 2, 2],
    [-2, 0, 1, 0, 0],
    [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2],
    [2, 0, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1],
    [0, 0, 1, 2, 1],
    [-2, 0, 2, 0, 0],
    [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2],
    [0, 0, 2, 2, 2],
    [0, 0, 2, 0, 0],
    [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0],
    [-2, 0, 0, 2, 0],
    [0, 0, -1, 2, 1],
    [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1],
    [-2, 2, 0, 2, 2],
    [0, 1, 0, 0, 1],
    [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1],
    [0, 0, 2, -2, 0],
    [2, 0, -1, 2, 1],
    [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2],
    [-2, 1, 1, 0, 0],
    [0, -1, 0, 2, 2],
    [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0],
    [-2, 0, 2, 2, 2],
    [-2, 0, 1, 2, 1],
    [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1],
    [0, -1, 1, 0, 0],
    [-2, -1, 0, 2, 1],
    [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1],
    [-2, 0, 2, 0, 1],
    [-2, 1, 0, 2, 1],
    [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0],
    [-2, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2],
    [-1, -1, 1, 0, 0],
    [0, 1, 1, 0, 0],
    [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2],
    [0, 0, 3, 2, 2],
    [2, -1, 0, 2, 2],
];

/// Convert heliocentric longitude to geocentric longitude
/// θ = L + 180°
pub fn calculate_geocentric_longitude(heliocentric_longitude: Angle<Degrees>) -> Angle<Degrees> {
    let theta = heliocentric_longitude.value() + 180.0;
    Angle::new(theta.rem_euclid(360.0))
}

/// Convert heliocentric latitude to geocentric latitude
/// β = -B
pub fn calculate_geocentric_latitude(heliocentric_latitude: Angle<Degrees>) -> Angle<Degrees> {
    Angle::new(-heliocentric_latitude.value())
}

/// Calculate mean elongation of the Moon from the Sun (X₀)
/// X₀ = 297.85036 + 445267.111480*JCE - 0.0019142*JCE² + JCE³/189474
fn mean_elongation(jce: JDay<JulianEphemerisCentury>) -> f64 {
    let jce_val = jce.value();
    297.85036 + 445267.111480 * jce_val - 0.0019142 * jce_val.powi(2) + jce_val.powi(3) / 189474.0
}

/// Calculate mean anomaly of the Sun (X₁)
/// X₁ = 357.52772 + 35999.050340*JCE - 0.0001603*JCE² - JCE³/300000
fn mean_anomaly_sun(jce: JDay<JulianEphemerisCentury>) -> f64 {
    let jce_val = jce.value();
    357.52772 + 35999.050340 * jce_val - 0.0001603 * jce_val.powi(2) - jce_val.powi(3) / 300000.0
}

/// Calculate mean anomaly of the Moon (X₂)
/// X₂ = 134.96298 + 477198.867398*JCE + 0.0086972*JCE² + JCE³/56250
fn mean_anomaly_moon(jce: JDay<JulianEphemerisCentury>) -> f64 {
    let jce_val = jce.value();
    134.96298 + 477198.867398 * jce_val + 0.0086972 * jce_val.powi(2) + jce_val.powi(3) / 56250.0
}

/// Calculate Moon's argument of latitude (X₃)
/// X₃ = 93.27191 + 483202.017538*JCE - 0.0036825*JCE² + JCE³/327270
fn moon_argument_latitude(jce: JDay<JulianEphemerisCentury>) -> f64 {
    let jce_val = jce.value();
    93.27191 + 483202.017538 * jce_val - 0.0036825 * jce_val.powi(2) + jce_val.powi(3) / 327270.0
}

/// Calculate longitude of ascending node of Moon's mean orbit (X₄)
/// X₄ = 125.04452 - 1934.136261*JCE + 0.0020708*JCE² + JCE³/450000
fn moon_ascending_longitude(jce: JDay<JulianEphemerisCentury>) -> f64 {
    let jce_val = jce.value();
    125.04452 - 1934.136261 * jce_val + 0.0020708 * jce_val.powi(2) + jce_val.powi(3) / 450000.0
}

/// Calculate nutation in longitude (Δψ) and obliquity (Δε)
/// Returns (Δψ, Δε) in degrees
pub fn calculate_nutation(
    jce: JDay<JulianEphemerisCentury>,
) -> (Nutation<NutationLongitude>, Nutation<NutationObliquity>) {
    let jce_val = jce.value();

    // Calculate mean parameters
    let x0 = mean_elongation(jce);
    let x1 = mean_anomaly_sun(jce);
    let x2 = mean_anomaly_moon(jce);
    let x3 = moon_argument_latitude(jce);
    let x4 = moon_ascending_longitude(jce);

    let mut delta_psi_sum = 0.0;
    let mut delta_eps_sum = 0.0;

    for i in 0..NUTATION_ABCD.len() {
        let [a, b, c, d] = NUTATION_ABCD[i];
        let [y0, y1, y2, y3, y4] = NUTATION_YTERMS[i];

        // Calculate argument: Σ(Yi * Xi) in radians
        let arg_deg =
            y0 as f64 * x0 + y1 as f64 * x1 + y2 as f64 * x2 + y3 as f64 * x3 + y4 as f64 * x4;
        let arg_rad = arg_deg.to_radians();

        // Δψ = Σ[(a + b*JCE) * sin(arg)]
        delta_psi_sum += (a + b * jce_val) * arg_rad.sin();

        // Δε = Σ[(c + d*JCE) * cos(arg)]
        delta_eps_sum += (c + d * jce_val) * arg_rad.cos();
    }

    // Convert from 0.0001 arcseconds to degrees
    let delta_psi = delta_psi_sum / 36000000.0;
    let delta_eps = delta_eps_sum / 36000000.0;

    (Nutation::new(delta_psi), Nutation::new(delta_eps))
}

/// Calculate mean ecliptic obliquity (ε₀) in arcseconds
/// ε₀ = 84381.448 - 4680.93U - 1.55U² + 1999.25U³ - 51.38U⁴ - 249.67U⁵
///      - 39.05U⁶ + 7.12U⁷ + 27.87U⁸ + 5.79U⁹ + 2.45U¹⁰
/// where U = JME/10
pub fn calculate_mean_ecliptic_obliquity(jme: JDay<JulianEphemerisMiliennium>) -> Angle<Degrees> {
    let u = jme.value() / 10.0;

    // Polynomial in arcseconds
    let e0_arcsec = 84381.448 - 4680.93 * u - 1.55 * u.powi(2) + 1999.25 * u.powi(3)
        - 51.38 * u.powi(4)
        - 249.67 * u.powi(5)
        - 39.05 * u.powi(6)
        + 7.12 * u.powi(7)
        + 27.87 * u.powi(8)
        + 5.79 * u.powi(9)
        + 2.45 * u.powi(10);

    // Convert arcseconds to degrees
    Angle::new(e0_arcsec / 3600.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geocentric_longitude_conversion() {
        let helio_lon = Angle::<Degrees>::new(100.0);
        let geo_lon = calculate_geocentric_longitude(helio_lon);

        assert!((geo_lon.value() - 280.0).abs() < 1e-10);
    }

    #[test]
    fn test_geocentric_longitude_wrapping() {
        let helio_lon = Angle::<Degrees>::new(200.0);
        let geo_lon = calculate_geocentric_longitude(helio_lon);

        // 200 + 180 = 380, wraps to 20
        assert!((geo_lon.value() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_geocentric_latitude_conversion() {
        let helio_lat = Angle::<Degrees>::new(5.0);
        let geo_lat = calculate_geocentric_latitude(helio_lat);

        assert!((geo_lat.value() - (-5.0)).abs() < 1e-10);
    }

    #[test]
    fn test_mean_parameters() {
        let jce = JDay::<JulianEphemerisCentury>::new(0.0);

        let x0 = mean_elongation(jce);
        let x1 = mean_anomaly_sun(jce);
        let x2 = mean_anomaly_moon(jce);
        let x3 = moon_argument_latitude(jce);
        let x4 = moon_ascending_longitude(jce);

        // All should be reasonable angles in degrees
        assert!(x0.abs() < 400.0);
        assert!(x1.abs() < 400.0);
        assert!(x2.abs() < 400.0);
        assert!(x3.abs() < 400.0);
        assert!(x4.abs() < 400.0);
    }

    #[test]
    fn test_nutation_calculation() {
        let jce = JDay::<JulianEphemerisCentury>::new(0.0);

        let (delta_psi, delta_eps) = calculate_nutation(jce);

        // Nutation values should be small (typically < 0.01 degrees)
        assert!(delta_psi.degrees().abs() < 0.1);
        assert!(delta_eps.degrees().abs() < 0.1);
    }

    #[test]
    fn test_mean_ecliptic_obliquity() {
        let jme = JDay::<JulianEphemerisMiliennium>::new(0.0);

        let e0 = calculate_mean_ecliptic_obliquity(jme);

        // At J2000.0, obliquity should be approximately 23.4 degrees
        assert!((e0.value() - 23.4).abs() < 1.0);
    }

    #[test]
    fn test_nutation_type_safety() {
        let jce = JDay::<JulianEphemerisCentury>::new(0.0);
        let (delta_psi, delta_eps) = calculate_nutation(jce);

        // Type system ensures we can't mix these up
        let _: Nutation<NutationLongitude> = delta_psi;
        let _: Nutation<NutationObliquity> = delta_eps;
    }
}
