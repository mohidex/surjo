use std::f64::consts::PI;
use std::marker::PhantomData;

use crate::types::angle::{Angle, Degrees};
use crate::types::julian_day::{JDay, JulianEphemerisMiliennium};

/// Marker type for Heliocentric Longitude
#[derive(Debug, Clone, Copy)]
pub struct HeliocentricLongitude;

/// Marker type for Heliocentric Latitude
#[derive(Debug, Clone, Copy)]
pub struct HeliocentricLatitude;

/// Marker type for Heliocentric Radius Vector (in AU)
#[derive(Debug, Clone, Copy)]
pub struct HeliocentricRadius;

/// Type-safe heliocentric coordinate
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Heliocentric<T>(f64, PhantomData<T>);

impl<T> Heliocentric<T> {
    /// Create a new heliocentric value
    pub fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    /// Get the underlying value
    pub fn value(&self) -> f64 {
        self.0
    }
}

/// Coefficient table for periodic terms
/// Each row: [A, B, C] where term = A * cos(B + C*JME)
struct CoefficientTable {
    coefficients: &'static [[f64; 3]],
}

impl CoefficientTable {
    const fn new(coefficients: &'static [[f64; 3]]) -> Self {
        Self { coefficients }
    }

    /// Evaluate Σ(A * cos(B + C*x))
    fn eval(&self, x: f64) -> f64 {
        self.coefficients
            .iter()
            .map(|&[a, b, c]| a * (b + c * x).cos())
            .sum()
    }
}

const L0: CoefficientTable = CoefficientTable::new(&[
    [175347046.0, 0.0, 0.0],
    [3341656.0, 4.6692568, 6283.07585],
    [34894.0, 4.6261, 12566.1517],
    [3497.0, 2.7441, 5753.3849],
    [3418.0, 2.8289, 3.5231],
    [3136.0, 3.6277, 77713.7715],
    [2676.0, 4.4181, 7860.4194],
    [2343.0, 6.1352, 3930.2097],
    [1324.0, 0.7425, 11506.7698],
    [1273.0, 2.0371, 529.691],
    [1199.0, 1.1096, 1577.3435],
    [990.0, 5.233, 5884.927],
    [902.0, 2.045, 26.298],
    [857.0, 3.508, 398.149],
    [780.0, 1.179, 5223.694],
    [753.0, 2.533, 5507.553],
    [505.0, 4.583, 18849.228],
    [492.0, 4.205, 775.523],
    [357.0, 2.92, 0.067],
    [317.0, 5.849, 11790.629],
    [284.0, 1.899, 796.298],
    [271.0, 0.315, 10977.079],
    [243.0, 0.345, 5486.778],
    [206.0, 4.806, 2544.314],
    [205.0, 1.869, 5573.143],
    [202.0, 2.458, 6069.777],
    [156.0, 0.833, 213.299],
    [132.0, 3.411, 2942.463],
    [126.0, 1.083, 20.775],
    [115.0, 0.645, 0.98],
    [103.0, 0.636, 4694.003],
    [102.0, 0.976, 15720.839],
    [102.0, 4.267, 7.114],
    [99.0, 6.21, 2146.17],
    [98.0, 0.68, 155.42],
    [86.0, 5.98, 161000.69],
    [85.0, 1.3, 6275.96],
    [85.0, 3.67, 71430.7],
    [80.0, 1.81, 17260.15],
    [79.0, 3.04, 12036.46],
    [75.0, 1.76, 5088.63],
    [74.0, 3.5, 3154.69],
    [74.0, 4.68, 801.82],
    [70.0, 0.83, 9437.76],
    [62.0, 3.98, 8827.39],
    [61.0, 1.82, 7084.9],
    [57.0, 2.78, 6286.6],
    [56.0, 4.39, 14143.5],
    [56.0, 3.47, 6279.55],
    [52.0, 0.19, 12139.55],
    [52.0, 1.33, 1748.02],
    [51.0, 0.28, 5856.48],
    [49.0, 0.49, 1194.45],
    [41.0, 5.37, 8429.24],
    [41.0, 2.4, 19651.05],
    [39.0, 6.17, 10447.39],
    [37.0, 6.04, 10213.29],
    [37.0, 2.57, 1059.38],
    [36.0, 1.71, 2352.87],
    [36.0, 1.78, 6812.77],
    [33.0, 0.59, 17789.85],
    [30.0, 0.44, 83996.85],
    [30.0, 2.74, 1349.87],
    [25.0, 3.16, 4690.48],
]);

const L1: CoefficientTable = CoefficientTable::new(&[
    [628331966747.0, 0.0, 0.0],
    [206059.0, 2.678235, 6283.07585],
    [4303.0, 2.6351, 12566.1517],
    [425.0, 1.59, 3.523],
    [119.0, 5.796, 26.298],
    [109.0, 2.966, 1577.344],
    [93.0, 2.59, 18849.23],
    [72.0, 1.14, 529.69],
    [68.0, 1.87, 398.15],
    [67.0, 4.41, 5507.55],
    [59.0, 2.89, 5223.69],
    [56.0, 2.17, 155.42],
    [45.0, 0.4, 796.3],
    [36.0, 0.47, 775.52],
    [29.0, 2.65, 7.11],
    [21.0, 5.34, 0.98],
    [19.0, 1.85, 5486.78],
    [19.0, 4.97, 213.3],
    [17.0, 2.99, 6275.96],
    [16.0, 0.03, 2544.31],
    [16.0, 1.43, 2146.17],
    [15.0, 1.21, 10977.08],
    [12.0, 2.83, 1748.02],
    [12.0, 3.26, 5088.63],
    [12.0, 5.27, 1194.45],
    [12.0, 2.08, 4694.0],
    [11.0, 0.77, 553.57],
    [10.0, 1.3, 6286.6],
    [10.0, 4.24, 1349.87],
    [9.0, 2.7, 242.73],
    [9.0, 5.64, 951.72],
    [8.0, 5.3, 2352.87],
    [6.0, 2.65, 9437.76],
    [6.0, 4.67, 4690.48],
]);

const L2: CoefficientTable = CoefficientTable::new(&[
    [52919.0, 0.0, 0.0],
    [8720.0, 1.0721, 6283.0758],
    [309.0, 0.867, 12566.152],
    [27.0, 0.05, 3.52],
    [16.0, 5.19, 26.3],
    [16.0, 3.68, 155.42],
    [10.0, 0.76, 18849.23],
    [9.0, 2.06, 77713.77],
    [7.0, 0.83, 775.52],
    [5.0, 4.66, 1577.34],
    [4.0, 1.03, 7.11],
    [4.0, 3.44, 5573.14],
    [3.0, 5.14, 796.3],
    [3.0, 6.05, 5507.55],
    [3.0, 1.19, 242.73],
    [3.0, 6.12, 529.69],
    [3.0, 0.31, 398.15],
    [3.0, 2.28, 553.57],
    [2.0, 4.38, 5223.69],
    [2.0, 3.75, 0.98],
]);

const L3: CoefficientTable = CoefficientTable::new(&[
    [289.0, 5.844, 6283.076],
    [35.0, 0.0, 0.0],
    [17.0, 5.49, 12566.15],
    [3.0, 5.2, 155.42],
    [1.0, 4.72, 3.52],
    [1.0, 5.3, 18849.23],
    [1.0, 5.97, 242.73],
]);

const L4: CoefficientTable = CoefficientTable::new(&[
    [114.0, PI, 0.0],
    [8.0, 4.13, 6283.08],
    [1.0, 3.84, 12566.15],
]);

const L5: CoefficientTable = CoefficientTable::new(&[[1.0, PI, 0.0]]);

const B0: CoefficientTable = CoefficientTable::new(&[
    [280.0, 3.199, 84334.662],
    [102.0, 5.422, 5507.553],
    [80.0, 3.88, 5223.69],
    [44.0, 3.7, 2352.87],
    [32.0, 4.0, 1577.34],
]);

const B1: CoefficientTable = CoefficientTable::new(&[[9.0, 3.9, 5507.55], [6.0, 1.73, 5223.69]]);

const R0: CoefficientTable = CoefficientTable::new(&[
    [100013989.0, 0.0, 0.0],
    [1670700.0, 3.0984635, 6283.07585],
    [13956.0, 3.05525, 12566.1517],
    [3084.0, 5.1985, 77713.7715],
    [1628.0, 1.1739, 5753.3849],
    [1576.0, 2.8469, 7860.4194],
    [925.0, 5.453, 11506.77],
    [542.0, 4.564, 3930.21],
    [472.0, 3.661, 5884.927],
    [346.0, 0.964, 5507.553],
    [329.0, 5.9, 5223.694],
    [307.0, 0.299, 5573.143],
    [243.0, 4.273, 11790.629],
    [212.0, 5.847, 1577.344],
    [186.0, 5.022, 10977.079],
    [175.0, 3.012, 18849.228],
    [110.0, 5.055, 5486.778],
    [98.0, 0.89, 6069.78],
    [86.0, 5.69, 15720.84],
    [86.0, 1.27, 161000.69],
    [65.0, 0.27, 17260.15],
    [63.0, 0.92, 529.69],
    [57.0, 2.01, 83996.85],
    [56.0, 5.24, 71430.7],
    [49.0, 3.25, 2544.31],
    [47.0, 2.58, 775.52],
    [45.0, 5.54, 9437.76],
    [43.0, 6.01, 6275.96],
    [39.0, 5.36, 4694.0],
    [38.0, 2.39, 8827.39],
    [37.0, 0.83, 19651.05],
    [37.0, 4.9, 12139.55],
    [36.0, 1.67, 12036.46],
    [35.0, 1.84, 2942.46],
    [33.0, 0.24, 7084.9],
    [32.0, 0.18, 5088.63],
    [32.0, 1.78, 398.15],
    [28.0, 1.21, 6286.6],
    [28.0, 1.9, 6279.55],
    [26.0, 4.59, 10447.39],
]);

const R1: CoefficientTable = CoefficientTable::new(&[
    [103019.0, 1.10749, 6283.07585],
    [1721.0, 1.0644, 12566.1517],
    [702.0, PI, 0.0],
    [32.0, 1.02, 18849.23],
    [31.0, 2.84, 5507.55],
    [25.0, 1.32, 5223.69],
    [18.0, 1.42, 1577.34],
    [10.0, 5.91, 10977.08],
    [9.0, 1.42, 6275.96],
    [9.0, 0.27, 5486.78],
]);

const R2: CoefficientTable = CoefficientTable::new(&[
    [4359.0, 5.7846, 6283.0758],
    [124.0, 5.579, 12566.152],
    [12.0, PI, 0.0],
    [9.0, 3.63, 77713.77],
    [6.0, 1.87, 5573.14],
    [3.0, 5.47, 18849.23],
]);

const R3: CoefficientTable =
    CoefficientTable::new(&[[145.0, 4.273, 6283.076], [7.0, 3.92, 12566.15]]);

const R4: CoefficientTable = CoefficientTable::new(&[[4.0, 2.56, 6283.08]]);

/// Calculate Earth's heliocentric longitude in degrees
///
/// L = (L0 + L1*JME + L2*JME² + L3*JME³ + L4*JME⁴ + L5*JME⁵) / 10⁸
///
/// Result is normalized to [0, 360) degrees
pub fn calculate_heliocentric_longitude(jme: JDay<JulianEphemerisMiliennium>) -> Angle<Degrees> {
    let jme_val = jme.value();

    let l0 = L0.eval(jme_val);
    let l1 = L1.eval(jme_val);
    let l2 = L2.eval(jme_val);
    let l3 = L3.eval(jme_val);
    let l4 = L4.eval(jme_val);
    let l5 = L5.eval(jme_val);

    // Polynomial: L = L0 + L1*JME + L2*JME² + L3*JME³ + L4*JME⁴ + L5*JME⁵
    let l_rad = (l0
        + l1 * jme_val
        + l2 * jme_val.powi(2)
        + l3 * jme_val.powi(3)
        + l4 * jme_val.powi(4)
        + l5 * jme_val.powi(5))
        / 1e8;

    let l_deg = l_rad.to_degrees();

    // Normalize to [0, 360)
    Angle::new(l_deg.rem_euclid(360.0))
}

/// Calculate Earth's heliocentric latitude in degrees
///
/// B = (B0 + B1*JME) / 10⁸
pub fn calculate_heliocentric_latitude(jme: JDay<JulianEphemerisMiliennium>) -> Angle<Degrees> {
    let jme_val = jme.value();

    let b0 = B0.eval(jme_val);
    let b1 = B1.eval(jme_val);

    // Polynomial: B = B0 + B1*JME
    let b_rad = (b0 + b1 * jme_val) / 1e8;

    let b_deg = b_rad.to_degrees();

    Angle::new(b_deg)
}

/// Calculate Earth's heliocentric radius vector in Astronomical Units (AU)
///
/// R = (R0 + R1*JME + R2*JME² + R3*JME³ + R4*JME⁴) / 10⁸
pub fn calculate_heliocentric_radius(
    jme: JDay<JulianEphemerisMiliennium>,
) -> Heliocentric<HeliocentricRadius> {
    let jme_val = jme.value();

    let r0 = R0.eval(jme_val);
    let r1 = R1.eval(jme_val);
    let r2 = R2.eval(jme_val);
    let r3 = R3.eval(jme_val);
    let r4 = R4.eval(jme_val);

    // Polynomial: R = R0 + R1*JME + R2*JME² + R3*JME³ + R4*JME⁴
    let r_au =
        (r0 + r1 * jme_val + r2 * jme_val.powi(2) + r3 * jme_val.powi(3) + r4 * jme_val.powi(4))
            / 1e8;

    Heliocentric::new(r_au)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::julian_day::{JulianDay, JulianEphemerisCentury, JulianEphemerisDay};

    #[test]
    fn test_coefficient_table_evaluation() {
        // Test simple coefficient table
        let table = CoefficientTable::new(&[
            [1.0, 0.0, 0.0], // 1 * cos(0) = 1
            [2.0, 0.0, 1.0], // 2 * cos(x)
        ]);

        let result = table.eval(0.0);
        assert!((result - 3.0).abs() < 1e-10); // 1 + 2*cos(0) = 3
    }

    #[test]
    fn test_heliocentric_longitude() {
        // Test at J2000.0 epoch (JME = 0)
        let jme = JDay::<JulianEphemerisMiliennium>::new(0.0);
        let longitude = calculate_heliocentric_longitude(jme);

        // Should be a valid angle in [0, 360)
        assert!(longitude.value() >= 0.0);
        assert!(longitude.value() < 360.0);
    }

    #[test]
    fn test_heliocentric_latitude() {
        // Test at J2000.0 epoch (JME = 0)
        let jme = JDay::<JulianEphemerisMiliennium>::new(0.0);
        let latitude = calculate_heliocentric_latitude(jme);

        // Latitude should be small (Earth's orbit is nearly in ecliptic plane)
        assert!(latitude.value().abs() < 1.0);
    }

    #[test]
    fn test_heliocentric_radius() {
        // Test at J2000.0 epoch (JME = 0)
        let jme = JDay::<JulianEphemerisMiliennium>::new(0.0);
        let radius = calculate_heliocentric_radius(jme);

        // Should be approximately 1 AU
        assert!((radius.value() - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_full_calculation_chain() {
        // Start from JulianDay, convert through the chain
        let jd = JDay::<JulianDay>::new(2451545.0); // J2000.0
        let jde: JDay<JulianEphemerisDay> = (jd, 64.0).into(); // Approximate delta T
        let jce: JDay<JulianEphemerisCentury> = jde.into();
        let jme: JDay<JulianEphemerisMiliennium> = jce.into();

        let longitude = calculate_heliocentric_longitude(jme);
        let latitude = calculate_heliocentric_latitude(jme);
        let radius = calculate_heliocentric_radius(jme);

        // All should be valid values
        assert!(longitude.value() >= 0.0 && longitude.value() < 360.0);
        assert!(latitude.value().abs() < 90.0);
        assert!(radius.value() > 0.0);
    }
}
