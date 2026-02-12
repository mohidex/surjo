use chrono::{DateTime, Datelike, Utc};

use super::error::DeltaTError;
use crate::types::julian_day::{JDay, JulianDay, JulianEphemerisDay};

/// Represents Delta T (ΔT) in seconds
/// The difference between Terrestrial Time (TT) and Universal Time (UT)
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DeltaT(f64);

impl DeltaT {
    /// Create a new DeltaT value in seconds
    pub fn new(seconds: f64) -> Self {
        Self(seconds)
    }

    /// Get the value in seconds
    pub fn seconds(&self) -> f64 {
        self.0
    }

    /// Get the value in days (for Julian Day calculations)
    pub fn days(&self) -> f64 {
        self.0 / 86400.0
    }
}

/// Polynomial evaluator for Delta T calculations
struct Polynomial {
    /// Coefficients in ascending order: c₀ + c₁x + c₂x² + ...
    coefficients: &'static [f64],
    /// Center point for the variable transformation
    center: f64,
    /// Scale factor for the variable (default: 1.0)
    scale: f64,
}

impl Polynomial {
    /// Create a new polynomial with coefficients, center point, and optional scale
    const fn new(coefficients: &'static [f64], center: f64, scale: f64) -> Self {
        Self {
            coefficients,
            center,
            scale,
        }
    }

    /// Evaluate the polynomial at the given year value
    /// Uses Horner's method for efficient evaluation
    fn eval(&self, y: f64) -> f64 {
        let u = (y - self.center) / self.scale;

        // Horner's method: p(x) = c₀ + x(c₁ + x(c₂ + x(...)))
        self.coefficients
            .iter()
            .rev()
            .fold(0.0, |acc, &coef| acc * u + coef)
    }
}

/// Parabolic formula: ΔT = -20 + 32u² where u = (y - 1820)/100
const PARABOLIC: Polynomial = Polynomial::new(&[-20.0, 0.0, 32.0], 1820.0, 100.0);

/// Year < -500: Parabolic extrapolation
const BEFORE_500BCE: Polynomial = PARABOLIC;

/// -500 ≤ Year < 500: 6th degree polynomial
/// ΔT = 10583.6 - 1014.41t + 33.78311t² - 5.952053t³ - 0.1798452t⁴ + 0.022174192t⁵ + 0.0090316521t⁶
/// where t = y/100
const ERA_500BCE_500CE: Polynomial = Polynomial::new(
    &[
        10583.6,
        -1014.41,
        33.78311,
        -5.952053,
        -0.1798452,
        0.022174192,
        0.0090316521,
    ],
    0.0,
    100.0,
);

/// 500 ≤ Year < 1600: 6th degree polynomial centered at 1000
/// ΔT = 1574.2 - 556.01u + 71.23472u² + 0.319781u³ - 0.8503463u⁴ - 0.005050998u⁵ + 0.0083572073u⁶
/// where u = (y - 1000)/100
const ERA_500_1600: Polynomial = Polynomial::new(
    &[
        1574.2,
        -556.01,
        71.23472,
        0.319781,
        -0.8503463,
        -0.005050998,
        0.0083572073,
    ],
    1000.0,
    100.0,
);

/// 1600 ≤ Year < 1700: 3rd degree polynomial
/// ΔT = 120 - 0.9808u - 0.01532u² + u³/7129
/// where u = y - 1600
const ERA_1600_1700: Polynomial =
    Polynomial::new(&[120.0, -0.9808, -0.01532, 1.0 / 7129.0], 1600.0, 1.0);

/// 1700 ≤ Year < 1800: 4th degree polynomial
/// ΔT = 8.83 + 0.1603u - 0.0059285u² + 0.00013336u³ - u⁴/1174000
/// where u = y - 1700
const ERA_1700_1800: Polynomial = Polynomial::new(
    &[8.83, 0.1603, -0.0059285, 0.00013336, -1.0 / 1174000.0],
    1700.0,
    1.0,
);

/// 1800 ≤ Year < 1860: 7th degree polynomial
/// ΔT = 13.72 - 0.332447u + 0.0068612u² + 0.0041116u³ - 0.00037436u⁴
///      + 0.0000121272u⁵ - 0.0000001699u⁶ + 0.000000000875u⁷
/// where u = y - 1800
const ERA_1800_1860: Polynomial = Polynomial::new(
    &[
        13.72,
        -0.332447,
        0.0068612,
        0.0041116,
        -0.00037436,
        0.0000121272,
        -0.0000001699,
        0.000000000875,
    ],
    1800.0,
    1.0,
);

/// 1860 ≤ Year < 1900: 5th degree polynomial
/// ΔT = 7.62 + 0.5737u - 0.251754u² + 0.01680668u³ - 0.0004473624u⁴ + u⁵/233174
/// where u = y - 1860
const ERA_1860_1900: Polynomial = Polynomial::new(
    &[
        7.62,
        0.5737,
        -0.251754,
        0.01680668,
        -0.0004473624,
        1.0 / 233174.0,
    ],
    1860.0,
    1.0,
);

/// 1900 ≤ Year < 1920: 4th degree polynomial
/// ΔT = -2.79 + 1.494119u - 0.0598939u² + 0.0061966u³ - 0.000197u⁴
/// where u = y - 1900
const ERA_1900_1920: Polynomial = Polynomial::new(
    &[-2.79, 1.494119, -0.0598939, 0.0061966, -0.000197],
    1900.0,
    1.0,
);

/// 1920 ≤ Year < 1941: 3rd degree polynomial
/// ΔT = 21.20 + 0.84493u - 0.076100u² + 0.0020936u³
/// where u = y - 1920
const ERA_1920_1941: Polynomial =
    Polynomial::new(&[21.20, 0.84493, -0.076100, 0.0020936], 1920.0, 1.0);

/// 1941 ≤ Year < 1961: 3rd degree polynomial centered at 1950
/// ΔT = 29.07 + 0.407u - u²/233 + u³/2547
/// where u = y - 1950
const ERA_1941_1961: Polynomial =
    Polynomial::new(&[29.07, 0.407, -1.0 / 233.0, 1.0 / 2547.0], 1950.0, 1.0);

/// 1961 ≤ Year < 1986: 3rd degree polynomial centered at 1975
/// ΔT = 45.45 + 1.067u - u²/260 - u³/718
/// where u = y - 1975
const ERA_1961_1986: Polynomial =
    Polynomial::new(&[45.45, 1.067, -1.0 / 260.0, -1.0 / 718.0], 1975.0, 1.0);

/// 1986 ≤ Year < 2005: 5th degree polynomial centered at 2000
/// ΔT = 63.86 + 0.3345u - 0.060374u² + 0.0017275u³ + 0.000651814u⁴ + 0.00002373599u⁵
/// where u = y - 2000
const ERA_1986_2005: Polynomial = Polynomial::new(
    &[
        63.86,
        0.3345,
        -0.060374,
        0.0017275,
        0.000651814,
        0.00002373599,
    ],
    2000.0,
    1.0,
);

/// 2005 ≤ Year < 2050: 2nd degree polynomial
/// ΔT = 62.92 + 0.32217u + 0.005589u²
/// where u = y - 2000
const ERA_2005_2050: Polynomial = Polynomial::new(&[62.92, 0.32217, 0.005589], 2000.0, 1.0);

/// Calculate the difference between Terrestrial Dynamical Time (TD)
/// and Universal Time (UT).
///
/// Equations taken from http://eclipse.gsfc.nasa.gov/SEcat5/deltatpoly.html
///
/// # Arguments
///
/// * `year` - The year (can be negative for BCE dates)
/// * `month` - The month (1-12)
///
/// # Returns
///
/// `Result<DeltaT, DeltaTError>` - Delta T in seconds or an error if year is out of range
///
/// # Errors
///
/// Returns `DeltaTError::YearOutOfRange` if year is before -1999 or after 3000.
pub fn calculate_deltat(year: i32, month: u8) -> Result<DeltaT, DeltaTError> {
    // Check valid range
    if !(-1999..=3000).contains(&year) {
        return Err(DeltaTError::YearOutOfRange(year));
    }

    let y = year as f64 + (month as f64 - 0.5) / 12.0;

    let seconds = match year {
        year if year < -500 => BEFORE_500BCE.eval(y),
        year if year < 500 => ERA_500BCE_500CE.eval(y),
        year if year < 1600 => ERA_500_1600.eval(y),
        year if year < 1700 => ERA_1600_1700.eval(y),
        year if year < 1800 => ERA_1700_1800.eval(y),
        year if year < 1860 => ERA_1800_1860.eval(y),
        year if year < 1900 => ERA_1860_1900.eval(y),
        year if year < 1920 => ERA_1900_1920.eval(y),
        year if year < 1941 => ERA_1920_1941.eval(y),
        year if year < 1961 => ERA_1941_1961.eval(y),
        year if year < 1986 => ERA_1961_1986.eval(y),
        year if year < 2005 => ERA_1986_2005.eval(y),
        year if year < 2050 => ERA_2005_2050.eval(y),
        year if year < 2150 => {
            // 2050 ≤ Year < 2150: Parabolic with linear correction
            // ΔT = -20 + 32u² - 0.5628(2150 - y)
            // where u = (y - 1820)/100
            PARABOLIC.eval(y) - 0.5628 * (2150.0 - y)
        }
        _ => {
            // Year ≥ 2150: Simple parabolic extrapolation
            PARABOLIC.eval(y)
        }
    };

    Ok(DeltaT::new(seconds))
}

/// TryFrom implementation for (year, month) tuple
impl TryFrom<(i32, u8)> for DeltaT {
    type Error = DeltaTError;

    fn try_from((year, month): (i32, u8)) -> Result<Self, Self::Error> {
        calculate_deltat(year, month)
    }
}

/// TryFrom implementation for DateTime
impl TryFrom<DateTime<Utc>> for DeltaT {
    type Error = DeltaTError;

    fn try_from(dt: DateTime<Utc>) -> Result<Self, Self::Error> {
        calculate_deltat(dt.year(), dt.month() as u8)
    }
}

/// Convert DeltaT and JulianDay to JulianEphemerisDay
impl From<(JDay<JulianDay>, DeltaT)> for JDay<JulianEphemerisDay> {
    fn from((jd, delta_t): (JDay<JulianDay>, DeltaT)) -> Self {
        JDay::new(jd.value() + delta_t.days())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deltat_newtype() {
        let dt = DeltaT::new(63.8);
        assert_eq!(dt.seconds(), 63.8);
        assert!((dt.days() - 63.8 / 86400.0).abs() < 1e-10);
    }

    #[test]
    fn test_year_2000() {
        let delta_t = calculate_deltat(2000, 1).unwrap();
        assert!((delta_t.seconds() - 63.8).abs() < 1.0);
    }

    #[test]
    fn test_year_1900() {
        let delta_t = calculate_deltat(1900, 1).unwrap();
        assert!((delta_t.seconds() - (-2.72)).abs() < 1.0);
    }

    #[test]
    fn test_year_out_of_range() {
        // Test year too early
        let result = calculate_deltat(-2000, 6);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), DeltaTError::YearOutOfRange(-2000));

        // Test year too late
        let result = calculate_deltat(3001, 6);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), DeltaTError::YearOutOfRange(3001));
    }

    #[test]
    fn test_try_from_tuple() {
        let delta_t: DeltaT = (2000, 1).try_into().unwrap();
        assert!((delta_t.seconds() - 63.8).abs() < 1.0);

        let result: Result<DeltaT, _> = (3001, 1).try_into();
        assert!(result.is_err());
    }

    #[test]
    fn test_try_from_datetime() {
        let dt = DateTime::parse_from_rfc3339("2000-01-01T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);

        let delta_t: DeltaT = dt.try_into().unwrap();
        assert!((delta_t.seconds() - 63.8).abs() < 1.0);
    }

    #[test]
    fn test_jd_to_jde_conversion() {
        let jd = JDay::<JulianDay>::new(2451545.0);
        let delta_t = DeltaT::new(64.0);

        let jde: JDay<JulianEphemerisDay> = (jd, delta_t).into();

        assert!((jde.value() - (2451545.0 + 64.0 / 86400.0)).abs() < 1e-10);
    }

    #[test]
    fn test_full_conversion_chain() {
        let dt = DateTime::parse_from_rfc3339("2000-01-01T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);

        let jd: JDay<JulianDay> = dt.into();
        let delta_t: DeltaT = dt.try_into().unwrap();
        let jde: JDay<JulianEphemerisDay> = (jd, delta_t).into();

        assert!(jde.value() > jd.value());
        assert!((jde.value() - jd.value()) < 0.001);
    }

    #[test]
    fn test_polynomial_evaluation() {
        // Test that polynomial evaluation gives same results
        let delta_t_2000 = calculate_deltat(2000, 1).unwrap();
        let delta_t_1900 = calculate_deltat(1900, 1).unwrap();
        let delta_t_1800 = calculate_deltat(1800, 1).unwrap();

        assert!((delta_t_2000.seconds() - 63.8).abs() < 1.0);
        assert!((delta_t_1900.seconds() - (-2.72)).abs() < 1.0);
        assert!((delta_t_1800.seconds() - 13.7).abs() < 1.0);
    }

    #[test]
    fn test_error_display() {
        let err = DeltaTError::YearOutOfRange(3001);
        let msg = format!("{}", err);
        assert!(msg.contains("3001"));
        assert!(msg.contains("outside valid range"));
    }
}
