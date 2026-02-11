use chrono::{DateTime, Utc};
use std::marker::PhantomData;

/// Marker type for Julian Day Number (JD)
/// Represents days since noon on January 1, 4713 BC (Julian calendar)
#[derive(Debug, Clone, Copy)]
pub struct JulianDay;

/// Marker type for Modified Julian Day (MJD)
/// Represents days since midnight on November 17, 1858 (Gregorian calendar)
/// MJD = JD - 2400000.5
#[derive(Debug, Clone, Copy)]
pub struct ModifiedJulianDay;

/// Marker type for Julian Ephemeris Day (JDE)
/// Julian Day adjusted for Delta T (difference between UT and TT)
#[derive(Debug, Clone, Copy)]
pub struct JulianEphemerisDay;

/// Marker type for Julian Century (JC)
/// Number of Julian centuries since J2000.0 epoch (2451545.0 JD)
#[derive(Debug, Clone, Copy)]
pub struct JulianCentury;

/// Marker type for Julian Ephemeris Century (JCE)
/// Number of Julian centuries based on Ephemeris Day
#[derive(Debug, Clone, Copy)]
pub struct JulianEphemerisCentury;

/// Marker type for Julian Ephemeris Millennium (JME)
/// Number of Julian millennia based on Ephemeris Century
#[derive(Debug, Clone, Copy)]
pub struct JulianEphemerisMiliennium;

/// Type-safe Julian Day representation
/// The generic parameter `T` ensures type safety between different day representations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct JDay<T>(f64, PhantomData<T>);

impl<T> JDay<T> {
    /// Create a new Julian day value
    pub fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    /// Get the underlying numeric value
    pub fn value(&self) -> f64 {
        self.0
    }
}

/// Conversion from Modified Julian Day to Julian Day
impl From<JDay<ModifiedJulianDay>> for JDay<JulianDay> {
    fn from(mjd: JDay<ModifiedJulianDay>) -> Self {
        JDay::new(mjd.value() + 2400000.5)
    }
}

/// Conversion from Julian Day to Modified Julian Day
impl From<JDay<JulianDay>> for JDay<ModifiedJulianDay> {
    fn from(jd: JDay<JulianDay>) -> Self {
        JDay::new(jd.value() - 2400000.5)
    }
}

/// Conversion from Unix timestamp to Julian Day
impl From<f64> for JDay<JulianDay> {
    fn from(unixtime: f64) -> Self {
        let jd = unixtime / 86400.0 + 2440587.5;
        JDay::new(jd)
    }
}

/// Conversion from Julian Day and Delta T (as tuple) to Julian Ephemeris Day
/// Delta T is in seconds
impl From<(JDay<JulianDay>, f64)> for JDay<JulianEphemerisDay> {
    fn from((jd, delta_t): (JDay<JulianDay>, f64)) -> Self {
        let jde = jd.value() + delta_t / 86400.0;
        JDay::new(jde)
    }
}

/// Conversion from Julian Day to Julian Century
/// Referenced to J2000.0 epoch (2451545.0 JD)
impl From<JDay<JulianDay>> for JDay<JulianCentury> {
    fn from(jd: JDay<JulianDay>) -> Self {
        let jc = (jd.value() - 2451545.0) / 36525.0;
        JDay::new(jc)
    }
}

/// Conversion from Julian Ephemeris Day to Julian Ephemeris Century
impl From<JDay<JulianEphemerisDay>> for JDay<JulianEphemerisCentury> {
    fn from(jde: JDay<JulianEphemerisDay>) -> Self {
        let jce = (jde.value() - 2451545.0) / 36525.0;
        JDay::new(jce)
    }
}

/// Conversion from Julian Ephemeris Century to Julian Ephemeris Millennium
impl From<JDay<JulianEphemerisCentury>> for JDay<JulianEphemerisMiliennium> {
    fn from(jce: JDay<JulianEphemerisCentury>) -> Self {
        let jme = jce.value() / 10.0;
        JDay::new(jme)
    }
}

impl From<DateTime<Utc>> for JDay<JulianDay> {
    fn from(dt: DateTime<Utc>) -> Self {
        let unixtime = dt.timestamp() as f64 + dt.timestamp_subsec_nanos() as f64 / 1e9;
        JDay::<JulianDay>::from(unixtime)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unix_timestamp_conversion() {
        // Unix epoch (1970-01-01 00:00:00) = JD 2440587.5
        let jd: JDay<JulianDay> = 0.0.into();
        assert!((jd.value() - 2440587.5).abs() < 1e-10);
    }

    #[test]
    fn test_datetime_to_jd() {
        // 2000-01-01 12:00:00 UTC = JD 2451545.0 (J2000.0 epoch)
        let dt = DateTime::parse_from_rfc3339("2000-01-01T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);
        let jd: JDay<JulianDay> = dt.into();
        assert!((jd.value() - 2451545.0).abs() < 1e-6);
    }

    #[test]
    fn test_jd_mjd_conversion() {
        let jd = JDay::<JulianDay>::new(2451545.0);
        let mjd: JDay<ModifiedJulianDay> = jd.into();
        assert!((mjd.value() - 51544.5).abs() < 1e-10);

        let jd_back: JDay<JulianDay> = mjd.into();
        assert!((jd_back.value() - 2451545.0).abs() < 1e-10);
    }

    #[test]
    fn test_ephemeris_day_conversion() {
        let jd = JDay::<JulianDay>::new(2451545.0);
        let delta_t = 64.0; // seconds
        let jde: JDay<JulianEphemerisDay> = (jd, delta_t).into();
        assert!((jde.value() - (2451545.0 + 64.0 / 86400.0)).abs() < 1e-10);
    }

    #[test]
    fn test_julian_century_conversion() {
        // At J2000.0 epoch, century should be 0
        let jd = JDay::<JulianDay>::new(2451545.0);
        let jc: JDay<JulianCentury> = jd.into();
        assert!(jc.value().abs() < 1e-10);

        // One century later
        let jd = JDay::<JulianDay>::new(2451545.0 + 36525.0);
        let jc: JDay<JulianCentury> = jd.into();
        assert!((jc.value() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_conversion_chain() {
        let dt = DateTime::parse_from_rfc3339("2000-01-01T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);

        let jd: JDay<JulianDay> = dt.into();
        let jde: JDay<JulianEphemerisDay> = (jd, 64.0).into();
        let jce: JDay<JulianEphemerisCentury> = jde.into();
        let jme: JDay<JulianEphemerisMiliennium> = jce.into();

        assert!((jd.value() - 2451545.0).abs() < 1e-6);
        assert!(jce.value().abs() < 1e-6);
        assert!(jme.value().abs() < 1e-6);
    }
}
