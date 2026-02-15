use chrono::{DateTime, Utc};

use super::deltat::DeltaT;
use super::error::DeltaTError;
use super::geocentric::{
    calculate_geocentric_latitude, calculate_geocentric_longitude,
    calculate_mean_ecliptic_obliquity, calculate_nutation,
};
use super::heliocentric::{
    calculate_heliocentric_latitude, calculate_heliocentric_longitude,
    calculate_heliocentric_radius,
};
use crate::types::angle::{Angle, Degrees, Radians};
use crate::types::julian_day::{
    JDay, JulianCentury, JulianDay, JulianEphemerisCentury, JulianEphemerisDay,
    JulianEphemerisMiliennium,
};

/// Solar position calculation result
#[derive(Debug, Clone, Copy)]
pub struct SolarPosition {
    zenith: Angle<Degrees>,
    zenith_no_refraction: Angle<Degrees>,
    elevation: Angle<Degrees>,
    elevation_no_refraction: Angle<Degrees>,
    azimuth: Angle<Degrees>,
    equation_of_time: f64, // minutes
}

impl SolarPosition {
    /// Create a new SolarPosition
    fn new(
        zenith: Angle<Degrees>,
        zenith_no_refraction: Angle<Degrees>,
        elevation: Angle<Degrees>,
        elevation_no_refraction: Angle<Degrees>,
        azimuth: Angle<Degrees>,
        equation_of_time: f64,
    ) -> Self {
        Self {
            zenith,
            zenith_no_refraction,
            elevation,
            elevation_no_refraction,
            azimuth,
            equation_of_time,
        }
    }

    /// Get topocentric zenith angle (with atmospheric refraction)
    pub fn zenith(&self) -> Angle<Degrees> {
        self.zenith
    }

    /// Get topocentric zenith angle (without atmospheric refraction)
    pub fn zenith_no_refraction(&self) -> Angle<Degrees> {
        self.zenith_no_refraction
    }

    /// Get topocentric elevation angle (with atmospheric refraction)
    pub fn elevation(&self) -> Angle<Degrees> {
        self.elevation
    }

    /// Get topocentric elevation angle (without atmospheric refraction)
    pub fn elevation_no_refraction(&self) -> Angle<Degrees> {
        self.elevation_no_refraction
    }

    /// Get topocentric azimuth angle (measured eastward from north)
    pub fn azimuth(&self) -> Angle<Degrees> {
        self.azimuth
    }

    /// Get equation of time in minutes
    pub fn equation_of_time(&self) -> f64 {
        self.equation_of_time
    }

    /// Check if the sun is above the horizon
    pub fn is_daylight(&self) -> bool {
        self.elevation.value() > 0.0
    }

    /// Check if it's civil twilight (sun between 0° and -6°)
    pub fn is_civil_twilight(&self) -> bool {
        let elev = self.elevation_no_refraction.value();
        elev > -6.0 && elev <= 0.0
    }

    /// Check if it's nautical twilight (sun between -6° and -12°)
    pub fn is_nautical_twilight(&self) -> bool {
        let elev = self.elevation_no_refraction.value();
        elev > -12.0 && elev <= -6.0
    }

    /// Check if it's astronomical twilight (sun between -12° and -18°)
    pub fn is_astronomical_twilight(&self) -> bool {
        let elev = self.elevation_no_refraction.value();
        elev > -18.0 && elev <= -12.0
    }

    /// Check if it's night (sun below -18°)
    pub fn is_night(&self) -> bool {
        self.elevation_no_refraction.value() <= -18.0
    }

    /// Get the atmospheric refraction correction applied
    pub fn refraction_correction(&self) -> Angle<Degrees> {
        Angle::new(self.elevation.value() - self.elevation_no_refraction.value())
    }

    /// Get the air mass (approximation using zenith angle)
    /// Uses Kasten and Young formula for air mass
    pub fn air_mass(&self) -> Option<f64> {
        if !self.is_daylight() {
            return None;
        }

        let zenith_rad: Angle<Radians> = self.zenith.into();
        let cos_zenith = zenith_rad.value().cos();

        // Kasten and Young formula
        let air_mass =
            1.0 / (cos_zenith + 0.50572 * (96.07995 - self.zenith.value()).powf(-1.6364));

        Some(air_mass)
    }
}

/// Observer location on Earth
#[derive(Debug, Clone, Copy)]
pub struct Observer {
    latitude: Angle<Degrees>,
    longitude: Angle<Degrees>,
    elevation: f64,
}

/// Error type for Observer construction
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ObserverError {
    InvalidLatitude(f64),
    InvalidLongitude(f64),
    InvalidElevation(f64),
}

impl std::fmt::Display for ObserverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObserverError::InvalidLatitude(lat) => {
                write!(f, "Latitude {} is outside valid range [-90, 90]", lat)
            }
            ObserverError::InvalidLongitude(lon) => {
                write!(f, "Longitude {} is outside valid range [-180, 180]", lon)
            }
            ObserverError::InvalidElevation(elev) => {
                write!(f, "Elevation {} is invalid (must be >= -500)", elev)
            }
        }
    }
}

impl std::error::Error for ObserverError {}

impl Observer {
    /// Create a new Observer with validation
    ///
    /// # Arguments
    ///
    /// * `latitude` - Latitude in degrees (range: -90 to 90, positive north)
    /// * `longitude` - Longitude in degrees (range: -180 to 180, positive east)
    /// * `elevation` - Elevation above sea level in meters (must be >= -500)
    ///
    /// # Errors
    ///
    /// Returns error if any parameter is outside valid range
    pub fn new(latitude: f64, longitude: f64, elevation: f64) -> Result<Self, ObserverError> {
        if !(-90.0..=90.0).contains(&latitude) {
            return Err(ObserverError::InvalidLatitude(latitude));
        }
        if !(-180.0..=180.0).contains(&longitude) {
            return Err(ObserverError::InvalidLongitude(longitude));
        }
        if elevation < -500.0 {
            return Err(ObserverError::InvalidElevation(elevation));
        }

        Ok(Self {
            latitude: Angle::new(latitude),
            longitude: Angle::new(longitude),
            elevation,
        })
    }

    /// Get latitude
    pub fn latitude(&self) -> Angle<Degrees> {
        self.latitude
    }

    /// Get longitude
    pub fn longitude(&self) -> Angle<Degrees> {
        self.longitude
    }

    /// Get elevation in meters above sea level
    pub fn elevation(&self) -> f64 {
        self.elevation
    }

    /// Calculate the u term for parallax correction
    /// u = atan(0.99664719 tan(φ))
    fn u_term(&self) -> Angle<Radians> {
        let lat_rad: Angle<Radians> = self.latitude.into();
        let result = (0.99664719 * lat_rad.value().tan()).atan();
        Angle::new(result)
    }

    /// Calculate the x term for parallax correction
    /// x = cos(u) + (elevation/6378140)cos(φ)
    fn x_term(&self) -> f64 {
        let u = self.u_term();
        let lat_rad: Angle<Radians> = self.latitude.into();
        u.value().cos() + (self.elevation / 6378140.0) * lat_rad.value().cos()
    }

    /// Calculate the y term for parallax correction
    /// y = 0.99664719 sin(u) + (elevation/6378140)sin(φ)
    fn y_term(&self) -> f64 {
        let u = self.u_term();
        let lat_rad: Angle<Radians> = self.latitude.into();
        0.99664719 * u.value().sin() + (self.elevation / 6378140.0) * lat_rad.value().sin()
    }
}

/// Atmospheric conditions for refraction correction
#[derive(Debug, Clone, Copy)]
pub struct Atmosphere {
    pressure: f64,
    temperature: f64,
}

/// Error type for Atmosphere construction
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AtmosphereError {
    InvalidPressure(f64),
    InvalidTemperature(f64),
}

impl std::fmt::Display for AtmosphereError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AtmosphereError::InvalidPressure(p) => {
                write!(f, "Pressure {} mb is invalid (must be > 0 and < 2000)", p)
            }
            AtmosphereError::InvalidTemperature(t) => {
                write!(
                    f,
                    "Temperature {} °C is invalid (must be between -273 and 100)",
                    t
                )
            }
        }
    }
}

impl std::error::Error for AtmosphereError {}

impl Atmosphere {
    /// Standard atmospheric conditions at sea level
    pub const STANDARD: Self = Self {
        pressure: 1013.25,
        temperature: 15.0,
    };

    /// Create new atmospheric conditions with validation
    ///
    /// # Arguments
    ///
    /// * `pressure` - Atmospheric pressure in millibars (must be > 0 and < 2000)
    /// * `temperature` - Temperature in degrees Celsius (must be between -273 and 100)
    ///
    /// # Errors
    ///
    /// Returns error if parameters are outside valid range
    pub fn new(pressure: f64, temperature: f64) -> Result<Self, AtmosphereError> {
        if pressure <= 0.0 || pressure >= 2000.0 {
            return Err(AtmosphereError::InvalidPressure(pressure));
        }
        if !(-273.0..=100.0).contains(&temperature) {
            return Err(AtmosphereError::InvalidTemperature(temperature));
        }

        Ok(Self {
            pressure,
            temperature,
        })
    }

    /// Get atmospheric pressure in millibars
    pub fn pressure(&self) -> f64 {
        self.pressure
    }

    /// Get temperature in degrees Celsius
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Calculate atmospheric refraction correction
    /// Δe = (P/1010)(283/(273+T))(1.02/(60 tan(e₀ + 10.3/(e₀ + 5.11))))
    pub fn refraction_correction(&self, elevation_no_refraction: Angle<Degrees>) -> Angle<Degrees> {
        let e0 = elevation_no_refraction.value();

        if e0 < -1.0 {
            // Sun is well below horizon
            return Angle::new(0.0);
        }

        let _e0_rad = e0.to_radians();
        let denominator = (e0 + 10.3 / (e0 + 5.11)).to_radians().tan();

        let correction = (self.pressure / 1010.0)
            * (283.0 / (273.0 + self.temperature))
            * (1.02 / (60.0 * denominator));

        Angle::new(correction)
    }
}

impl Default for Atmosphere {
    fn default() -> Self {
        Self::STANDARD
    }
}

// ============================================================================
// Helper Calculations
// ============================================================================

/// Calculate true ecliptic obliquity
/// ε = ε₀ + Δε
fn true_ecliptic_obliquity(epsilon0: Angle<Degrees>, delta_epsilon: f64) -> Angle<Degrees> {
    Angle::new(epsilon0.value() + delta_epsilon)
}

/// Calculate aberration correction in degrees
/// Δτ = -20.4898 / (3600 * R)
fn aberration_correction(radius: f64) -> f64 {
    -20.4898 / (3600.0 * radius)
}

/// Calculate apparent sun longitude
/// λ = Θ + Δψ + Δτ
fn apparent_sun_longitude(theta: Angle<Degrees>, delta_psi: f64, delta_tau: f64) -> Angle<Degrees> {
    Angle::new(theta.value() + delta_psi + delta_tau)
}

/// Calculate mean sidereal time at Greenwich in degrees
/// ν₀ = 280.46061837 + 360.98564736629(JD - 2451545) + 0.000387933JC² - JC³/38710000
fn mean_sidereal_time(jd: JDay<JulianDay>, jc: JDay<JulianCentury>) -> f64 {
    let jd_val = jd.value();
    let jc_val = jc.value();

    let v0 = 280.46061837 + 360.98564736629 * (jd_val - 2451545.0) + 0.000387933 * jc_val.powi(2)
        - jc_val.powi(3) / 38710000.0;

    v0.rem_euclid(360.0)
}

/// Calculate apparent sidereal time at Greenwich in degrees
/// ν = ν₀ + Δψ cos(ε)
fn apparent_sidereal_time(v0: f64, delta_psi: f64, epsilon: Angle<Degrees>) -> f64 {
    let epsilon_rad: Angle<Radians> = epsilon.into();
    v0 + delta_psi * epsilon_rad.value().cos()
}

/// Calculate geocentric sun right ascension in degrees
/// α = atan2(sin(λ)cos(ε) - tan(β)sin(ε), cos(λ))
fn geocentric_sun_right_ascension(
    lambda: Angle<Degrees>,
    epsilon: Angle<Degrees>,
    beta: Angle<Degrees>,
) -> f64 {
    let lambda_rad: Angle<Radians> = lambda.into();
    let epsilon_rad: Angle<Radians> = epsilon.into();
    let beta_rad: Angle<Radians> = beta.into();

    let numerator = lambda_rad.value().sin() * epsilon_rad.value().cos()
        - beta_rad.value().tan() * epsilon_rad.value().sin();
    let denominator = lambda_rad.value().cos();

    numerator.atan2(denominator).to_degrees().rem_euclid(360.0)
}

/// Calculate geocentric sun declination in degrees
/// δ = asin(sin(β)cos(ε) + cos(β)sin(ε)sin(λ))
fn geocentric_sun_declination(
    lambda: Angle<Degrees>,
    epsilon: Angle<Degrees>,
    beta: Angle<Degrees>,
) -> f64 {
    let lambda_rad: Angle<Radians> = lambda.into();
    let epsilon_rad: Angle<Radians> = epsilon.into();
    let beta_rad: Angle<Radians> = beta.into();

    let arg = beta_rad.value().sin() * epsilon_rad.value().cos()
        + beta_rad.value().cos() * epsilon_rad.value().sin() * lambda_rad.value().sin();

    arg.asin().to_degrees()
}

/// Calculate sun mean longitude in degrees
/// M = 280.4664567 + 360007.6982779JME + 0.03032028JME² + JME³/49931 - JME⁴/15300 - JME⁵/2000000
fn sun_mean_longitude(jme: JDay<JulianEphemerisMiliennium>) -> f64 {
    let jme_val = jme.value();

    let m = 280.4664567
        + 360007.6982779 * jme_val
        + 0.03032028 * jme_val.powi(2)
        + jme_val.powi(3) / 49931.0
        - jme_val.powi(4) / 15300.0
        - jme_val.powi(5) / 2000000.0;

    m.rem_euclid(360.0)
}

/// Calculate equation of time in minutes
/// E = 4(M - 0.0057183 - α + Δψ cos(ε))
fn equation_of_time(m: f64, alpha: f64, delta_psi: f64, epsilon: Angle<Degrees>) -> f64 {
    let epsilon_rad: Angle<Radians> = epsilon.into();
    4.0 * (m - 0.0057183 - alpha + delta_psi * epsilon_rad.value().cos())
}

/// Calculate local hour angle in degrees
/// H = ν + longitude - α
fn local_hour_angle(v: f64, longitude: f64, alpha: f64) -> f64 {
    (v + longitude - alpha).rem_euclid(360.0)
}

/// Calculate equatorial horizontal parallax in degrees
/// ξ = 8.794 / (3600 * R)
fn equatorial_horizontal_parallax(radius: f64) -> f64 {
    8.794 / (3600.0 * radius)
}

/// Calculate u term for parallax correction
/// u = atan(0.99664719 tan(φ))
fn uterm(latitude: f64) -> f64 {
    let lat_rad = latitude.to_radians();
    (0.99664719 * lat_rad.tan()).atan()
}

/// Calculate x term for parallax correction
/// x = cos(u) + (elevation/6378140)cos(φ)
fn xterm(u: f64, latitude: f64, elevation: f64) -> f64 {
    let lat_rad = latitude.to_radians();
    u.cos() + (elevation / 6378140.0) * lat_rad.cos()
}

/// Calculate y term for parallax correction
/// y = 0.99664719 sin(u) + (elevation/6378140)sin(φ)
fn yterm(u: f64, latitude: f64, elevation: f64) -> f64 {
    let lat_rad = latitude.to_radians();
    0.99664719 * u.sin() + (elevation / 6378140.0) * lat_rad.sin()
}

/// Calculate parallax in sun right ascension in degrees
/// Δα = atan2(-x sin(ξ) sin(H), cos(δ) - x sin(ξ) cos(H))
fn parallax_sun_right_ascension(x: f64, xi: f64, h: f64, delta: f64) -> f64 {
    let xi_rad = xi.to_radians();
    let h_rad = h.to_radians();
    let delta_rad = delta.to_radians();

    let numerator = -x * xi_rad.sin() * h_rad.sin();
    let denominator = delta_rad.cos() - x * xi_rad.sin() * h_rad.cos();

    numerator.atan2(denominator).to_degrees()
}

/// Calculate topocentric sun declination in degrees
/// δ' = atan2((sin(δ) - y sin(ξ)) cos(Δα), cos(δ) - x sin(ξ) cos(H))
fn topocentric_sun_declination(
    delta: f64,
    x: f64,
    y: f64,
    xi: f64,
    delta_alpha: f64,
    h: f64,
) -> f64 {
    let delta_rad = delta.to_radians();
    let xi_rad = xi.to_radians();
    let delta_alpha_rad = delta_alpha.to_radians();
    let h_rad = h.to_radians();

    let numerator = (delta_rad.sin() - y * xi_rad.sin()) * delta_alpha_rad.cos();
    let denominator = delta_rad.cos() - x * xi_rad.sin() * h_rad.cos();

    numerator.atan2(denominator).to_degrees()
}

/// Calculate topocentric local hour angle in degrees
/// H' = H - Δα
fn topocentric_local_hour_angle(h: f64, delta_alpha: f64) -> f64 {
    h - delta_alpha
}

/// Calculate topocentric elevation angle without atmospheric refraction in degrees
/// e₀ = asin(sin(φ)sin(δ') + cos(φ)cos(δ')cos(H'))
fn topocentric_elevation_angle_without_atmosphere(
    latitude: f64,
    delta_prime: f64,
    h_prime: f64,
) -> f64 {
    let lat_rad = latitude.to_radians();
    let delta_prime_rad = delta_prime.to_radians();
    let h_prime_rad = h_prime.to_radians();

    let arg = lat_rad.sin() * delta_prime_rad.sin()
        + lat_rad.cos() * delta_prime_rad.cos() * h_prime_rad.cos();

    arg.asin().to_degrees()
}

/// Calculate atmospheric refraction correction in degrees
/// Δe = (P/1010)(283/(273+T))(1.02/(60 tan(e₀ + 10.3/(e₀ + 5.11))))
fn atmospheric_refraction_correction(pressure: f64, temperature: f64, e0: f64) -> f64 {
    if e0 < -1.0 {
        // Sun is well below horizon
        return 0.0;
    }

    let _e0_rad = e0.to_radians();
    let denominator = (e0 + 10.3 / (e0 + 5.11)).to_radians().tan();

    (pressure / 1010.0) * (283.0 / (273.0 + temperature)) * (1.02 / (60.0 * denominator))
}

/// Calculate topocentric elevation angle with atmospheric refraction in degrees
/// e = e₀ + Δe
fn topocentric_elevation_angle(e0: f64, delta_e: f64) -> f64 {
    e0 + delta_e
}

/// Calculate topocentric zenith angle in degrees
/// Θ = 90 - e
fn topocentric_zenith_angle(e: f64) -> f64 {
    90.0 - e
}

/// Calculate topocentric astronomers azimuth in degrees
/// Γ = atan2(sin(H'), cos(H')sin(φ) - tan(δ')cos(φ))
fn topocentric_astronomers_azimuth(h_prime: f64, delta_prime: f64, latitude: f64) -> f64 {
    let h_prime_rad = h_prime.to_radians();
    let delta_prime_rad = delta_prime.to_radians();
    let lat_rad = latitude.to_radians();

    let numerator = h_prime_rad.sin();
    let denominator = h_prime_rad.cos() * lat_rad.sin() - delta_prime_rad.tan() * lat_rad.cos();

    numerator.atan2(denominator).to_degrees()
}

/// Calculate topocentric azimuth angle (measured eastward from north) in degrees
/// Φ = Γ + 180 (mod 360)
fn topocentric_azimuth_angle(gamma: f64) -> f64 {
    (gamma + 180.0).rem_euclid(360.0)
}

// ============================================================================
// Main Solar Position Function
// ============================================================================

/// Calculate solar position for a given time and observer location
///
/// # Arguments
///
/// * `datetime` - UTC datetime
/// * `observer` - Observer location (latitude, longitude, elevation)
/// * `atmosphere` - Atmospheric conditions for refraction correction
///
/// # Returns
///
/// `Result<SolarPosition, DeltaTError>` - Solar position or error if date is out of valid range
///
pub fn calculate_solar_position(
    datetime: DateTime<Utc>,
    observer: Observer,
    atmosphere: Atmosphere,
) -> Result<SolarPosition, DeltaTError> {
    // Convert datetime to Julian Day
    let jd: JDay<JulianDay> = datetime.into();

    // Calculate Delta T and Julian Ephemeris Day
    let delta_t: DeltaT = datetime.try_into()?;
    let jde: JDay<JulianEphemerisDay> = (jd, delta_t).into();

    // Calculate Julian centuries and millennia
    let jc: JDay<JulianCentury> = jd.into();
    let jce: JDay<JulianEphemerisCentury> = jde.into();
    let jme: JDay<JulianEphemerisMiliennium> = jce.into();

    // Calculate heliocentric coordinates
    let radius = calculate_heliocentric_radius(jme);
    let helio_lon = calculate_heliocentric_longitude(jme);
    let helio_lat = calculate_heliocentric_latitude(jme);

    // Calculate geocentric coordinates
    let theta = calculate_geocentric_longitude(helio_lon);
    let beta = calculate_geocentric_latitude(helio_lat);

    // Calculate nutation
    let (delta_psi, delta_epsilon) = calculate_nutation(jce);

    // Calculate ecliptic obliquity
    let epsilon0 = calculate_mean_ecliptic_obliquity(jme);
    let epsilon = true_ecliptic_obliquity(epsilon0, delta_epsilon.degrees());

    // Calculate aberration and apparent sun longitude
    let delta_tau = aberration_correction(radius.value());
    let lambda = apparent_sun_longitude(theta, delta_psi.degrees(), delta_tau);

    // Calculate sidereal time
    let v0 = mean_sidereal_time(jd, jc);
    let v = apparent_sidereal_time(v0, delta_psi.degrees(), epsilon);

    // Calculate geocentric sun position
    let alpha = geocentric_sun_right_ascension(lambda, epsilon, beta);
    let delta = geocentric_sun_declination(lambda, epsilon, beta);

    // Calculate equation of time
    let m = sun_mean_longitude(jme);
    let eot = equation_of_time(m, alpha, delta_psi.degrees(), epsilon);

    // Calculate local hour angle
    let h = local_hour_angle(v, observer.longitude().value(), alpha);

    // Calculate parallax corrections using Observer methods
    let xi = equatorial_horizontal_parallax(radius.value());
    let x = observer.x_term();
    let y = observer.y_term();

    let delta_alpha = parallax_sun_right_ascension(x, xi, h, delta);
    let delta_prime = topocentric_sun_declination(delta, x, y, xi, delta_alpha, h);
    let h_prime = topocentric_local_hour_angle(h, delta_alpha);

    // Calculate topocentric elevation and zenith
    let e0 = topocentric_elevation_angle_without_atmosphere(
        observer.latitude().value(),
        delta_prime,
        h_prime,
    );
    let delta_e = atmosphere
        .refraction_correction(Angle::<Degrees>::new(e0))
        .value();
    let e = topocentric_elevation_angle(e0, delta_e);

    let theta_zenith = topocentric_zenith_angle(e);
    let theta0_zenith = topocentric_zenith_angle(e0);

    // Calculate azimuth
    let gamma = topocentric_astronomers_azimuth(h_prime, delta_prime, observer.latitude().value());
    let phi = topocentric_azimuth_angle(gamma);

    Ok(SolarPosition::new(
        Angle::<Degrees>::new(theta_zenith),
        Angle::<Degrees>::new(theta0_zenith),
        Angle::<Degrees>::new(e),
        Angle::<Degrees>::new(e0),
        Angle::<Degrees>::new(phi),
        eot,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_observer_validation() {
        // Valid observer
        assert!(Observer::new(40.0, -105.0, 1830.0).is_ok());

        // Invalid latitude
        assert!(Observer::new(91.0, 0.0, 0.0).is_err());
        assert!(Observer::new(-91.0, 0.0, 0.0).is_err());

        // Invalid longitude
        assert!(Observer::new(0.0, 181.0, 0.0).is_err());
        assert!(Observer::new(0.0, -181.0, 0.0).is_err());

        // Invalid elevation
        assert!(Observer::new(0.0, 0.0, -600.0).is_err());
    }

    #[test]
    fn test_atmosphere_validation() {
        // Valid atmosphere
        assert!(Atmosphere::new(1013.25, 15.0).is_ok());

        // Invalid pressure
        assert!(Atmosphere::new(0.0, 15.0).is_err());
        assert!(Atmosphere::new(2001.0, 15.0).is_err());

        // Invalid temperature
        assert!(Atmosphere::new(1013.25, -300.0).is_err());
        assert!(Atmosphere::new(1013.25, 150.0).is_err());
    }

    #[test]
    fn test_solar_position_calculation() {
        // Test at noon on January 15, 2024 in Boulder, CO
        let dt = DateTime::parse_from_rfc3339("2024-01-15T19:00:00Z") // 12:00 MST
            .unwrap()
            .with_timezone(&Utc);

        let observer = Observer::new(40.0, -105.0, 1830.0).unwrap();
        let atmosphere = Atmosphere::STANDARD;

        let pos = calculate_solar_position(dt, observer, atmosphere).unwrap();

        // Sun should be relatively high in the sky at noon
        assert!(pos.elevation().value() > 20.0);
        assert!(pos.elevation().value() < 50.0);

        // Zenith should be complement of elevation
        assert!((pos.zenith().value() + pos.elevation().value() - 90.0).abs() < 0.1);

        // Azimuth should be roughly south (180°) at solar noon
        assert!(pos.azimuth().value() > 150.0 && pos.azimuth().value() < 210.0);

        // Equation of time should be reasonable (within ±20 minutes)
        assert!(pos.equation_of_time().abs() < 20.0);

        // Should be daylight
        assert!(pos.is_daylight());
        assert!(!pos.is_night());
    }

    #[test]
    fn test_atmospheric_refraction() {
        let observer = Observer::new(0.0, 0.0, 0.0).unwrap();

        // At high elevation angles, refraction should be small
        let atmos = Atmosphere::STANDARD;
        let dt = DateTime::parse_from_rfc3339("2024-01-15T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);

        let pos = calculate_solar_position(dt, observer, atmos).unwrap();

        // Difference between refracted and non-refracted should be small
        let refraction_effect = pos.refraction_correction();
        assert!(refraction_effect.value() >= 0.0);
        assert!(refraction_effect.value() < 1.0); // Less than 1 degree
    }

    #[test]
    fn test_equation_of_time() {
        // Test that equation of time is within expected range
        let dt = DateTime::parse_from_rfc3339("2024-01-15T12:00:00Z")
            .unwrap()
            .with_timezone(&Utc);

        let observer = Observer::new(0.0, 0.0, 0.0).unwrap();

        let pos = calculate_solar_position(dt, observer, Atmosphere::default()).unwrap();

        // Equation of time varies throughout the year but should be within ±20 minutes
        assert!(pos.equation_of_time().abs() < 20.0);
    }

    #[test]
    fn test_twilight_detection() {
        // Create a position at different elevations to test twilight
        let pos_day = SolarPosition::new(
            Angle::<Degrees>::new(45.0),
            Angle::<Degrees>::new(45.5),
            Angle::<Degrees>::new(45.0),
            Angle::<Degrees>::new(45.5),
            Angle::<Degrees>::new(180.0),
            0.0,
        );
        assert!(pos_day.is_daylight());
        assert!(!pos_day.is_civil_twilight());

        let pos_civil = SolarPosition::new(
            Angle::<Degrees>::new(93.0),
            Angle::<Degrees>::new(93.5),
            Angle::<Degrees>::new(-3.0),
            Angle::<Degrees>::new(-3.5),
            Angle::<Degrees>::new(180.0),
            0.0,
        );
        assert!(!pos_civil.is_daylight());
        assert!(pos_civil.is_civil_twilight());

        let pos_night = SolarPosition::new(
            Angle::<Degrees>::new(110.0),
            Angle::<Degrees>::new(110.5),
            Angle::<Degrees>::new(-20.0),
            Angle::<Degrees>::new(-20.5),
            Angle::<Degrees>::new(180.0),
            0.0,
        );
        assert!(pos_night.is_night());
    }

    #[test]
    fn test_air_mass_calculation() {
        let pos = SolarPosition::new(
            Angle::<Degrees>::new(30.0),
            Angle::<Degrees>::new(30.5),
            Angle::<Degrees>::new(60.0),
            Angle::<Degrees>::new(60.5),
            Angle::<Degrees>::new(180.0),
            0.0,
        );
        let air_mass = pos.air_mass();
        assert!(air_mass.is_some());
        assert!(air_mass.unwrap() > 1.0);
        assert!(air_mass.unwrap() < 2.0);

        // Below horizon - no air mass
        let pos_below = SolarPosition::new(
            Angle::<Degrees>::new(95.0),
            Angle::<Degrees>::new(95.5),
            Angle::<Degrees>::new(-5.0),
            Angle::<Degrees>::new(-5.5),
            Angle::<Degrees>::new(180.0),
            0.0,
        );
        assert!(pos_below.air_mass().is_none());
    }
}
