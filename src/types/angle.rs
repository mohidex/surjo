use std::f64::consts::PI;
use std::marker::PhantomData;

use super::error::AngleError;

/// Marker type for radians
#[derive(Debug, Clone, Copy)]
pub struct Radians;

/// Marker type for degrees
#[derive(Debug, Clone, Copy)]
pub struct Degrees;

/// Generic Angle struct that can represent either degrees or radians
#[derive(Debug, Clone, Copy)]
pub struct Angle<U>(f64, PhantomData<U>);
impl<U> Angle<U> {
    pub fn new(value: f64) -> Self {
        Self(value, PhantomData)
    }

    pub fn value(&self) -> f64 {
        self.0
    }
}

/// Type conversion between Angle<Degrees> and Angle<Radians>
impl From<Angle<Degrees>> for Angle<Radians> {
    fn from(degrees: Angle<Degrees>) -> Self {
        let radians_value = degrees.value() * PI / 180.0;
        Self::new(radians_value)
    }
}

impl From<Angle<Radians>> for Angle<Degrees> {
    fn from(radians: Angle<Radians>) -> Self {
        let degrees_value = radians.value() * 180.0 / PI;
        Self::new(degrees_value)
    }
}

/// TryFrom implementation to validate degree values in [0, 360]
impl TryFrom<f64> for Angle<Degrees> {
    type Error = AngleError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if (0.0..=360.0).contains(&value) {
            Ok(Self(value, PhantomData))
        } else {
            Err(AngleError::InvalidDegrees(value))
        }
    }
}

/// TryFrom implementation to validate radian values in [0, 2π]
impl TryFrom<f64> for Angle<Radians> {
    type Error = AngleError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if (0.0..=2.0 * PI).contains(&value) {
            Ok(Self(value, PhantomData))
        } else {
            Err(AngleError::InvalidRadians(value))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryInto;

    #[test]
    fn test_degree_validation() {
        // Valid boundary values
        let deg: Result<Angle<Degrees>, _> = 0.0.try_into();
        assert!(deg.is_ok());
        let deg: Result<Angle<Degrees>, _> = 180.0.try_into();
        assert!(deg.is_ok());
        let deg: Result<Angle<Degrees>, _> = 360.0.try_into();
        assert!(deg.is_ok());

        // Invalid values
        let deg: Result<Angle<Degrees>, _> = (-1.0).try_into();
        assert!(deg.is_err());
        matches!(deg, Err(AngleError::InvalidDegrees(-1.0)));
        assert!(
            deg.err()
                .expect("REASON")
                .to_string()
                .contains("Invalid degree value")
        );

        let deg: Result<Angle<Degrees>, _> = 361.0.try_into();
        assert!(deg.is_err());
        matches!(deg, Err(AngleError::InvalidDegrees(361.0)));
        assert!(
            deg.err()
                .expect("REASON")
                .to_string()
                .contains("Invalid degree value")
        );
    }

    #[test]
    fn test_radian_validation() {
        // Valid boundary values
        let rad: Result<Angle<Radians>, _> = 0.0.try_into();
        assert!(rad.is_ok());
        let rad: Result<Angle<Radians>, _> = PI.try_into();
        assert!(rad.is_ok());
        let rad: Result<Angle<Radians>, _> = (2.0 * PI).try_into();
        assert!(rad.is_ok());

        // Invalid values
        let rad: Result<Angle<Radians>, _> = (-0.1).try_into();
        assert!(rad.is_err());
        matches!(rad, Err(AngleError::InvalidRadians(-0.1)));
        assert!(
            rad.err()
                .expect("REASON")
                .to_string()
                .contains("Invalid radian value")
        );

        let rad: Result<Angle<Radians>, _> = (2.0 * PI + 0.1).try_into();
        assert!(rad.is_err());
        matches!(rad, Err(AngleError::InvalidRadians(v)) if (v - (2.0 * PI + 0.1)).abs() < 1e-10);
        assert!(
            rad.err()
                .expect("REASON")
                .to_string()
                .contains("Invalid radian value")
        );
    }

    #[test]
    fn test_angle_conversions() {
        // Degrees to Radians
        let deg_90: Angle<Degrees> = 90.0.try_into().unwrap();
        let rad_from_deg: Angle<Radians> = deg_90.into();
        assert!((rad_from_deg.value() - PI / 2.0).abs() < 1e-10);

        // Radians to Degrees
        let rad_pi_4: Angle<Radians> = (PI / 4.0).try_into().unwrap();
        let deg_from_rad: Angle<Degrees> = rad_pi_4.into();
        assert!((deg_from_rad.value() - 45.0).abs() < 1e-10);

        // Full rotation
        let deg_360: Angle<Degrees> = 360.0.try_into().unwrap();
        let rad_2pi: Angle<Radians> = deg_360.into();
        assert!((rad_2pi.value() - 2.0 * PI).abs() < 1e-10);
    }
}
