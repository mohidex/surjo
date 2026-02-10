use std::fmt;

#[derive(Debug, Clone, Copy)]
pub enum AngleError {
    InvalidDegrees(f64),
    InvalidRadians(f64),
}

impl fmt::Display for AngleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AngleError::InvalidDegrees(val) => {
                write!(
                    f,
                    "Invalid degree value: {}. Must be in range [0, 360]",
                    val
                )
            }
            AngleError::InvalidRadians(val) => {
                write!(f, "Invalid radian value: {}. Must be in range [0, 2π]", val)
            }
        }
    }
}

impl std::error::Error for AngleError {}
