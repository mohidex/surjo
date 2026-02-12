/// Errors that can occur during Delta T calculations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeltaTError {
    /// Year is outside the valid range (-1999 to 3000)
    YearOutOfRange(i32),
}

impl std::fmt::Display for DeltaTError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DeltaTError::YearOutOfRange(year) => {
                write!(
                    f,
                    "Year {} is outside valid range (-1999 to 3000). \
                     Delta T calculations are not accurate for this year.",
                    year
                )
            }
        }
    }
}

impl std::error::Error for DeltaTError {}
