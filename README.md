# surjo
Type-safe, high-precision solar position calculations based on the NREL SPA algorithm.

## Features

- **Type-Safe**: Uses types to prevent mixing degrees/radians, Julian days, and other units
- **High Precision**: Implements the full NREL SPA algorithm with sub-arcsecond accuracy
- **Validated Inputs**: Compile-time and runtime validation of coordinates and atmospheric conditions
- **Rich API**: Query solar position, twilight states, air mass, and more
- **Zero-Cost Abstractions**: Type safety with no runtime overhead

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
chrono = "0.4"
```

## Quick Start

```rust
use chrono::{DateTime, Utc};
use solar_position::{calculate_solar_position, Observer, Atmosphere};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a datetime
    let dt = DateTime::parse_from_rfc3339("2024-06-21T12:00:00Z")?
        .with_timezone(&Utc);
    
    // Define observer location (Boulder, Colorado)
    let observer = Observer::new(
        40.0,      // latitude (degrees, positive north)
        -105.0,    // longitude (degrees, positive east)
        1830.0,    // elevation (meters above sea level)
    )?;
    
    // Use standard atmospheric conditions
    let atmosphere = Atmosphere::STANDARD;
    
    // Calculate solar position
    let position = calculate_solar_position(dt, observer, atmosphere)?;
    
    // Access results with type-safe angles
    println!("Solar Elevation: {:.2}°", position.elevation().value());
    println!("Solar Azimuth: {:.2}°", position.azimuth().value());
    println!("Solar Zenith: {:.2}°", position.zenith().value());
    
    // Check if sun is above horizon
    if position.is_daylight() {
        println!("The sun is up!");
        
        // Get air mass (optical path length through atmosphere)
        if let Some(air_mass) = position.air_mass() {
            println!("Air mass: {:.2}", air_mass);
        }
    }
    
    Ok(())
}
```

## Core Concepts

### Type-Safe Angles

The library uses phantom types to distinguish between degrees and radians at compile time:

```rust
use solar_position::angle::{Angle, Degrees, Radians};

// Create an angle in degrees
let elevation: Angle<Degrees> = position.elevation();

// Access the raw value
let degrees: f64 = elevation.value();

// Convert to radians (type-safe conversion)
let elevation_rad: Angle<Radians> = elevation.into();
let radians: f64 = elevation_rad.value();

// This won't compile (type mismatch):
// let bad: Angle<Radians> = elevation; // Error!
```

### Observer Location

The `Observer` struct represents a location on Earth with validation:

```rust
use solar_position::Observer;

// Valid observer
let observer = Observer::new(40.0, -105.0, 1830.0)?;

// These will return errors:
let invalid_lat = Observer::new(91.0, 0.0, 0.0);   // Latitude out of range
let invalid_lon = Observer::new(0.0, 181.0, 0.0);  // Longitude out of range
let invalid_elev = Observer::new(0.0, 0.0, -600.0); // Elevation too low

// Access validated values
let lat: Angle<Degrees> = observer.latitude();
let lon: Angle<Degrees> = observer.longitude();
let elev: f64 = observer.elevation();
```

### Atmospheric Conditions

Control atmospheric refraction correction:

```rust
use solar_position::Atmosphere;

// Standard atmosphere (1013.25 mb, 15°C)
let standard = Atmosphere::STANDARD;

// Custom conditions
let custom = Atmosphere::new(
    1000.0,  // pressure in millibars
    25.0,    // temperature in Celsius
)?;

// Access values
println!("Pressure: {} mb", custom.pressure());
println!("Temperature: {} °C", custom.temperature());
```

### Solar Position Results

The `SolarPosition` struct provides comprehensive solar position data:

```rust
let position = calculate_solar_position(dt, observer, atmosphere)?;

// Angular positions (all return Angle<Degrees>)
let elevation = position.elevation();
let azimuth = position.azimuth();
let zenith = position.zenith();

// Versions without atmospheric refraction
let elev_no_refraction = position.elevation_no_refraction();
let zenith_no_refraction = position.zenith_no_refraction();

// How much refraction was applied
let refraction: Angle<Degrees> = position.refraction_correction();

// Equation of time (solar time vs clock time)
let eot: f64 = position.equation_of_time(); // minutes
```

### Twilight Detection

Built-in methods for detecting twilight states:

```rust
// Check sun position relative to horizon
if position.is_daylight() {
    println!("Sun is above horizon");
} else if position.is_civil_twilight() {
    println!("Civil twilight (0° to -6°)");
} else if position.is_nautical_twilight() {
    println!("Nautical twilight (-6° to -12°)");
} else if position.is_astronomical_twilight() {
    println!("Astronomical twilight (-12° to -18°)");
} else if position.is_night() {
    println!("Night (below -18°)");
}
```

### Air Mass Calculation

Get the atmospheric optical path length:

```rust
if let Some(air_mass) = position.air_mass() {
    println!("Air mass: {:.2}", air_mass);
    // Air mass = 1.0 at zenith, increases toward horizon
    // Returns None if sun is below horizon
}
```

## Type System Architecture

The library uses phantom types extensively for compile-time safety:

### Angle Types

```rust
pub struct Angle<U>(f64, PhantomData<U>);
pub struct Degrees;
pub struct Radians;

// Automatic conversion via From/Into traits
let deg: Angle<Degrees> = Angle::new(45.0);
let rad: Angle<Radians> = deg.into();
```

### Julian Day Types

```rust
pub struct JDay<T>(f64, PhantomData<T>);
pub struct JulianDay;
pub struct JulianEphemerisDay;
pub struct JulianCentury;
pub struct JulianEphemerisCentury;
pub struct JulianEphemerisMiliennium;

// Type-safe conversions ensure correct astronomical time scales
let jd: JDay<JulianDay> = datetime.into();
let jce: JDay<JulianEphemerisCentury> = jde.into();
```

### Heliocentric Coordinates

```rust
pub struct Heliocentric<T>(f64, PhantomData<T>);
pub struct HeliocentricLongitude;
pub struct HeliocentricLatitude;
pub struct HeliocentricRadius;

let longitude: Angle<Degrees> = calculate_heliocentric_longitude(jme);
let latitude: Angle<Degrees> = calculate_heliocentric_latitude(jme);
let radius: Heliocentric<HeliocentricRadius> = calculate_heliocentric_radius(jme);
```

## Advanced Usage

### Working with Different Time Scales

```rust
use solar_position::jday::{JDay, JulianDay, JulianEphemerisDay};
use solar_position::deltat::DeltaT;

let dt = DateTime::parse_from_rfc3339("2024-01-15T12:00:00Z")?
    .with_timezone(&Utc);

// Convert to Julian Day
let jd: JDay<JulianDay> = dt.into();

// Calculate Delta T (difference between TT and UT)
let delta_t: DeltaT = dt.try_into()?;

// Convert to Julian Ephemeris Day
let jde: JDay<JulianEphemerisDay> = (jd, delta_t).into();

println!("Delta T: {:.2} seconds", delta_t.seconds());
```

### Custom Atmospheric Models

```rust
// High altitude location with low pressure
let high_altitude = Atmosphere::new(700.0, 5.0)?;

// Sea level on a hot day
let sea_level_hot = Atmosphere::new(1013.25, 35.0)?;

// Compare refraction effects
let pos1 = calculate_solar_position(dt, observer, high_altitude)?;
let pos2 = calculate_solar_position(dt, observer, sea_level_hot)?;

let refraction1 = pos1.refraction_correction();
let refraction2 = pos2.refraction_correction();

println!("Refraction difference: {:.4}°", 
         (refraction1.value() - refraction2.value()).abs());
```

### Batch Calculations

```rust
use chrono::Duration;

let start = DateTime::parse_from_rfc3339("2024-06-21T00:00:00Z")?
    .with_timezone(&Utc);
let observer = Observer::new(40.0, -105.0, 1830.0)?;
let atmosphere = Atmosphere::STANDARD;

// Calculate solar position every hour for a day
for hour in 0..24 {
    let dt = start + Duration::hours(hour);
    let pos = calculate_solar_position(dt, observer, atmosphere)?;
    
    println!("Hour {}: Elevation = {:.2}°, Azimuth = {:.2}°",
             hour, pos.elevation().value(), pos.azimuth().value());
}
```

## Error Handling

The library uses Result types for all fallible operations:

```rust
use solar_position::{Observer, ObserverError, DeltaTError};

// Observer validation errors
match Observer::new(95.0, 0.0, 0.0) {
    Ok(obs) => println!("Valid observer"),
    Err(ObserverError::InvalidLatitude(lat)) => {
        println!("Invalid latitude: {}", lat);
    }
    Err(e) => println!("Other error: {}", e),
}

// Delta T calculation errors (year out of range)
let far_future = DateTime::parse_from_rfc3339("3500-01-01T12:00:00Z")?
    .with_timezone(&Utc);

match calculate_solar_position(far_future, observer, atmosphere) {
    Ok(pos) => println!("Position calculated"),
    Err(DeltaTError::YearOutOfRange(year)) => {
        println!("Year {} is outside valid range", year);
    }
}
```

## Accuracy

This implementation follows the NREL Solar Position Algorithm (SPA) which provides:

- **Uncertainty**: ±0.0003 degrees (≈1 arcsecond) for years 2000-6000
- **Valid range**: Years -1999 to 3000 (Delta T calculations)
- **Atmospheric refraction**: Accurate to ±0.01 degrees for elevation > 5°

## Performance

- **Zero-cost abstractions**: Phantom types compile away completely
- **Efficient calculations**: Uses pre-computed coefficient tables
- **No allocations**: All calculations on the stack
- **Fast conversions**: Type conversions are compile-time only

Typical performance on modern hardware:
- Single calculation: ~10-20 microseconds
- 1000 calculations: ~10-20 milliseconds

## Algorithm Details

The library implements the complete NREL SPA calculation chain:

1. **Time Conversion**: UTC → Julian Day → Julian Ephemeris Day
2. **Delta T**: Calculate TT-UT difference using polynomial approximations
3. **Heliocentric Position**: Earth's position relative to the Sun
4. **Geocentric Position**: Sun's position relative to Earth
5. **Nutation**: Corrections for Earth's axial precession
6. **Aberration**: Light-time correction
7. **Atmospheric Refraction**: Corrections for atmospheric bending
8. **Parallax**: Corrections for observer's position on Earth's surface
9. **Topocentric Position**: Final azimuth and elevation angles

### Mathematical Formulas

The following sections detail the mathematical formulas used at each step of the calculation.

#### 1. Time Conversions

**Julian Day (JD)** from Unix timestamp:
```
JD = unixtime / 86400 + 2440587.5
```

**Julian Ephemeris Day (JDE)** incorporates Delta T:
```
JDE = JD + ΔT / 86400
```

**Julian Century (JC)** from J2000.0 epoch:
```
JC = (JD - 2451545.0) / 36525
```

**Julian Ephemeris Century (JCE)**:
```
JCE = (JDE - 2451545.0) / 36525
```

**Julian Ephemeris Millennium (JME)**:
```
JME = JCE / 10
```

#### 2. Delta T Calculation

Delta T (ΔT) is the difference between Terrestrial Time (TT) and Universal Time (UT). Different polynomial approximations are used for different time periods.

For years 2005-2050:
```
y = year + (month - 0.5) / 12
ΔT = 62.92 + 0.32217(y - 2000) + 0.005589(y - 2000)²
```

For years 1986-2005:
```
u = y - 2000
ΔT = 63.86 + 0.3345u - 0.060374u² + 0.0017275u³ + 0.000651814u⁴ + 0.00002373599u⁵
```

See the code for complete polynomial formulas for all time periods.

#### 3. Earth Heliocentric Coordinates

**Heliocentric Longitude (L)**:
```
L = (L₀ + L₁·JME + L₂·JME² + L₃·JME³ + L₄·JME⁴ + L₅·JME⁵) / 10⁸
```

Where each Lᵢ term is:
```
Lᵢ = Σ Aᵢⱼ cos(Bᵢⱼ + Cᵢⱼ·JME)
```

**Heliocentric Latitude (B)**:
```
B = (B₀ + B₁·JME) / 10⁸
```

**Heliocentric Radius Vector (R)** in AU:
```
R = (R₀ + R₁·JME + R₂·JME² + R₃·JME³ + R₄·JME⁴) / 10⁸
```

All coefficient tables (L₀-L₅, B₀-B₁, R₀-R₄) are pre-computed constants.

#### 4. Geocentric Coordinates

**Geocentric Longitude (Θ)**:
```
Θ = L + 180° (mod 360°)
```

**Geocentric Latitude (β)**:
```
β = -B
```

#### 5. Nutation Calculations

First, calculate mean parameters:

**Mean Elongation of Moon (X₀)**:
```
X₀ = 297.85036 + 445267.111480·JCE - 0.0019142·JCE² + JCE³/189474
```

**Mean Anomaly of Sun (X₁)**:
```
X₁ = 357.52772 + 35999.050340·JCE - 0.0001603·JCE² - JCE³/300000
```

**Mean Anomaly of Moon (X₂)**:
```
X₂ = 134.96298 + 477198.867398·JCE + 0.0086972·JCE² + JCE³/56250
```

**Moon's Argument of Latitude (X₃)**:
```
X₃ = 93.27191 + 483202.017538·JCE - 0.0036825·JCE² + JCE³/327270
```

**Longitude of Moon's Ascending Node (X₄)**:
```
X₄ = 125.04452 - 1934.136261·JCE + 0.0020708·JCE² + JCE³/450000
```

Then calculate **Nutation in Longitude (Δψ)** and **Nutation in Obliquity (Δε)**:
```
For each term i in the nutation table:
  arg = Y₀ᵢX₀ + Y₁ᵢX₁ + Y₂ᵢX₂ + Y₃ᵢX₃ + Y₄ᵢX₄

Δψ = Σ[(aᵢ + bᵢ·JCE) sin(arg)] / 36000000
Δε = Σ[(cᵢ + dᵢ·JCE) cos(arg)] / 36000000
```

#### 6. Ecliptic Obliquity

**Mean Ecliptic Obliquity (ε₀)** in arcseconds:
```
U = JME / 10
ε₀ = 84381.448 - 4680.93U - 1.55U² + 1999.25U³ - 51.38U⁴ 
     - 249.67U⁵ - 39.05U⁶ + 7.12U⁷ + 27.87U⁸ + 5.79U⁹ + 2.45U¹⁰
```

**True Ecliptic Obliquity (ε)**:
```
ε = ε₀ + Δε
```

#### 7. Aberration and Apparent Sun Longitude

**Aberration Correction (Δτ)**:
```
Δτ = -20.4898 / (3600·R)
```

**Apparent Sun Longitude (λ)**:
```
λ = Θ + Δψ + Δτ
```

#### 8. Sidereal Time

**Mean Sidereal Time at Greenwich (ν₀)**:
```
ν₀ = 280.46061837 + 360.98564736629(JD - 2451545) + 0.000387933·JC² - JC³/38710000
```

**Apparent Sidereal Time (ν)**:
```
ν = ν₀ + Δψ cos(ε)
```

#### 9. Geocentric Sun Position

**Geocentric Right Ascension (α)**:
```
α = atan2(sin(λ)cos(ε) - tan(β)sin(ε), cos(λ))
```

**Geocentric Declination (δ)**:
```
δ = asin(sin(β)cos(ε) + cos(β)sin(ε)sin(λ))
```

#### 10. Equation of Time

**Sun Mean Longitude (M)**:
```
M = 280.4664567 + 360007.6982779·JME + 0.03032028·JME² 
    + JME³/49931 - JME⁴/15300 - JME⁵/2000000
```

**Equation of Time (E)** in minutes:
```
E = 4(M - 0.0057183 - α + Δψ cos(ε))
```

#### 11. Local Hour Angle

**Local Hour Angle (H)**:
```
H = ν + longitude - α
```

#### 12. Parallax Corrections

**Equatorial Horizontal Parallax (ξ)**:
```
ξ = 8.794 / (3600·R)
```

**Observer's geocentric position**:
```
u = atan(0.99664719 tan(φ))
x = cos(u) + (elevation/6378140)cos(φ)
y = 0.99664719 sin(u) + (elevation/6378140)sin(φ)
```

**Parallax in Sun Right Ascension (Δα)**:
```
Δα = atan2(-x sin(ξ) sin(H), cos(δ) - x sin(ξ) cos(H))
```

**Topocentric Sun Declination (δ')**:
```
δ' = atan2((sin(δ) - y sin(ξ)) cos(Δα), cos(δ) - x sin(ξ) cos(H))
```

**Topocentric Local Hour Angle (H')**:
```
H' = H - Δα
```

#### 13. Topocentric Position

**Topocentric Elevation Angle without Refraction (e₀)**:
```
e₀ = asin(sin(φ)sin(δ') + cos(φ)cos(δ')cos(H'))
```

**Atmospheric Refraction Correction (Δe)** in degrees:
```
Δe = (P/1010)(283/(273+T))(1.02/(60 tan(e₀ + 10.3/(e₀ + 5.11))))
```
Where P is pressure in millibars and T is temperature in Celsius.

**Topocentric Elevation Angle (e)**:
```
e = e₀ + Δe
```

**Topocentric Zenith Angle (θ)**:
```
θ = 90° - e
```

**Topocentric Azimuth Angle (Φ)** measured eastward from north:
```
Γ = atan2(sin(H'), cos(H')sin(φ) - tan(δ')cos(φ))
Φ = Γ + 180° (mod 360°)
```

#### 14. Air Mass Calculation

For solar applications, air mass is calculated using the Kasten-Young formula:
```
AM = 1 / (cos(θ) + 0.50572(96.07995 - θ)⁻¹·⁶³⁶⁴)
```

Where θ is the zenith angle in degrees. Valid only when the sun is above the horizon.

### Formula Notation

- **φ** = Observer latitude
- **λ** = Apparent sun longitude (or observer longitude in context)
- **ε** = True ecliptic obliquity
- **δ** = Declination
- **α** = Right ascension
- **θ** = Zenith angle
- **Φ** = Azimuth angle
- **H** = Hour angle
- **R** = Heliocentric radius (AU)
- **JME** = Julian Ephemeris Millennium
- **JCE** = Julian Ephemeris Century
- **ΔT** = Delta T (TT - UT)
- **Δψ** = Nutation in longitude
- **Δε** = Nutation in obliquity

All angles are in degrees unless otherwise specified. Intermediate calculations use radians where trigonometric functions are applied.

## References

- Reda, I.; Andreas, A. (2003). Solar Position Algorithm for Solar Radiation Applications. 
  NREL Report No. TP-560-34302. [NREL SPA](https://www.nrel.gov/docs/fy08osti/34302.pdf)

- Meeus, J. (1998). Astronomical Algorithms (2nd ed.). Willmann-Bell, Inc.

## Examples

See the `examples/` directory for complete examples:

- `basic.rs` - Simple solar position calculation
- `twilight.rs` - Detecting different twilight phases
- `daily_path.rs` - Plotting the sun's path through the day
- `equation_of_time.rs` - Calculating solar time corrections
- `air_mass.rs` - Computing atmospheric optical depth

## Contributing

Contributions are welcome! Please ensure:

- All tests pass (`cargo test`)
- Code is formatted (`cargo fmt`)
- No clippy warnings (`cargo clippy`)
- Documentation is updated

## Acknowledgments

- NREL for the Solar Position Algorithm
- Jean Meeus for Astronomical Algorithms
- The Rust community for excellent tools and libraries