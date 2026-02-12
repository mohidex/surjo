use chrono::Utc;
use surjo::*;

fn main() {
    let angle_in_degrees = Angle::<Degrees>::new(180.0);
    println!("Angle: {:?}", angle_in_degrees);

    let angle_in_radians: Angle<Radians> = angle_in_degrees.into();
    println!("Angle: {:?}", angle_in_radians);

    // Convert back to degrees
    let converted_back: Angle<Degrees> = angle_in_radians.into();
    println!("Angle: {:?}", converted_back);

    let dt = Utc::now();
    println!("Current UTC time: {:?}", dt);
    let jd: JDay<JulianDay> = dt.into();
    println!("Julian Day: {:?}", jd);

    let deltat: DeltaT = dt.try_into().unwrap();
    println!("Delta T: {:?}", deltat);
}
