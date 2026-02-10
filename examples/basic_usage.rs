use surjo::*;

fn main() {
    let angle_in_degrees = Angle::<Degrees>::new(180.0);
    println!("Angle: {:?}", angle_in_degrees);

    let angle_in_radians: Angle<Radians> = angle_in_degrees.into();
    println!("Angle: {:?}", angle_in_radians);

    // Convert back to degrees
    let converted_back: Angle<Degrees> = angle_in_radians.into();
    println!("Angle: {:?}", converted_back);
}
