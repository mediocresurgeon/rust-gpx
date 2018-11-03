pub mod great_circle;
pub mod vincenty;


use std::f64::consts::PI;


/// Converts degrees to radians.
pub fn to_radians(degrees: f64) -> f64 {
    (PI / 180.0) * degrees
}