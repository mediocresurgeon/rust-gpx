/// Calculates the haversine of an angle (in radians).
fn haversine(radians: f64) -> f64 {
    // #1 Divide the angle by two.
    // #2 Get the sin of #1.
    // #3 Square the result of #2.
    // #4 Return the result of #3.
    (radians / 2.0).sin()
                   .powf(2.0)
}


/// Returns the angle (in radians) between two coordinates on a sphere.
/// The latitude and longitute must be in radians (not degrees).
/// 
/// # Hint:
/// Multiply the returned value by the radius of the sphere to get the distance between the two points.
pub fn get_angle((lat1, lon1): &(f64, f64),
                 (lat2, lon2): &(f64, f64)) -> f64 {
    // Pure formula:
    // distance = 2 radius arcsin( sqrt( hav(φ2 - φ1) + cos(φ1) * cos(φ2) * hav(λ2 - λ1) ) )
    // https://en.wikipedia.org/wiki/Haversine_formula#The_haversine_formula
    let cos_lat1 = lat1.cos();
    let cos_lat2 = lat2.cos();
    let hav_delta_lat = haversine(lat2 - lat1);
    let hav_delta_lon = haversine(lon2 - lon1);

    // Computers a better at calculating small distances if atan2 is used instead of asin sqrt
    // This requires intermediate steps and cannot be written as cleanly.
    // http://mathforum.org/library/drmath/view/51879.html

    // ir is an intermediate result (so we don't have to do this calculation more than once)
    let ir = hav_delta_lat + (cos_lat1 * cos_lat2 * hav_delta_lon);
    let y = ir.sqrt();
    let x = (1.0 - ir).sqrt();
    2.0 * y.atan2(x)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;


    // Tests should have a tolerance of +/- 1cm on a sphere the size of Earth
    const ANGLE_OF_ONE_CM_ON_EARTH: f64 = 0.000000002;


    #[test]
    /// Tests that the angle between the north pole and equatorial prime meridian is a right angle.
    fn get_angle_northpole_pmequator() {
        // Arrange
        let north_pole:             (f64, f64) = (PI/2.0, 0.0);
        let equator_prime_medidian: (f64, f64) = (   0.0, 0.0);
        
        // Act
        let result = get_angle(&north_pole, &equator_prime_medidian);
        
        // Assert
        let exact_answer = PI/2.0;
        assert!((result - exact_answer).abs() < ANGLE_OF_ONE_CM_ON_EARTH);
    }


    #[test]
    /// Tests that the angle between 45N 45E and 45S 45W is 120deg
    fn get_angle_pmequator_90equator() {
        // Arrange
        let n_e: (f64, f64) = ( PI/4.0,  PI/4.0); //45N, 45E
        let s_w: (f64, f64) = (-PI/4.0, -PI/4.0); //45S, 45W
        
        // Act
        let result = get_angle(&n_e, &s_w);
        
        // Assert
        let exact_answer = (2.0/3.0)* PI;
        assert!((result - exact_answer).abs() < ANGLE_OF_ONE_CM_ON_EARTH);
    }
}