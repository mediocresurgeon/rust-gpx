/* Implementation of Thadeus Vincenty's algorithm,
 * originally published in 1975:
 * https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
 * 
 * Written with help from
 * https://nathanrooy.github.io/posts/2016-12-18/vincenty-formula-with-python/
 * and
 * https://www.movable-type.co.uk/scripts/latlong-vincenty.html
 * 
 * The notation in private functions is intended to 
 * match the notation presented by Vincenty.
 * As such, it may violate certain Rust conventions.
 * Public functions should still adhere to Rust conventions.
 */


use std::f64::EPSILON;


/// Finds the square of the cosine of an angle.
/// 
/// cos²α = 1 − sin²α
fn get_cos_sq_alpha(sin_alpha: f64) -> f64 {

    // Begin with trigonometric identity:
    // cos²α + sin²α = 1

    // make cos²α the subject:
    // cos²α = 1 − sin²α

    1.0 - sin_alpha.powf(2.0)
}


/// Returns the cosine of the reduced latitude on an ellipsoid.
/// 
/// `tan_U` - The tangent of the reduced latitude.
fn get_cos_U(tan_U: f64) -> f64 {

    // Start with trigonometric identity:
    // sec²θ = 1 + tan²θ
    
    // take the square root of both sides:
    // secθ = √(1 + tan²θ)

    // remember that secθ is the inverse of cosθ
    // taking the inverse of both sides gives:
    // cosθ = 1 / √(1 + tan²θ)

    1.0 / (1.0 + tan_U.powf(2.0)).sqrt()
}


/// Returns the sine of the reduced latitude on an ellipsoid.
/// 
/// `tan_U` - The tangent of the reduced latitude.
/// `cos_U` - The cosine of the reduced latitude.
fn get_sin_U(tan_U: f64, cos_U: f64) -> f64 {

    // Start with trigonometric identity:
    // tanθ = sinθ / cosθ

    // rearrange to get:
    // sinθ = tanθ cosθ

    tan_U * cos_U
}


/// Returns the tangent of the reduced latitude on an ellipsoid.
/// 
/// `f`   - The flattening of the ellipsoid.
/// `phi` - The latitude (in radians).
/// 
/// As shown in "Notation" section of Vincenty's paper.
/// tanU = (1 - f) tanφ
fn get_tan_U(f: f64, phi: &f64) -> f64 {
    (1.0 - f) * phi.tan()
}


/// Returns the flattening of an ellipsoid.
/// 
/// `a` - The major semi-axis of the ellipsoid.
/// `b` - The minor semi-axis of the ellipsoid.
/// 
/// As shown in "Notation" section of Vincenty's paper.
/// f, flattening = (a - b) / a
fn get_flattening(a: &f64, b: &f64) -> f64 {
    (a - b) / a
}


/// Equation (3) of Vincenty's direct formula.
/// 
/// A = 1 + (u² / 16384) {4096 + u² [-768 + u² (320 - 175u²)]}
fn get_A(u_sq: f64) -> f64 {

    // A = 1 + (u² / 16384) {4096 + u² [-768 + u² (320 - 175u²)]}
    //         |----------| |-----------------------------------|
    //            term4         term3
    //                                 |-----------------------|
    //                                     term2
    //                                            |-----------|
    //                                                term1

    let term1 = 320.0 - 175.0 * u_sq;
    let term2 = -768.0 + u_sq * term1;
    let term3 = 4096.0 + u_sq * term2;
    let term4 = u_sq / 16384.0;

    1.0 + term4 * term3
}


/// Equation (4) of Vincenty's direct formula.
/// 
/// B = (u² / 1024) {256 + u² [−128 + u² (74 − 47u²)]}
fn get_B(u_sq: f64) -> f64 {

    // B = (u² / 1024) {256 + u² [−128 + u² (74 − 47u²)]}
    //     |---------| |--------------------------------|
    //        term4       term3
    //                           |---------------------|
    //                               term2
    //                                      |---------|
    //                                         term1

    let term1 = 74.0 - 47.0 * u_sq;
    let term2 = -128.0 + u_sq * term1;
    let term3 = 256.0 + u_sq * term2;
    let term4 = u_sq / 1024.0;

    term4 * term3
}


/// Equation (6) of Vincenty's direct formula.
///
/// Δσ = B sinσ {cos2σm + (B / 4) [cosσ (−1 + 2cos²2σm) − (B / 6) cos2σm (−3 + 4sin²σ) (−3 + 4cos²2σm)]}
fn get_delta_sigma(B:             f64,
                   sin_sigma:     f64,
                   cos_sigma:     f64,
                   cos_2_sigma_m: f64) -> f64 {
    
    // Δσ = B sinσ {cos2σm + (B / 4) [cosσ (−1 + 2cos²2σm) − (B / 6) cos2σm (−3 + 4sin²σ) (−3 + 4cos²2σm)]}
    //      |----| |--------------------------------------------------------------------------------------|
    //       term8    term7
    //                       |-----| |-------------------------------------------------------------------|
    //                        term6   term5
    //                                     |-------------|   |-----|        |-----------| |-------------|
    //                                          term4         term3             term2          term1

    let cos_sq_2_sigma_m = cos_2_sigma_m.powf(2.0); // This term is used multiple times

    let term1 = -3.0 + 4.0 * cos_sq_2_sigma_m;
    let term2 = -3.0 + 4.0 * sin_sigma.powf(2.0);
    let term3 = B / 6.0;
    let term4 = -1.0 + 2.0 * cos_sq_2_sigma_m;
    let term5 = cos_sigma * term4 - term3 * cos_2_sigma_m * term2 * term1;
    let term6 = B / 4.0;
    let term7 = cos_2_sigma_m + term6 * term5;
    let term8 = B * sin_sigma;

    term8 * term7
}


/// Equation (10) of Vincenty's direct formula.
/// 
/// C = (f / 16) cos²α [4 + f (4 - 3 cos²α)]
fn get_C(f: f64, cos_sq_alpha: f64) -> f64 {

    // C = (f / 16) cos²α [4 + f (4 - 3 cos²α)]
    //     |------|       |-------------------|
    //      term3           term2
    //                           |-----------|
    //                               term1

    let term1 = 4.0 - 3.0 * cos_sq_alpha;
    let term2 = 4.0 + f * term1;
    let term3 = f / 16.0;

    term3 * cos_sq_alpha * term2
}


/// Equation (11) of Vincenty's direct formula
/// 
/// λ = L + (1 − C) f sinα {σ + C sinσ [cos2σm + C cosσ (−1 + 2cos²2σm)]}
/// 
/// Note: This function is modified to make λ the subject, not L.
fn get_lambda(f:             f64,
              L:             f64,
              C:             f64,
              sigma:         f64,
              sin_sigma:     f64,
              cos_sigma:     f64,
              sin_alpha:     f64,
              cos_2_sigma_m: f64,) -> f64 {

    // Begin with equation (11):
    // L = λ - (1 − C) f sinα {σ + C sinσ [cos2σm + C cosσ (−1 + 2cos²2σm)]}

    // Subtract (λ + L) from both sides;
    // -λ = -L - (1 − C) f sinα {σ + C sinσ [cos2σm + C cosσ (−1 + 2cos²2σm)]}

    // Multiply both sides by -1
    // λ = L + (1 − C) f sinα {σ + C sinσ [cos2σm + C cosσ (−1 + 2cos²2σm)]}
    //         |-----|        |--------------------------------------------|
    //          term4             term3
    //                                    |-------------------------------|
    //                                          term2
    //                                                     |-------------|
    //                                                          term1

    let term1 = -1.0 + 2.0 * cos_2_sigma_m.powf(2.0);
    let term2 = cos_2_sigma_m + C * cos_sigma * term1;
    let term3 = sigma + C * sin_sigma * term2;
    let term4 = 1.0 - C;

    L + term4 * f * sin_alpha * term3
}


/// Equation (14) of Vincenty's inverse formula.
/// 
/// sin²σ = (cosU₂ sinλ)² + (cosU₁ sinU₂ - sinU₁ cosU₂ cosλ)²
///
/// Note: Modified to return sinσ instead of sin²σ.
fn get_sin_sigma(sin_lambda: f64,
             cos_lambda: f64,
             sin_U1:     f64,
             sin_U2:     f64,
             cos_U1:     f64,
             cos_U2:     f64) -> f64 {

    // sin²σ = (cosU₂ sinλ)² + (cosU₁ sinU₂ - sinU₁ cosU₂ cosλ)²
    //         |----------|    |------------------------------|
    //            term1                      term2

    let term1 = cos_U2 * sin_lambda;
    let term2 = cos_U1 * sin_U2 - sin_U1 * cos_U2 * cos_lambda;

    (term1.powf(2.0) + term2.powf(2.0)).sqrt()
}


/// Equation (15) of Vincenty's inverse formula.
/// 
/// cosσ = sinU₁ sinU₂ + cosU₁ cosU₂ cosλ
fn get_cos_sigma(cos_lambda: f64,
             sin_U1:     f64,
             sin_U2:     f64,
             cos_U1:     f64,
             cos_U2:     f64) -> f64 {

    // cosσ = sinU₁ sinU₂ + cosU₁ cosU₂ cosλ
    //        |-------|     |--------------|
    //          term1             term2

    let term1 = sin_U1 * sin_U2;
    let term2 = cos_U1 * cos_U2 * cos_lambda;
    
    term1 + term2
}


/// Equation (16) of Vincenty's inverse formula.
/// 
/// tanσ = sinσ / cosσ
/// 
/// Note: This function has been modified to return σ instead of tanσ.
fn get_sigma(sin_sigma: f64, cos_sigma: f64) -> f64 {
    
    // Begin with trigonometric identity:
    // tanσ = sinσ / cosσ 

    // make σ the subject:
    // σ = tan⁻¹(sinσ / cosσ)

    // Consider the case where sinσ and cosσ are both negative.
    // In this case, it is the same as taking the ratio of their absolute values,
    // which gives an erroneous result.
    // Instead of using an atan function (which does not always preserve the quadrant),
    // we'll use atan2 (which *does* always preserve the quadrant.)

    sin_sigma.atan2(cos_sigma)
}


/// Equation (17) of Vincenty's inverse formula.
/// 
/// sinα = cosU₁ cosU₂ sinλ / sinσ
fn get_sin_alpha(sin_lambda: f64,
             sin_sigma:  f64,
             cos_U1:     f64,
             cos_U2:     f64) -> f64 {

    // sinα = cosU₁ cosU₂ sinλ / sinσ

    cos_U1 * cos_U2 * sin_lambda / sin_sigma
}


/// Equation (18) of Vincenty's inverse formula.
/// 
/// cos2σm = cosσ - 2sinU₁ sinU₂ / cos²α
fn get_cos_2_sigma_m(cos_sigma:    f64,
                     cos_sq_alpha: f64,
                     sin_U1:       f64,
                     sin_U2:       f64) -> f64 {

    // Vincenty notes that this equation will be indeterminate
    // over equatorial lines (when cos²α ≃ 0).
    // When this occurs, zero should be returned from this function.
    if cos_sq_alpha < EPSILON {
        0.0
    } else {

        // cos2σm = cosσ - 2sinU₁ sinU₂ / cos²α

        cos_sigma - 2.0 * sin_U1 * sin_U2 / cos_sq_alpha
    }
}


/// Equation (19) of Vincenty's inverse formula.
/// 
/// s = b A (σ − Δσ)
fn get_s(b:           &f64,
         A:           f64,
         sigma:       f64,
         delta_sigma: f64) -> f64 {
    
    // s = b A (σ − Δσ)

    b * A * (sigma - delta_sigma)
}


/// As shown in Vincenty's "Notation" section.
/// 
/// `cos_sq_alpha` - The square of the cosine of the azimuth at the equator
/// `a`            - The major semi-axis of the ellipsoid.
/// `b`            - The minor semi-axis of the ellipsoid.
/// 
/// u² = cos²α (a² - b²) / b²
/// 
/// Note: `u` should not be confused with `U`.
/// Also, `α` should not be confused with `a`.
fn get_u_sq(cos_sq_alpha: f64, a: &f64, b: &f64) -> f64 {

    // u² = cos²α (a² - b²) / b²

    let b_sq = b.powf(2.0); // term b² is used twice

    cos_sq_alpha * (a.powf(2.0) - b_sq) / b_sq
}


/// Intermediate values of interest
/// which are associated with
/// a particular lambda value.
struct ConvergedValues {
    sin_sigma:     f64,
    cos_sigma:     f64,
    sigma:         f64,
    cos_sq_alpha:  f64,
    cos_2_sigma_m: f64,
}


/// Calculates a new lambda value
/// and its associated intermediate values.
fn converge(f:      f64,
            L:      f64,
            lambda: f64,
            sin_U1: f64,
            sin_U2: f64,
            cos_U1: f64,
            cos_U2: f64) -> (f64, ConvergedValues) {

    // Trigonometry is expensive!
    // Remeber these values so we only have to compute them a single time.
    let sin_lambda = lambda.sin();
    let cos_lambda = lambda.cos();

    let sin_sigma = get_sin_sigma(sin_lambda, cos_lambda, sin_U1, sin_U2, cos_U1, cos_U2);
    let cos_sigma = get_cos_sigma(cos_lambda, sin_U1, sin_U2, cos_U1, cos_U2);
    let sigma     = get_sigma(sin_sigma, cos_sigma);

    let sin_alpha    = get_sin_alpha(sin_lambda, sin_sigma, cos_U1, cos_U2);
    let cos_sq_alpha = get_cos_sq_alpha(sin_alpha);

    let cos_2_sigma_m = get_cos_2_sigma_m(cos_sigma, cos_sq_alpha, sin_U1, sin_U2);
    let C             = get_C(f, cos_sq_alpha);

    let new_lambda = get_lambda(f, L, C, sigma, sin_sigma, cos_sigma, sin_alpha, cos_2_sigma_m);

    (new_lambda, ConvergedValues { 
        sin_sigma,
        cos_sigma,
        sigma,
        cos_sq_alpha,
        cos_2_sigma_m
    })
}


/// Returns the distance between two locations on the surface of an ellipsoid.
/// The distance is in the same units as the major and minor semi-axis.
/// In some cases (such as antipodal points) it may not be possible to produce a result.
/// 
/// For Earth, consider calling `get_distance_wgs84()` instead.
/// 
/// `(major_semi_axis, minor_semi_axis)` - The major and minor semi-axis which define the shape of the ellipsoid.
/// `(lat1, lat1)`                       - The latitude and longitude of the first location.  Must be in radians.
/// `(lat2, lat2)`                       - The latitude and longitude of the second location.  Must be in radians.
/// `tolerance`                          - Allowable error. Larger values require fewer operations; smaller values give more precise results.
pub fn get_distance((major_semi_axis, minor_semi_axis): &(f64, f64),
                    (lat1, lon1):                       &(f64, f64),
                    (lat2, lon2):                       &(f64, f64),
                    tolerance:                          f64) -> Option<f64> {

    // f is the "flattening" of the ellipsoid, which defines its shape.
    let f = get_flattening(major_semi_axis, minor_semi_axis);

    // L is the difference in longitude between the two points.
    let L = lon2 - lon1;

    // cache some intermediate values which we'll use later
    let tan_U1 = get_tan_U(f, lat1);
    let cos_U1 = get_cos_U(tan_U1);
    let sin_U1 = get_sin_U(tan_U1, cos_U1);

    let tan_U2 = get_tan_U(f, lat2);
    let cos_U2 = get_cos_U(tan_U2);
    let sin_U2 = get_sin_U(tan_U2, cos_U2);

    // Our first estimate for the distance will be the difference in longitude.
    // As we start looping to find better and better estimates,
    // we will (hopefully) converge on a solution.
    let mut lambda = L;
    
    // Put a cap on the number of loops we perform.
    // In general, fewer than 6 loops are required.
    let mut iterations = 0_u8;
    let intermediate_values = loop {

        // Remember what the old lambda value is so we can compare against it later
        let prev_lambda = lambda;

        // Calculate a new lambda and its associated intermediate values
        let (new_lambda, iv) = converge(f, L, prev_lambda, sin_U1, sin_U2, cos_U1, cos_U2);
        lambda = new_lambda;

        // did our lambda values converge?
        // if they did, we have  eveyrthing we need to exit the loop!
        if (prev_lambda - new_lambda).abs() < tolerance {
            break iv;
        } else {
            iterations += 1;
            // If the results have failed to converge, return None.
            // Vincenty notes that this can happen for nearly antipodal points.
            if iterations > 100 {
                return None
            }
        }
    };

    let u_sq        = get_u_sq(intermediate_values.cos_sq_alpha, major_semi_axis, minor_semi_axis);
    let A           = get_A(u_sq);
    let B           = get_B(u_sq);
    let delta_sigma = get_delta_sigma(B, intermediate_values.sin_sigma, intermediate_values.cos_sigma, intermediate_values.cos_2_sigma_m);
    let s           = get_s(minor_semi_axis, A, intermediate_values.sigma, delta_sigma);

    Some(s)
}


/// Returns the distance (in meters)
/// between two locations on Earth
/// as defined by the WGS-84 standard.
/// In some cases (such as antipodal points)
/// it may not be possible to produce a result.
/// 
/// `(lat1, lat1)` - The latitude and longitude of the first location.  Must be in radians.
/// `(lat2, lat2)` - The latitude and longitude of the second location.  Must be in radians.
pub fn get_distance_wgs84((lat1, lon1): &(f64, f64),
                          (lat2, lon2): &(f64, f64)) -> Option<f64> {
    
    const WGS84:     (f64, f64) = (6_378_137.0, 6_356_752.314_245);
    const TOLERANCE:  f64       = 0.000_000_000_001; // 1mm - 6mm on Earth

    let result = get_distance(&WGS84, &(*lat1, *lon1), &(*lat2, *lon2), TOLERANCE);

    match result {
        Some(distance) => {
            // Rounds the result to the nearest millimeter,
            // since accuracy beyond that is not guaranteed.
            Some((distance * 1000.0).round() / 1000.0)
        },
        None => None
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use super::super::to_radians;

    const WGS84:            (f64, f64) = (6_378_137.0, 6_356_752.314_245);
    const WGS84_FLATTENING:  f64       = 1.0 / 298.257_223_563;


    #[test]
    /// Ensures flattening() matches results from WGS84 standard.
    fn flattening_earth() {
        // Arrange
        let (major, minor) = WGS84;
        
        // Act
        let f = get_flattening(&major, &minor);

        // Assert
        assert!((f - WGS84_FLATTENING).abs() < 0.000_000_000_000_1); // Accurate to 13 decimal places
    }


    #[test]
    /// Ensures that the distance from
    /// the north pole to the equator
    /// is calculated correctly.
    fn north_pole_to_equator() {
        // Arrange
        let p1: (f64, f64) = (to_radians(90.0), 0.0);
        let p2: (f64, f64) = (            0.0,  0.0);

        // Act
        let distance = get_distance_wgs84(&p1, &p2);

        // Assert
        assert_eq!(distance.unwrap(), 10_001_965.729);
    }


    #[test]
    /// Ensures that the distance a quarter turn around the equator
    /// is calculated correctly.
    fn quarter_around_equator() {
        // Arrange
        let p1: (f64, f64) = (0.0, to_radians(90.0));
        let p2: (f64, f64) = (0.0,             0.0);

        // Act
        let distance = get_distance_wgs84(&p1, &p2);

        // Assert
        assert_eq!(distance.unwrap(), 10_018_754.171);
    }


    #[test]
    /// Ensures that the distance from
    /// 45N 45E to 45S 45W
    /// is calculated correctly.
    fn fourty_five_degrees() {
        // Arrange
        let p1: (f64, f64) = (to_radians( 45.0), to_radians( 45.0));
        let p2: (f64, f64) = (to_radians(-45.0), to_radians(-45.0));

        // Act
        let distance = get_distance_wgs84(&p1, &p2);

        // Assert
        assert_eq!(distance.unwrap(), 13_324_945.436);
    }


    #[test]
    /// Ensures that very, very small distances are calculated correctly.
    /// (This is a typical use-case.)
    fn small_angle() {
        // Arrange
        let p1: (f64, f64) = (to_radians(39.152501), to_radians(-84.412977));
        let p2: (f64, f64) = (to_radians(39.152505), to_radians(-84.412946));

        // Act
        let distance = get_distance_wgs84(&p1, &p2);

        // Assert
        assert_eq!(distance.unwrap(), 2.716);
    }
}