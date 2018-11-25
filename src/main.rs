// https://stackoverflow.com/questions/26946646/rust-package-with-both-a-library-and-a-binary#answer-50402684


extern crate serialization;
extern crate geodesics;


use serialization::io;
use geodesics::to_radians;
use geodesics::vincenty;


fn main() {
    let bothell = (to_radians( 47.8009363), to_radians(-122.231981));
    let albany  = (to_radians(-36.7259),    to_radians( 174.6936));

    let distance = vincenty::get_distance_wgs84(&bothell, &albany);

    match distance {
        Some(d) => println!("Distance: {}", d / 1000.0),
        None    => println!("Unable to calculate distance.")
    }
}