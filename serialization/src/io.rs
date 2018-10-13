extern crate serde;
extern crate serde_xml_rs;


#[derive(Debug, Deserialize)]
#[serde(rename = "gpx")]
pub struct WaypointsRoot {
    #[serde(rename = "wpt", default)]
    pub waypoints: Vec<Waypoint>
}


#[derive(Debug, Deserialize)]
#[serde(rename = "wpt")]
pub struct Waypoint {
    #[serde(rename="lat")]
    pub latitude:  f64,
    #[serde(rename="lon")]
    pub longitude: f64,
    pub name:      String,
}


pub fn do_something() {
    println!("Hello from io!")
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        // Arrange
        let foo = r##"
            <gpx version="1.1" creator="Xcode">
                <wpt lat="-36.7328" lon="174.7006">
                    <name>Aperture Memorial</name>
                    <time>2018-09-01T00:00:00Z</time>
                </wpt>
                <wpt lat="-36.7331" lon="174.7011">
                    <name>The Atrium - Massey University</name>
                    <time>2018-09-01T00:00:20Z</time>
                </wpt>
            </gpx>
        "##;
        
        // Act
        let result: WaypointsRoot = serde_xml_rs::from_str(foo).unwrap();

        // Assert
        println!("{:?}", result);
        assert_eq!(result.waypoints.len(), 2);
        
        assert_eq!(result.waypoints[0].name, "Aperture Memorial");
        assert_eq!(result.waypoints[0].latitude, -36.7328);
        assert_eq!(result.waypoints[0].longitude, 174.7006);
        
        assert_eq!(result.waypoints[1].name, "The Atrium - Massey University");
        assert_eq!(result.waypoints[1].latitude, -36.7331);
        assert_eq!(result.waypoints[1].longitude, 174.7011);
    }
}