// https://stackoverflow.com/questions/26946646/rust-package-with-both-a-library-and-a-binary#answer-50402684

extern crate serialization;


fn main() {
    println!("Hello, world!");
    serialization::io::do_something()
}
