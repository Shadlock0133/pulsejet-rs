use std::fs;

use pulsejet_rs::{decode, Std};

fn main() {
    eprintln!("hello");
    let pulsejet_stream = fs::read("encoded.pj").unwrap();
    eprintln!("loaded encoded");
    let decoded = unsafe { decode::<Std>(pulsejet_stream.as_ptr()) };
    eprintln!("encoded decoded");
    let raw = decoded
    .into_iter()
    .flat_map(|x| x.to_le_bytes())
    .collect::<Vec<_>>();
    eprintln!("decoded transcoded");
    fs::write("rust.raw", raw).unwrap();
    eprintln!("transcoded saved");
}
