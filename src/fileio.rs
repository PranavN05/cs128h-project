use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

use num_complex::Complex64;

pub fn complex_vec_from_file(path: &str) -> Vec<Complex64> {
    // Takes a file path as an input, returns vector of Complex64
    // Expecting data points separated by whitespace
    // Each data point represented by two whitespace separated values
    // for real and imaginary part
    // If odd number of values in file, last complex part will be 0

    //Read String from File
    let mut file = match File::open(Path::new(path)) {
        Err(_) => panic!("couldn't open {}", path),
        Ok(file) => file,
    };

    let mut s = String::new();
    match file.read_to_string(&mut s) {
        Err(err) => panic!("couldn't read {}: {}", path, err),
        _ => (),
    }
    
    //Parse String into Vec<Complex64>
    let mut toreturn: Vec<Complex64> = Vec::new();
    let mut iter = s.split_ascii_whitespace();
    loop {
        let mut num = Complex64::new(0.0, 0.0);
        //real part
        match iter.next() {
            Some(x) => {
                match x.parse::<f64>() {
                    Ok(x) => num.re = x,
                    Err(_) => panic!("couldn't convert {} into a decimal value", x),
                }
            }
            None => break,
        }
        //imaginary part
        match iter.next() {
            Some(x) => {
                match x.parse::<f64>() {
                    Ok(x) => num.im = x,
                    Err(_) => panic!("couldn't convert {} into a decimal value", x),
                }
            }
            None => num.im = 0.0,
        }
        toreturn.push(num);
    }
    toreturn
}

