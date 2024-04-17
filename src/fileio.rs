use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::fs::OpenOptions;

use num_complex::Complex64;

pub fn complex_vec_from_file(path: &str) -> Vec<Complex64> {
    // Takes a file path as an input, returns vector of Complex64
    // Expecting data points separated by whitespace
    // Each data point represented by two whitespace separated values
    // for real and imaginary part
    // If odd number of values in file, last complex part will be 0
    // Panics if path invalid or if file has non-decimal text

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

pub fn complex_vec_to_file(path: &str, vec: &Vec<Complex64>) {
    //Takes a vector of complex numbers as input, writes them to file specified by path
    //Each complex number will be written on its own line
    //Real and imaginary parts will be separated by whitespace
    //If path doesn't exist will create a file, otherwise will append to existing file
    //Panics if file cannot be opened (i.e. nonexistent directory in path, no permission)
    let mut file = match OpenOptions::new().create(true).append(true).open(path) {
        Ok(x) => x,
        Err(err) => panic!("Cannot open {}, {}", path, err),
    };
    let mut nums: String = "\n".to_string();
    for num in vec {
        nums += format!("{} {}\n", num.re, num.im).as_str();
    }
    match file.write_all(nums.as_bytes()) {
        Ok(_) => (),
        Err(err) => panic!("write output to file failed, {}", err),
    }
}