use cs128h_project::fft::fft_in_place;
use cs128h_project::fileio;
use std::env;

fn main() {
    //Get filename from command line
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Correct usage: cargo run -- <filename>\nEnsure <filename> is in input folder");
        std::process::exit(1);
    }
    let mut filepath = String::from("./input/");
    filepath += args[1].as_str();
    //Read input from file
    let mut inp = fileio::complex_vec_from_file(&filepath);
    //Perform fft
    fft_in_place(&mut inp);
    //Read output to file
    filepath = filepath.replace("./input/", "./output/");
    fileio::complex_vec_to_file(&filepath, &inp);
}