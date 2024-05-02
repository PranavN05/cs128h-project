use cs128h_project::czt::czt;
use cs128h_project::fft::{fft_optimized, precompute_weights};
use cs128h_project::fileio;
use std::env;

fn main() {
    //Get filename from command line
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 || (args[2] != "real" && args[2] != "complex") {
        eprintln!("Correct usage: cargo run -- <filename> <type>");
        eprintln!("Ensure <filename> is in input folder");
        eprintln!("<type> must be either \"real\" or \"complex\"");
        std::process::exit(1);
    }
    let mut filepath = String::from("./input/");
    filepath += args[1].as_str();
    //Read input from file
    let mut inp = match args[2].as_str() {
        "real" => fileio::complex_vec_from_real_file(&filepath),
        "complex" => fileio::complex_vec_from_complex_file(&filepath),
        _ => panic!("Should not reach this"),
    };
    //Perform fft
    match inp.len() & inp.len() - 1 == 0 {
        true => {
            let weights = precompute_weights(inp.len().ilog2() as usize);
            fft_optimized(&mut inp, &weights);
        }
        false => {
            inp = czt(inp.clone(), inp.len(), None, None);
        }
    }

    //Read output to file
    filepath = filepath.replace("./input/", "./output/");
    fileio::complex_vec_to_file(&filepath, &inp);
}
