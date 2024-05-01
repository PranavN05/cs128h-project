### How to use this project
1. In your desired directory, run `git clone https://github.com/PranavN05/cs128h-project.git`
2. You will need a file of data points to perform the fft  
  If your inputs are real numbers only, separate the values by whitespace  
  If your inputs are complex, each data point should be formatted as a real part and imaginary part(even if 0.0), separated by whitespace  
  The data points themselves should be separated by whitespace  
4. Drag this file into the `./input/` directory
5. Run `cargo run -- <filename> <input_type>`  
    Filename is NOT a path - just the name of the file in `./input/`  
    If your file is formatted as only real values, input_type is `real`  
    else input_type is `complex`  
