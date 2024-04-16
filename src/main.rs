#[cfg(test)]
mod tests;
mod fileio;

use std::f64::consts::PI;

use num_complex::{Complex64, ComplexFloat};

fn bit_reversal(ind: usize, nbits: usize) -> usize {
    let mut out: usize = 0;
    for i in 0..nbits {
        out += ((ind >> (nbits - i - 1)) & 1) << i;
    }
    out
}


fn fft(inp: &Vec<Complex64>) -> Vec<Complex64> {
    assert!(inp.len() & (inp.len() - 1) == 0); // checks if power of 2
    let n = (inp.len() as f64).log2() as usize; // inp.len() = 2^n
    let mut v: Vec<Complex64> = (0..inp.len()).map(|x| inp[bit_reversal(x, n)]).collect();
    for i in 0..n {
        // n - 1 iterations
        // println!("{:?}", v);
        for j in 0..(1 << (n - i - 1)) {
            // 2^(n - i - 1) chunks for each smaller fft/butterfly
            let root_unity = Complex64::cis(-2f64 * PI / ((1 << (i + 1)) as f64)); // is it -2 or 2
            let mut cur_root_unity = Complex64::new(1f64, 0f64);
            let a = v[((1 << (i + 1)) * j)..((1 << (i + 1)) * j + (1 << i))].to_vec();
            let b =
                v[((1 << (i + 1)) * j + (1 << i))..((1 << (i + 1)) * j + (1 << (i + 1)))].to_vec();
            // println!("{:?} {:?}", a, b);
            for k in 0..(1 << i) {
                // we have 2^(i + 1) size of each chunk
                v[(1 << (i + 1)) * j + k] = a[k] + cur_root_unity * b[k];
                v[(1 << (i + 1)) * j + (1 << i) + k] = a[k] - cur_root_unity * b[k];
                cur_root_unity *= root_unity;
            }
        }
    }
    v
}

fn main() {
    println!("Hello, world!");
    let v: Vec<Complex64> = vec![
        Complex64::new(2f64, 0f64),
        Complex64::new(1f64, 0f64),
        Complex64::new(-1f64, 0f64),
        Complex64::new(5f64, 0f64),
        Complex64::new(0f64, 0f64),
        Complex64::new(3f64, 0f64),
        Complex64::new(0f64, 0f64),
        Complex64::new(-4f64, 0f64),
    ];
    let h = fft(&v);
    println!("{:?}", h);
}
