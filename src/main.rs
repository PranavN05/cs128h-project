use std::f32::consts::PI;

use num_complex::{Complex32, ComplexFloat};

fn bit_reversal(ind: usize, nbits: usize) -> usize {
    let mut out: usize = 0;
    for i in 0..nbits {
        out += ((ind >> (nbits - i - 1)) & 1) << i;
    }
    out
}

fn fft_butterfly(inp: Vec<Complex32>) -> Vec<Complex32> {
    todo!()
}

fn fft(inp: Vec<Complex32>) -> Vec<Complex32> {
    assert!(inp.len() & (inp.len() - 1) == 0); // checks if power of 2
    let n = (inp.len() as f32).log2() as usize; // inp.len() = 2^n
    let mut v: Vec<Complex32> = (0..inp.len()).map(|x| inp[bit_reversal(x, n)]).collect();
    for i in 0..n {
        // n - 1 iterations
        // println!("{:?}", v);
        for j in 0..(1 << (n - i - 1)) {
            // 2^(n - i - 1) chunks for each smaller fft/butterfly
            let root_unity = Complex32::cis(-2f32 * PI / ((1 << (i + 1)) as f32)); // is it -2 or 2
            let mut cur_root_unity = Complex32::new(1f32, 0f32);
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
    let v: Vec<Complex32> = vec![
        Complex32::new(2f32, 0f32),
        Complex32::new(1f32, 0f32),
        Complex32::new(-1f32, 0f32),
        Complex32::new(5f32, 0f32),
        Complex32::new(0f32, 0f32),
        Complex32::new(3f32, 0f32),
        Complex32::new(0f32, 0f32),
        Complex32::new(-4f32, 0f32),
    ];
    let h = fft(v);
    println!("{:?}", h);
}
