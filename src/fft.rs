use num_complex::Complex64;
use std::f64::consts::PI;

fn bit_reversal(ind: usize, nbits: usize) -> usize {
    let mut out: usize = 0;
    for i in 0..nbits {
        out += ((ind >> (nbits - i - 1)) & 1) << i;
    }
    out
}

fn fft_butterfly(mut inp: Vec<Complex64>, rad: usize, lev: usize) -> Vec<Complex64> {
    let n = inp.len().ilog2() as usize;
    for i in 0..(1 << (n - lev - rad)) {
        for j in 0..(1 << lev) {
            let ind: Vec<usize> = (0..(1 << rad))
                .map(|x| i * (1 << (lev + rad)) + x * (1 << lev) + j)
                .collect();
            let a: Vec<Complex64> = ind
                .iter()
                .map(|x| {
                    let temp = inp[*x];
                    inp[*x] = Complex64::new(0f64, 0f64);
                    temp
                })
                .collect();
            for k in 0..(1 << rad) {
                let root_unity = Complex64::cis(
                    -2f64 * PI * bit_reversal(k, rad) as f64 / (1 << (lev + rad)) as f64,
                );
                let mut cur_root_unity = Complex64::new(1f64, 0f64);
                for l in 0..(1 << rad) {
                    inp[ind[l]] += a[k] * cur_root_unity;
                    cur_root_unity *= root_unity;
                }
            }
        }
    }
    inp
}

pub fn fft(inp: &Vec<Complex64>) -> Vec<Complex64> {
    assert!(inp.len() & (inp.len() - 1) == 0); // checks if power of 2
    let n = inp.len().ilog2() as usize; // inp.len() = 2^n
    let mut v: Vec<Complex64> = (0..inp.len()).map(|x| inp[bit_reversal(x, n)]).collect();
    for i in 0..n {
        // n - 1 iterations
        // println!("{:?}", v);
        // v = fft_butterfly(v, 2, i * 2);
        let root_unity = Complex64::cis(-2f64 * PI / ((1 << (i + 1)) as f64)); // is it -2 or 2
        for j in 0..(1 << (n - i - 1)) {
            // 2^(n - i - 1) chunks for each smaller fft/butterfly
            let mut cur_root_unity = Complex64::new(1f64, 0f64);
            // println!("{:?} {:?}", a, b);
            for k in 0..(1 << i) {
                // we have 2^(i + 1) size of each chunk
                let a = v[(1 << (i + 1)) * j + k];
                let b = v[(1 << (i + 1)) * j + (1 << i) + k];

                v[(1 << (i + 1)) * j + k] = a + cur_root_unity * b;
                v[(1 << (i + 1)) * j + (1 << i) + k] = a - cur_root_unity * b;
                cur_root_unity *= root_unity;
            }
        }
    }
    v
}
