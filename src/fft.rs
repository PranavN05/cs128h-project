use num_complex::Complex64;
use std::f64::consts::PI;

#[inline(always)]
fn bit_reversal(ind: usize, nbits: usize) -> usize {
    ind.reverse_bits() >> (usize::BITS as usize - nbits)
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

fn root_unity(N: usize) -> Vec<Complex64> {
    (0..(1 << (N - 1)))
        .map(|x| Complex64::cis(-2f64 * PI * (x as f64) / ((1 << N) as f64)))
        .collect()
}

fn fft2_butterfly(inp: &mut Vec<Complex64>, lev: usize, n: usize) {
    let root_unity = Complex64::cis(-2f64 * PI / ((1 << lev + 1) as f64));
    for i in 0..(1 << (n - lev - 1)) {
        let mut cur_root_unity = Complex64::new(1f64, 0f64);
        for j in 0..(1 << lev) {
            let a1 = inp[(1 << (lev + 1)) * i + j];
            let a2 = inp[(1 << (lev + 1)) * i + (1 << lev) + j];
            inp[(1 << (lev + 1)) * i + j] = a1 + cur_root_unity * a2;
            inp[(1 << (lev + 1)) * i + (1 << lev) + j] = a1 - cur_root_unity * a2;
            cur_root_unity *= root_unity;
        }
    }
}

fn fft2_butterfly_weights(
    inp: &mut Vec<Complex64>,
    weights: &Vec<Complex64>,
    lev: usize,
    n: usize,
) {
    for i in 0..(1 << (n - lev - 1)) {
        for j in 0..(1 << lev) {
            let a1 = inp[(1 << (lev + 1)) * i + j];
            let a2 = inp[(1 << (lev + 1)) * i + (1 << lev) + j];
            inp[(1 << (lev + 1)) * i + j] = a1 + weights[j * (1 << (n - lev - 1))] * a2;
            inp[(1 << (lev + 1)) * i + (1 << lev) + j] =
                a1 - weights[j * (1 << (n - lev - 1))] * a2;
        }
    }
}

fn fft4_butterfly(inp: &mut Vec<Complex64>, lev: usize, n: usize) {
    for i in 0..(1 << (n - lev - 2)) {
        for j in 0..(1 << lev) {
            // let a: Vec<Complex64> =
        }
    }
}

pub fn fft(inp: &Vec<Complex64>) -> Vec<Complex64> {
    assert!(inp.len() & (inp.len() - 1) == 0); // checks if power of 2
    let n = inp.len().ilog2() as usize; // inp.len() = 2^n
    let mut v: Vec<Complex64> = (0..inp.len()).map(|x| inp[bit_reversal(x, n)]).collect();
    // let w = root_unity(n);
    for i in 0..n {
        fft2_butterfly(&mut v, i, n);
        // fft2_butterfly_weights(&mut v, &w, i, n);
    }
    v
}

pub fn fft_in_place(inp: &mut Vec<Complex64>) {
    let n = inp.len().ilog2() as usize;
    for i in 0..(1 << n) {
        if i < bit_reversal(i, n) {
            inp.swap(i, bit_reversal(i, n));
        }
    }
    for i in 0..n {
        fft2_butterfly(inp, i, n);
    }
}
