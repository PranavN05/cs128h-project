use num_complex::Complex64;
use std::f64::consts::PI;

#[inline(always)]
fn bit_reversal(ind: usize, nbits: usize) -> usize {
    ind.reverse_bits() >> ((usize::BITS as usize - nbits) % usize::BITS as usize)
}

pub fn precompute_weights(n: usize) -> Vec<Complex64> {
    (0..(1 << n))
        .map(|x| Complex64::cis(-2f64 * PI * (x as f64) / ((1 << n) as f64)))
        .collect()
}

fn fft_butterfly(inp: &mut Vec<Complex64>, w: &Vec<Complex64>, rad: usize, lev: usize, n: usize) {
    for (k0, chunk) in inp.chunks_mut(1 << (n - lev)).enumerate() {
        for k2 in 0..(1 << (n - lev - rad)) {
            let a: Vec<Complex64> = chunk
                .iter()
                .skip(k2)
                .step_by(1 << (n - lev - rad))
                .map(|x| *x)
                .collect();
            for k1 in 0..(1 << rad) {
                let mut sum = Complex64::new(0f64, 0f64);
                for j1 in 0..(1 << rad) {
                    sum += a[j1]
                        * w[j1 * bit_reversal(k1, rad) * (1 << (n - rad)) % (1 << n)]
                        * w[j1 * bit_reversal(k0, lev) * (1 << (n - lev - rad)) % (1 << n)];
                }
                chunk[(1 << (n - lev - rad)) * k1 + k2] = sum;
            }
        }
    }
}

fn fft2_butterfly(inp: &mut Vec<Complex64>, lev: usize) {
    //n = log(inp.len)
    //lev = log(size of chunk) - 1
    //root_unity = factor between "twiddle factors" = e^(-2pi/(size of chunk)) = e^(-2pi / 2^(lev+1))
    let root_unity = Complex64::cis(-2f64 * PI / ((1 << lev + 1) as f64));
    //Iterate through every chunk
    //num chunks is total size / size_per_chunk = 2^n/2^(lev + 1) = 2^(n-lev-1)
    for chunk in inp.chunks_mut(1 << (lev + 1)) {
        //curr_root_unity = current "twiddle factor"
        let mut cur_root_unity = Complex64::new(1f64, 0f64);
        //iterate through top and bottom halves of chunk in parallel
        //chunk is 2^(lev + 1), half is 2^lev
        for j in 0..(1 << lev) {
            //a1 = jth element of top half, a2 = jth element of bottom half
            let a1 = chunk[j];
            let a2 = chunk[(1 << lev) + j];
            //modify a1 and a2 using twiddle factors and each other
            chunk[j] = a1 + cur_root_unity * a2;
            chunk[(1 << lev) + j] = a1 - cur_root_unity * a2;
            //"increment" twiddle factor
            cur_root_unity *= root_unity;
        }
    }
}

pub fn fft_naive(inp: &Vec<Complex64>) -> Vec<Complex64> {
    assert!(inp.len() & (inp.len() - 1) == 0); // checks if power of 2
    let n = inp.len().ilog2() as usize; // inp.len() = 2^n
    let mut v: Vec<Complex64> = (0..inp.len()).map(|x| inp[bit_reversal(x, n)]).collect();
    for i in 0..n {
        fft2_butterfly(&mut v, i);
    }
    v
}

fn fft2_butterfly_weights(inp: &mut Vec<Complex64>, w: &Vec<Complex64>, lev: usize, n: usize) {
    for (k0, chunk) in inp.chunks_mut(1 << (n - lev)).enumerate() {
        for k2 in 0..(1 << (n - lev - 1)) {
            let a1: Complex64 = chunk[k2];
            let a2: Complex64 = chunk[(1 << (n - lev - 1)) + k2];
            chunk[k2] = a1 + w[bit_reversal(k0, lev) * (1 << (n - lev - 1)) % (1 << n)] * a2;
            chunk[(1 << (n - lev - 1)) + k2] =
                a1 - w[bit_reversal(k0, lev) * (1 << (n - lev - 1)) % (1 << n)] * a2;
        }
    }
}

fn fft4_butterfly_simd(inp: &mut Vec<Complex64>, w: &Vec<Complex64>, lev: usize, n: usize) {
    for (k0, chunk) in inp.chunks_mut(1 << (n - lev)).enumerate() {
        let w1r: f64 = w[1 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        let w1i: f64 = w[1 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im / w1r;
        let w2r: f64 = w[2 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        let w2i: f64 = w[2 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im / w2r;
        let w3r: f64 = w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re / w1r;
        let w3i: f64 = w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im
            / w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        for k2a in 0..(1 << (n - lev - 4)) {
            for k2b in 0..4 {
                let a0 = chunk[k2a * 4 + k2b];
                let a1 = chunk[(1 << (n - lev - 2)) + k2a * 4 + k2b];
                let a2 = chunk[2 * (1 << (n - lev - 2)) + k2a * 4 + k2b];
                let a3 = chunk[3 * (1 << (n - lev - 2)) + k2a * 4 + k2b];
                let b1r = a1.re - a1.im * w1i;
                let b1i = a1.im + a1.re * w1i;
                let b2r = a2.re - a2.im * w2i;
                let b2i = a2.im + a2.re * w2i;
                let b3r = a3.re - a3.im * w3i;
                let b3i = a3.im + a3.re * w3i;
                let c0r = a0.re + b2r * w2r;
                let c0i = a0.im + b2i * w2r;
                let c2r = a0.re - b2r * w2r;
                let c2i = a0.im - b2i * w2r;
                let c1r = b1r + b3r * w3r;
                let c1i = b1i + b3i * w3r;
                let c3r = b1r - b3r * w3r;
                let c3i = b1i - b3i * w3r;
                chunk[k2a * 4 + k2b] = Complex64::new(c0r + c1r * w1r, c0i + c1i * w1r);
                chunk[1 * (1 << (n - lev - 2)) + k2a * 4 + k2b] =
                    Complex64::new(c0r - c1r * w1r, c0i - c1i * w1r);
                chunk[3 * (1 << (n - lev - 2)) + k2a * 4 + k2b] =
                    Complex64::new(c2r - c3i * w1r, c2i + c3r * w1r);
                chunk[2 * (1 << (n - lev - 2)) + k2a * 4 + k2b] =
                    Complex64::new(c2r + c3i * w1r, c2i - c3r * w1r);
            }
        }
    }
}

fn fft4_butterfly_final(inp: &mut Vec<Complex64>, w: &Vec<Complex64>, lev: usize, n: usize) {
    for (k0, chunk) in inp.chunks_mut(1 << (n - lev)).enumerate() {
        let w1r: f64 = w[1 * bit_reversal(k0, lev) % (1 << n)].re;
        let w1i: f64 = w[1 * bit_reversal(k0, lev) % (1 << n)].im / w1r;
        let w2r: f64 = w[2 * bit_reversal(k0, lev) % (1 << n)].re;
        let w2i: f64 = w[2 * bit_reversal(k0, lev) % (1 << n)].im / w2r;
        let w3r: f64 = w[3 * bit_reversal(k0, lev) % (1 << n)].re / w1r;
        let w3i: f64 =
            w[3 * bit_reversal(k0, lev) % (1 << n)].im / w[3 * bit_reversal(k0, lev) % (1 << n)].re;

        let a0 = chunk[0];
        let a1 = chunk[1];
        let a2 = chunk[2];
        let a3 = chunk[3];
        let b1r = a1.re - a1.im * w1i;
        let b1i = a1.im + a1.re * w1i;
        let b2r = a2.re - a2.im * w2i;
        let b2i = a2.im + a2.re * w2i;
        let b3r = a3.re - a3.im * w3i;
        let b3i = a3.im + a3.re * w3i;
        let c0r = a0.re + b2r * w2r;
        let c0i = a0.im + b2i * w2r;
        let c2r = a0.re - b2r * w2r;
        let c2i = a0.im - b2i * w2r;
        let c1r = b1r + b3r * w3r;
        let c1i = b1i + b3i * w3r;
        let c3r = b1r - b3r * w3r;
        let c3i = b1i - b3i * w3r;
        chunk[0] = Complex64::new(c0r + c1r * w1r, c0i + c1i * w1r);
        chunk[1] = Complex64::new(c0r - c1r * w1r, c0i - c1i * w1r);
        chunk[3] = Complex64::new(c2r - c3i * w1r, c2i + c3r * w1r);
        chunk[2] = Complex64::new(c2r + c3i * w1r, c2i - c3r * w1r);
    }
    /*
    let c0 = 1 << (n - lev - 2);
    for (k0, chunk) in inp.chunks_mut(1 << (n - lev)).enumerate() {
        let w1r: f64 = w[1 * bit_reversal(k0, lev) * c0 % (1 << n)].re;
        let w1i: f64 = w[1 * bit_reversal(k0, lev) * c0 % (1 << n)].im / w1r;
        let w2r: f64 = w[2 * bit_reversal(k0, lev) * c0 % (1 << n)].re;
        let w2i: f64 = w[2 * bit_reversal(k0, lev) * c0 % (1 << n)].im / w2r;
        let w3r: f64 = w[3 * bit_reversal(k0, lev) * c0 % (1 << n)].re / w1r;
        let w3i: f64 = w[3 * bit_reversal(k0, lev) * c0 % (1 << n)].im
            / w[3 * bit_reversal(k0, lev) * c0 % (1 << n)].re;

        for k2 in 0..c0 {
            let a0 = chunk[k2];
            let a1 = chunk[c0 + k2];
            let a2 = chunk[2 * c0 + k2];
            let a3 = chunk[3 * c0 + k2];
            let b1r = a1.re - a1.im * w1i;
            let b1i = a1.im + a1.re * w1i;
            let b2r = a2.re - a2.im * w2i;
            let b2i = a2.im + a2.re * w2i;
            let b3r = a3.re - a3.im * w3i;
            let b3i = a3.im + a3.re * w3i;
            let c0r = a0.re + b2r * w2r;
            let c0i = a0.im + b2i * w2r;
            let c2r = a0.re - b2r * w2r;
            let c2i = a0.im - b2i * w2r;
            let c1r = b1r + b3r * w3r;
            let c1i = b1i + b3i * w3r;
            let c3r = b1r - b3r * w3r;
            let c3i = b1i - b3i * w3r;
            chunk[k2] = Complex64::new(c0r + c1r * w1r, c0i + c1i * w1r);
            chunk[1 * c0 + k2] = Complex64::new(c0r - c1r * w1r, c0i - c1i * w1r);
            chunk[3 * c0 + k2] = Complex64::new(c2r - c3i * w1r, c2i + c3r * w1r);
            chunk[2 * c0 + k2] = Complex64::new(c2r + c3i * w1r, c2i - c3r * w1r);
        }
    }
    */
    /*
    for k0 in 0..(1 << lev) {
        let w1r: f64 = w[1 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        let w1i: f64 = w[1 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im / w1r;
        let w2r: f64 = w[2 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        let w2i: f64 = w[2 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im / w2r;
        let w3r: f64 = w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re / w1r;
        let w3i: f64 = w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].im
            / w[3 * bit_reversal(k0, lev) * (1 << (n - lev - 2)) % (1 << n)].re;
        for k2 in 0..(1 << (n - lev - 2)) {
            let a0 = inp[k0 * (1 << (n - lev)) + k2];
            let a1 = inp[k0 * (1 << (n - lev)) + (1 << (n - lev - 2)) + k2];
            let a2 = inp[k0 * (1 << (n - lev)) + 2 * (1 << (n - lev - 2)) + k2];
            let a3 = inp[k0 * (1 << (n - lev)) + 3 * (1 << (n - lev - 2)) + k2];
            let b1r = a1.re - a1.im * w1i;
            let b1i = a1.im + a1.re * w1i;
            let b2r = a2.re - a2.im * w2i;
            let b2i = a2.im + a2.re * w2i;
            let b3r = a3.re - a3.im * w3i;
            let b3i = a3.im + a3.re * w3i;
            let c0r = a0.re + b2r * w2r;
            let c0i = a0.im + b2i * w2r;
            let c2r = a0.re - b2r * w2r;
            let c2i = a0.im - b2i * w2r;
            let c1r = b1r + b3r * w3r;
            let c1i = b1i + b3i * w3r;
            let c3r = b1r - b3r * w3r;
            let c3i = b1i - b3i * w3r;

            inp[k0 * (1 << (n - lev)) + k2] = Complex64::new(c0r + c1r * w1r, c0i + c1i * w1r);
            inp[k0 * (1 << (n - lev)) + 1 * (1 << (n - lev - 2)) + k2] =
                Complex64::new(c0r - c1r * w1r, c0i - c1i * w1r);
            inp[k0 * (1 << (n - lev)) + 3 * (1 << (n - lev - 2)) + k2] =
                Complex64::new(c2r - c3i * w1r, c2i + c3r * w1r);
            inp[k0 * (1 << (n - lev)) + 2 * (1 << (n - lev - 2)) + k2] =
                Complex64::new(c2r + c3i * w1r, c2i - c3r * w1r);

            /*
            inp[k0 * (1 << (n - lev)) + k2].re = c0r + c1r * w1r;
            inp[k0 * (1 << (n - lev)) + k2].im = c0i + c1i * w1r;
            inp[k0 * (1 << (n - lev)) + 1 * (1 << (n - lev - 2)) + k2].re = c0r - c1r * w1r;
            inp[k0 * (1 << (n - lev)) + 1 * (1 << (n - lev - 2)) + k2].im = c0i - c1i * w1r;
            inp[k0 * (1 << (n - lev)) + 3 * (1 << (n - lev - 2)) + k2].re = c2r - c3i * w1r;
            inp[k0 * (1 << (n - lev)) + 3 * (1 << (n - lev - 2)) + k2].im = c2i + c3r * w1r;
            inp[k0 * (1 << (n - lev)) + 2 * (1 << (n - lev - 2)) + k2].re = c2r + c3i * w1r;
            inp[k0 * (1 << (n - lev)) + 2 * (1 << (n - lev - 2)) + k2].im = c2i - c3r * w1r;
            */
        }
    }
    */
}

fn fft4_butterfly0(inp: &mut Vec<Complex64>, n: usize) {
    for k2 in 0..(1 << (n - 2)) {
        let a0 = inp[k2];
        let a1 = inp[(1 << (n - 2)) + k2];
        let a2 = inp[2 * (1 << (n - 2)) + k2];
        let a3 = inp[3 * (1 << (n - 2)) + k2];
        let c0r = a0.re + a2.re;
        let c0i = a0.im + a2.im;
        let c2r = a0.re - a2.re;
        let c2i = a0.im - a2.im;
        let c1r = a1.re + a3.re;
        let c1i = a1.im + a3.im;
        let c3r = a1.re - a3.re;
        let c3i = a1.im - a3.im;
        inp[k2] = Complex64::new(c0r + c1r, c0i + c1i);
        inp[(1 << (n - 2)) + k2] = Complex64::new(c0r - c1r, c0i - c1i);
        inp[3 * (1 << (n - 2)) + k2] = Complex64::new(c2r - c3i, c2i + c3r);
        inp[2 * (1 << (n - 2)) + k2] = Complex64::new(c2r + c3i, c2i - c3r);
    }
}

pub fn fft_optimized(inp: &mut Vec<Complex64>, w: &Vec<Complex64>) {
    let n = inp.len().ilog2() as usize;
    let mut lev = 0;
    if n >= 2 {
        fft4_butterfly0(inp, n);
        lev += 2;
    }
    if n & 1 != 0 {
        fft2_butterfly_weights(inp, w, lev, n);
        lev += 1;
    }
    for i in (lev..(n - 2)).step_by(2) {
        fft4_butterfly_simd(inp, w, i, n);
    }
    if lev != n {
        fft4_butterfly_final(inp, w, n - 2, n);
    }
    for i in 0..(1 << n) {
        if i < bit_reversal(i, n) {
            inp.swap(i, bit_reversal(i, n));
        }
    }
}

pub fn ifft(inp: &Vec<Complex64>) -> Vec<Complex64> {
    let n = inp.len().ilog2() as usize;
    let mut c: Vec<Complex64> = inp.clone();
    (&mut c[1..]).reverse();
    fft_optimized(&mut c, &precompute_weights(n));
    c.iter_mut().for_each(|x| *x /= inp.len() as f64);

    c
}
