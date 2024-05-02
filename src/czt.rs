use crate::fft::{ifft, fft_optimized, precompute_weights};
use num_complex::Complex64;
use std::f64::consts::PI;



pub fn czt(
    inp: Vec<Complex64>,
    m: usize,
    a: Option<Complex64>,
    w: Option<Complex64>,
) -> Vec<Complex64> {
    let n = inp.len();
    let a = a.unwrap_or(Complex64::new(1.0, 0.0));
    let w = w.unwrap_or(Complex64::new(0.0, -2.0 * PI / m as f64).exp());

    let next_pow2 = (n + m - 1).next_power_of_two();

    let mut padded_inp = vec![Complex64::new(0.0, 0.0); next_pow2];
    for i in 0..n {
        //let exp = (Complex64::new(0.0, -PI * i as f64 * i as f64) * w.powi(i as i32)).exp();
        let i_flt = i as f64;
        let exp = w.powf(i_flt * i_flt / 2.0);
        padded_inp[i] = inp[i] / a.powu(i as u32) * exp;
    }

    fft_optimized(&mut padded_inp, &precompute_weights(next_pow2.ilog2() as usize));

    let mut chirp_post = Vec::with_capacity(next_pow2);
    for i in 0..next_pow2 {
        //let exp = (Complex64::new(0.0, -PI * i as f64 * i as f64) * w.powi(-1 * i as i32)).exp();
        if i < m {
            let i_flt = i as f64;
            chirp_post.push(w.powf(i_flt * i_flt / -2.0));
        } else if i < next_pow2 {
            let other_flt = (next_pow2 - i) as f64;
            chirp_post.push(w.powf(other_flt * other_flt / -2.0));
        }
    }

    fft_optimized(&mut chirp_post, &precompute_weights(next_pow2.ilog2() as usize));

    for i in 0..next_pow2 {
        padded_inp[i] = padded_inp[i] * chirp_post[i];
    }

    padded_inp = ifft(&padded_inp);

    let mut out = Vec::with_capacity(m);
    for i in 0..m {
        let i_flt = i as f64;
        let exp = w.powf(i_flt * i_flt / 2.0);
        out.push(padded_inp[i] * exp);
    }

    out
}

pub fn iczt(
    inp: Vec<Complex64>,
    m: usize,
    a: Option<Complex64>,
    w: Option<Complex64>,
) -> Vec<Complex64> {
    let n = inp.len();
    let a = a.unwrap_or(Complex64::new(1.0, 0.0));
    let w = w.unwrap_or(Complex64::new(0.0, -2.0 * PI / m as f64).exp());

    let mut chirp_pm = Vec::with_capacity(n);
    for i in 0..n {
        let exp = (Complex64::new(0.0, PI * i as f64 * i as f64) * w.powi(i as i32)).exp();
        chirp_pm.push(inp[i] * a.powi(i as i32) * exp);
    }

    let next_pow2 = (n + m - 1).next_power_of_two();

    let mut padded_inp = vec![Complex64::new(0.0, 0.0); next_pow2];
    padded_inp[..n].copy_from_slice(&chirp_pm);

    fft_optimized(&mut padded_inp, &precompute_weights(next_pow2));

    let mut chirp_post = Vec::with_capacity(next_pow2);
    for i in 0..next_pow2 {
        let exp = (Complex64::new(0.0, PI * i as f64 * i as f64) * w.powi(-1 * i as i32)).exp();
        chirp_post.push(exp);
    }

    for i in 0..next_pow2 {
        padded_inp[i] = padded_inp[i] * chirp_post[i];
    }

    padded_inp = ifft(&padded_inp);

    let mut out = Vec::with_capacity(m);
    for i in 0..m {
        let exp = (Complex64::new(0.0, PI * i as f64 * i as f64) * w.powi(i as i32)).exp();
        out.push(padded_inp[i] * exp);
    }

    out
}
