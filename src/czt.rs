use num_complex::Complex64;
use crate::fft::{fft_optimized, precompute_weights};
use std::f64::consts::PI;

pub fn ifft(inp: &Vec<Complex64>) -> Vec<Complex64> {
    let mut n = inp.len();
    let mut c: Vec<Complex64> = inp.iter().map(|&x| x.conj()).collect(); 
    fft_optimized(&mut c, &precompute_weights(n)); 
    c.iter_mut().for_each(|x| *x /= n as f64);

    c
}

pub fn czt(inp: Vec<Complex64>, m: usize, a: Option<Complex64>, w: Option<Complex64>) -> Vec<Complex64> {
    let n = inp.len();
    let a = a.unwrap_or(Complex64::new(1.0, 0.0)); 
    let w = w.unwrap_or(Complex64::new(0.0, -2.0 * PI / m as f64).exp()); 

    let mut chirp_pm = Vec::with_capacity(n);
    for i in 0..n {
        let exp = (Complex64::new(0.0, -PI * i as f64 * i as f64) * w.powi(i as i32)).exp();
        chirp_pm.push(inp[i] * a.powi(i as i32) * exp);
    }

      let next_pow2 = (n + m - 1).next_power_of_two();

      let mut padded_inp = vec![Complex64::new(0.0, 0.0); next_pow2];
      padded_inp[..n].copy_from_slice(&chirp_pm);

      fft_optimized(&mut padded_inp, &precompute_weights(next_pow2));

      let mut chirp_post = Vec::with_capacity(next_pow2);
      for i in 0..next_pow2 {
          let exp = (Complex64::new(0.0, -PI * i as f64 * i as f64) * w.powi(-1 * i as i32)).exp();
          chirp_post.push(exp);
      }

      for i in 0..next_pow2 {
          padded_inp[i] = padded_inp[i] * chirp_post[i];
      }

      padded_inp = ifft(&mut padded_inp);

      let mut out = Vec::with_capacity(m);
      for i in 0..m {
          let exp = (Complex64::new(0.0, -PI * i as f64 * i as f64) * w.powi(i as i32)).exp();
          out.push(padded_inp[i] * exp);
      }

      out
  }

pub fn iczt(inp: Vec<Complex64>, m: usize, a: Option<Complex64>, w: Option<Complex64>) -> Vec<Complex64> {
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