use cs128h_project::fft::{fft_optimized, precompute_weights};
use cs128h_project::fileio;
use cs128h_project::czt::{ifft, czt};
use rand::Rng;
use rustfft::{num_complex::Complex64, FftPlanner};
use czt::{transform, c64};
use num_complex::Complex;
use std::f64::consts::PI;

#[test]
fn base2fft_accuracy_randvals() {
    let numvals = 1 << 16;
    let mut rng = rand::thread_rng();
    let vals: Vec<Complex64> = (0..numvals)
        .map(|_| {
            Complex64::new(
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0,
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0,
            )
        })
        .collect();
    // println!("{:?}", vals);
    let mut our_truth = vals.clone();
    let mut their_truth = vals.clone();
    let w = precompute_weights(16);
    fft_optimized(&mut our_truth, &w);
    FftPlanner::new()
        .plan_fft_forward(numvals)
        .process(&mut their_truth);
    // println!("{:?}", our_truth);
    // println!("{:?}", their_truth);
    for i in 0..numvals {
        let diff = our_truth[i] - their_truth[i];
        assert!(diff.norm_sqr() < 0.01);
    }
}

#[test]
fn file_input() {
    let ground_truth1 = vec![
        Complex64::new(1.6, 4.0),
        Complex64::new(7.0, 23.0),
        Complex64::new(41.2, 1.0),
        Complex64::new(194.134, 5.1294),
        Complex64::new(2.127, 5.0),
    ];
    assert_eq!(
        fileio::complex_vec_from_complex_file("./testfiles/test1.txt"),
        ground_truth1
    );
    let ground_truth2 = vec![
        Complex64::new(1.0, 4.0),
        Complex64::new(3.0, 2.0),
        Complex64::new(2.0, 4.0),
        Complex64::new(3.0, 3.0),
        Complex64::new(1.0, 1.0),
        Complex64::new(0.0, 3.0),
        Complex64::new(78.0, 7.0),
        Complex64::new(-182.0, 0.0),
    ];
    assert_eq!(
        fileio::complex_vec_from_complex_file("./testfiles/test2.txt"),
        ground_truth2
    );
}

#[test]
#[should_panic(expected = "couldn't convert these into a decimal value")]
fn file_input_non_decimal() {
    fileio::complex_vec_from_complex_file("./testfiles/invalid.txt");
}

#[test]
#[should_panic(expected = "couldn't open doesntexist")]
fn file_input_bad_path() {
    fileio::complex_vec_from_complex_file("doesntexist");
}

#[test]
fn ifft_accuracy_randvals() {
    let numvals = 1 << 20;
    let mut rng = rand::thread_rng();
    let vals: Vec<Complex64> = (0..numvals)
        .map(|_| {
            Complex64::new(
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0,
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0,
            )
        })
        .collect();
    let mut our_truth = vals.clone();
    let mut their_truth = vals.clone();
    ifft(&mut our_truth);
    FftPlanner::new()
        .plan_fft_forward(numvals)
        .process(&mut their_truth);
    for i in 0..numvals {
        let diff = our_truth[i] - their_truth[i];
        assert!(diff.norm_sqr() < 0.01);
    }
}

#[test]
fn czt_accuracy_randvals() {
    let numvals = 1 << 10;
    let mut rng = rand::thread_rng();

    let mut vals: Vec<Complex64> = (0..numvals)
          .map(|_| {
              Complex64::new(
                  rng.gen::<f64>() * 20.0 - 10.0,
                  rng.gen::<f64>() * 20.0 - 10.0,
              )
          })
          .collect();

    let mut vals2: Vec<c64> = (0..numvals)
          .map(|_| {
            c64::new{
                rng.gen::<f64>() * 20.0 - 10.0,
                rng.gen::<f64>() * 20.0 - 10.0,
            }
          }).collect();
      
      let our_truth = czt(vals, numvals, None, None);
      let w = c64::new(0.0, -2.0 * PI / numvals as f64).exp();
      let a = c64::new(1.0, 0.0);



      vals = transform(&vals2, numvals, w, a);

      for i in 0..numvals {
          let diff = our_truth[i] - vals[i];
          assert!(diff.norm_sqr() < 0.1);
      }
  }
