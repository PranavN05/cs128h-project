

extern crate rand;
use core::num;

use rand::Rng;
use rustfft::{num_complex::Complex32, FftPlanner};

#[test]
fn test_accuracy_random_values() {
    let numvals = (1 << 14);
    let mut rng = rand::thread_rng();
    let mut vals: Vec<Complex32> = (0..numvals)
        .map(|_| {
            Complex32::new(
                rng.gen::<f32>() * 20.0 - 10.0,
                rng.gen::<f32>() * 20.0 - 10.0,
            )
        })
        .collect();
    let our_truth = crate::fft(&vals);
    FftPlanner::new().plan_fft_forward(numvals).process(&mut vals);
    for i in 0..numvals {
        let diff = our_truth[i] - vals[i];
        assert!(diff.norm_sqr() < 0.1);
    }
}
