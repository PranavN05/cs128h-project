extern crate rand;
use rand::Rng;
use rustfft::{num_complex::Complex32, FftPlanner};

#[test]
fn test_accuracy_random_values() {
    let mut rng = rand::thread_rng();
    let mut vals: Vec<Complex32> = (0..(1 << 16))
        .map(|_| {
            Complex32::new(
                rng.gen::<f32>() * 20.0 - 10.0,
                rng.gen::<f32>() * 20.0 - 10.0,
            )
        })
        .collect();
    let our_truth = crate::fft(&vals);
    FftPlanner::new().plan_fft_forward((1 << 16)).process(&mut vals);
    for i in 0..(1 << 16) {
        let diff = our_truth[i] - vals[i];
        assert!(diff.norm_sqr() < 0.1);
    }
}
