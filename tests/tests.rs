use cs128h_project::fft::fft;
use cs128h_project::fileio;
use rand::Rng;
use rustfft::{num_complex::Complex64, FftPlanner};

#[test]
fn base2fft_accuracy_randvals() {
    let numvals = 1 << 25;
    let mut rng = rand::thread_rng();
    let mut vals: Vec<Complex64> = (0..numvals)
        .map(|_| {
            Complex64::new(
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0,
                ((rng.gen::<f64>() * 200.0 - 100.0) * 1000.0).trunc() / 1000.0
            )
        })
        .collect();
    let our_truth = fft(&vals);
    FftPlanner::new()
        .plan_fft_forward(numvals)
        .process(&mut vals);
    for i in 0..numvals {
        let diff = our_truth[i] - vals[i];
        assert!(diff.norm_sqr() < 0.1);
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
        fileio::complex_vec_from_file("testfiles\\test1.txt"),
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
        fileio::complex_vec_from_file("testfiles\\test2.txt"),
        ground_truth2
    );
}

#[test]
#[should_panic(expected = "couldn't convert these into a decimal value")]
fn file_input_non_decimal() {
    fileio::complex_vec_from_file("testfiles\\invalid.txt");
}

#[test]
#[should_panic(expected = "couldn't open doesntexist")]
fn file_input_bad_path() {
    fileio::complex_vec_from_file("doesntexist");
}
