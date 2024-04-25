use cs128h_project::fft::fft;
use cs128h_project::fileio;
use num_complex::Complex64;

fn main() {
    let v: Vec<Complex64> = vec![
        Complex64::new(2f64, 0f64),
        Complex64::new(1f64, 0f64),
        Complex64::new(-1f64, 0f64),
        Complex64::new(5f64, 0f64),
        Complex64::new(0f64, 0f64),
        Complex64::new(3f64, 0f64),
        Complex64::new(0f64, 0f64),
        Complex64::new(-4f64, 0f64),
    ];
    let h = fft(&v);
    println!("{:?}", h);
}
