use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cs128h_project::fft::{fft, fft_in_place};
use rand::Rng;
use rustfft::{num_complex::Complex64, FftPlanner};

fn base2fft_benchmark(c: &mut Criterion) {
    let numvals = 1 << 10;
    let mut rng = rand::thread_rng();
    let vals: Vec<Complex64> = (0..numvals)
        .map(|_| {
            Complex64::new(
                rng.gen::<f64>() * 200000.0 - 100000.0,
                rng.gen::<f64>() * 200000.0 - 100000.0,
            )
        })
        .collect();
    let mut group = c.benchmark_group("FFT");
    group.bench_function(BenchmarkId::new("Ours", "ComplexVec"), |b| {
        b.iter(|| {
            let our = vals.clone();
            fft(&our)
        })
    });
    group.bench_function(BenchmarkId::new("OursInPlace", "ComplexVec"), |b| {
        b.iter(|| {
            let mut our = vals.clone();
            fft_in_place(&mut our)
        })
    });
    group.bench_function(BenchmarkId::new("RustFFT", "ComplexVec"), |b| {
        b.iter(|| {
            let mut their = vals.clone();
            FftPlanner::new()
                .plan_fft_forward(numvals)
                .process(&mut their)
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = base2fft_benchmark
}
criterion_main!(benches);
