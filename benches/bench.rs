use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cs128h_project::fft::{fft_naive, fft_optimized, precompute_weights};
use rand::Rng;
use rustfft::{num_complex::Complex64, FftPlanner};

fn fft_benchmark_compare(c: &mut Criterion) {
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
    let planner = FftPlanner::new().plan_fft_forward(numvals);
    let weights = precompute_weights(10);
    let mut group = c.benchmark_group("FFT");
    group.bench_function(BenchmarkId::new("OursNaive", "ComplexVec"), |b| {
        b.iter(|| {
            fft_naive(&vals);
        })
    });
    group.bench_function(BenchmarkId::new("OursOptimized", "ComplexVec"), |b| {
        b.iter(|| {
            let mut vals_copy = vals.clone();
            fft_optimized(&mut vals_copy, &weights);
        })
    });
    group.bench_function(BenchmarkId::new("RustFFT", "ComplexVec"), |b| {
        b.iter(|| {
            let mut vals_copy = vals.clone();
            planner.process(&mut vals_copy);
        })
    });
}

fn fft_optimized_bench(c: &mut Criterion) {
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
    let weights = precompute_weights(10);
    c.bench_function("Our FFT Optimized", |b| {
        b.iter(|| {
            let mut vals_copy = vals.clone();
            fft_optimized(&mut vals_copy, &weights);
        })
    });
}

criterion_group! {
    name = compare;
    config = Criterion::default();
    targets = fft_benchmark_compare
}

criterion_group! {
    name = optimize;
    config = Criterion::default();
    targets = fft_optimized_bench
}

criterion_main!(optimize);
