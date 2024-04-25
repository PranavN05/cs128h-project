use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use cs128h_project::fft::fft;
use rand::Rng;
use rustfft::{num_complex::Complex64, FftPlanner};

fn base2fft_benchmark(c: &mut Criterion) {
    let numvals = 1 << 23;
    let mut rng = rand::thread_rng();
    let mut vals: Vec<Complex64> = (0..numvals)
        .map(|_| {
            Complex64::new(
                rng.gen::<f64>() * 200000.0 - 100000.0,
                rng.gen::<f64>() * 200000.0 - 100000.0,
            )
        })
        .collect();
    let mut group = c.benchmark_group("FFT");
    group.bench_function(BenchmarkId::new("Ours", "ComplexVec"), |b| {
        b.iter(|| fft(&vals))
    });
    group.bench_function(BenchmarkId::new("RustFFT", "ComplexVec"), |b| {
        b.iter(|| {
            FftPlanner::new()
                .plan_fft_forward(numvals)
                .process(&mut vals)
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = base2fft_benchmark
}
criterion_main!(benches);
