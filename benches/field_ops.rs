use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use num_bigint::BigUint;
use zerosync::arithmetic::{
    field::Fp,
    simd::{field_mul_avx2, field_add_avx2, has_avx2},
};
use rand::Rng;

fn field_arithmetic_benchmark(c: &mut Criterion) {
    // BN254 base field modulus
    let modulus = BigUint::from_str_radix(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583",
        10
    ).unwrap();

    let mut group = c.benchmark_group("Field Operations");
    
    // Standard arithmetic benchmarks
    {
        let a = Fp::new(BigUint::from(12345u32), modulus.clone());
        let b = Fp::new(BigUint::from(67890u32), modulus.clone());

        group.bench_function("standard/addition", |bencher| {
            bencher.iter(|| {
                black_box(a.clone() + b.clone())
            })
        });

        group.bench_function("standard/multiplication", |bencher| {
            bencher.iter(|| {
                black_box(a.clone() * b.clone())
            })
        });

        group.bench_function("standard/squaring", |bencher| {
            bencher.iter(|| {
                black_box(a.clone() * a.clone())
            })
        });
    }

    // SIMD arithmetic benchmarks
    if has_avx2() {
        let mut rng = rand::thread_rng();
        let a: [u64; 4] = std::array::from_fn(|_| rng.gen());
        let b: [u64; 4] = std::array::from_fn(|_| rng.gen());
        let modulus: [u64; 4] = [
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0x0FFFFFFFFFFFFFFF,
        ];

        group.bench_function("simd/addition", |bencher| {
            bencher.iter(|| unsafe {
                black_box(field_add_avx2(&a, &b, &modulus))
            })
        });

        group.bench_function("simd/multiplication", |bencher| {
            bencher.iter(|| unsafe {
                black_box(field_mul_avx2(&a, &b, &modulus))
            })
        });
    }

    // Batch operation benchmarks
    {
        let size = 1000;
        let mut rng = rand::thread_rng();
        let elements: Vec<Fp> = (0..size)
            .map(|_| Fp::new(BigUint::from(rng.gen::<u32>()), modulus.clone()))
            .collect();

        group.bench_function("batch/sum", |bencher| {
            bencher.iter(|| {
                elements.iter().fold(
                    Fp::new(BigUint::from(0u32), modulus.clone()),
                    |acc, x| acc + x.clone()
                )
            })
        });

        group.bench_function("batch/product", |bencher| {
            bencher.iter(|| {
                elements.iter().fold(
                    Fp::new(BigUint::from(1u32), modulus.clone()),
                    |acc, x| acc * x.clone()
                )
            })
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .sample_size(100)
        .measurement_time(std::time::Duration::from_secs(5));
    targets = field_arithmetic_benchmark
}
criterion_main!(benches); 