use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zerosync::arithmetic::field::Fp;
use zerosync::curve::bn254::{BN254, G1Affine, G2Affine};
use num_bigint::BigUint;
use std::str::FromStr;

fn field_operations_benchmark(c: &mut Criterion) {
    let modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    
    let a = Fp::new(BigUint::from(5u32), modulus.clone());
    let b = Fp::new(BigUint::from(3u32), modulus.clone());
    
    let mut group = c.benchmark_group("Field Operations");
    
    // Field addition
    group.bench_function("Field Addition", |bench| {
        bench.iter(|| {
            black_box(a.clone() + b.clone());
        });
    });
    
    // Field multiplication
    group.bench_function("Field Multiplication", |bench| {
        bench.iter(|| {
            black_box(a.clone() * b.clone());
        });
    });
    
    // Field inversion
    group.bench_function("Field Inversion", |bench| {
        bench.iter(|| {
            black_box(a.clone().inverse());
        });
    });
    
    group.finish();
}

fn curve_operations_benchmark(c: &mut Criterion) {
    let curve = BN254::new();
    let g1 = curve.g1_generator();
    let g2 = curve.g2_generator();
    
    let mut group = c.benchmark_group("Curve Operations");
    
    // G1 point addition
    group.bench_function("G1 Point Addition", |bench| {
        bench.iter(|| {
            black_box(g1.clone() + g1.clone());
        });
    });
    
    // G1 scalar multiplication
    group.bench_function("G1 Scalar Multiplication", |bench| {
        bench.iter(|| {
            black_box(g1.clone() * 2);
        });
    });
    
    // G1 windowed scalar multiplication
    group.bench_function("G1 Windowed Scalar Multiplication", |bench| {
        bench.iter(|| {
            black_box(g1.windowed_mul(&BigUint::from(2u32)));
        });
    });
    
    // G2 point addition
    group.bench_function("G2 Point Addition", |bench| {
        bench.iter(|| {
            black_box(g2.clone() + g2.clone());
        });
    });
    
    // G2 scalar multiplication
    group.bench_function("G2 Scalar Multiplication", |bench| {
        bench.iter(|| {
            black_box(g2.clone() * 2);
        });
    });
    
    // G2 windowed scalar multiplication
    group.bench_function("G2 Windowed Scalar Multiplication", |bench| {
        bench.iter(|| {
            black_box(g2.windowed_mul(&BigUint::from(2u32)));
        });
    });
    
    group.finish();
}

fn gas_analysis_benchmark(c: &mut Criterion) {
    let modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    
    let curve = BN254::new();
    let g1 = curve.g1_generator();
    let g2 = curve.g2_generator();
    
    let mut group = c.benchmark_group("Gas Analysis");
    
    // Field operations gas cost
    group.bench_function("Field Addition Gas", |bench| {
        bench.iter(|| {
            let a = Fp::new(BigUint::from(5u32), modulus.clone());
            let b = Fp::new(BigUint::from(3u32), modulus.clone());
            black_box(a + b);
        });
    });
    
    group.bench_function("Field Multiplication Gas", |bench| {
        bench.iter(|| {
            let a = Fp::new(BigUint::from(5u32), modulus.clone());
            let b = Fp::new(BigUint::from(3u32), modulus.clone());
            black_box(a * b);
        });
    });
    
    // Curve operations gas cost
    group.bench_function("G1 Point Addition Gas", |bench| {
        bench.iter(|| {
            black_box(g1.clone() + g1.clone());
        });
    });
    
    group.bench_function("G1 Scalar Multiplication Gas", |bench| {
        bench.iter(|| {
            black_box(g1.clone() * 2);
        });
    });
    
    group.bench_function("G2 Point Addition Gas", |bench| {
        bench.iter(|| {
            black_box(g2.clone() + g2.clone());
        });
    });
    
    group.bench_function("G2 Scalar Multiplication Gas", |bench| {
        bench.iter(|| {
            black_box(g2.clone() * 2);
        });
    });
    
    group.finish();
}

criterion_group!(
    benches,
    field_operations_benchmark,
    curve_operations_benchmark,
    gas_analysis_benchmark
);
criterion_main!(benches); 