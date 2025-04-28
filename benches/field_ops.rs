use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use num_bigint::BigUint;
use zerosync::arithmetic::{
    field::Fp,
    simd::{field_mul_avx2, field_add_avx2, has_avx2},
};
use rand::Rng;
use zerosync::arithmetic::traits::Field;
use num_traits::{Zero, One};
use std::str::FromStr;
use std::time::Duration;

/// BN254 prime field modulus
const BN254_MODULUS: &str = "21888242871839275222246405745257275088696311157297823662689037894645226208583";

/// Generate a random Fp element
fn random_fp() -> Fp {
    let mut rng = rand::thread_rng();
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    
    // Generate random value less than modulus
    let limbs: [u64; 4] = [
        rng.gen(),
        rng.gen(),
        rng.gen(),
        rng.gen_range(0..0x0FFFFFFFFFFFFFFF),
    ];
    
    let bytes: Vec<u8> = limbs.iter().flat_map(|x| x.to_le_bytes()).collect();
    let value = BigUint::from_bytes_le(&bytes) % &modulus;
    
    Fp::new(value, modulus)
}

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

fn bench_field_add(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let a = Fp::new(BigUint::from(12345u32), modulus.clone());
    let b = Fp::new(BigUint::from(67890u32), modulus);
    
    let mut group = c.benchmark_group("Field Addition");
    group.measurement_time(Duration::from_secs(5));
    
    group.bench_function("simple_addition", |bench| {
        bench.iter(|| {
            black_box(a.clone() + b.clone())
        })
    });
    
    // Benchmark addition with different input sizes
    let test_cases = vec![
        ("small_values", (10u32, 20u32)), 
        ("medium_values", (u32::MAX, u32::MAX)),
        ("large_values", (12345678u64, 87654321u64)),
    ];
    
    for (name, (v1, v2)) in test_cases {
        let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
        let a = Fp::new(BigUint::from(v1), modulus.clone());
        let b = Fp::new(BigUint::from(v2), modulus);
        
        group.bench_with_input(BenchmarkId::new("addition", name), &(a, b), |bench, (a, b)| {
            bench.iter(|| {
                black_box(a.clone() + b.clone())
            })
        });
    }
    
    group.finish();
}

fn bench_field_mul(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let a = Fp::new(BigUint::from(12345u32), modulus.clone());
    let b = Fp::new(BigUint::from(67890u32), modulus.clone());
    
    let mut group = c.benchmark_group("Field Multiplication");
    group.measurement_time(Duration::from_secs(5));
    
    group.bench_function("simple_multiplication", |bench| {
        bench.iter(|| {
            black_box(a.clone() * b.clone())
        })
    });
    
    // Benchmark multiplication with different input sizes
    let test_cases = vec![
        ("small_values", (10u32, 20u32)), 
        ("medium_values", (u32::MAX, u32::MAX)),
        ("large_values", (12345678u64, 87654321u64)),
    ];
    
    for (name, (v1, v2)) in test_cases {
        let a = Fp::new(BigUint::from(v1), modulus.clone());
        let b = Fp::new(BigUint::from(v2), modulus.clone());
        
        group.bench_with_input(BenchmarkId::new("multiplication", name), &(a, b), |bench, (a, b)| {
            bench.iter(|| {
                black_box(a.clone() * b.clone())
            })
        });
    }
    
    // Benchmark squaring (should be faster than general multiplication)
    let values = vec![
        ("small", 42u32),
        ("medium", u32::MAX),
        ("large", 12345678901u64),
    ];
    
    for (name, val) in values {
        let a = Fp::new(BigUint::from(val), modulus.clone());
        
        group.bench_with_input(BenchmarkId::new("squaring", name), &a, |bench, a| {
            bench.iter(|| {
                black_box(a.clone() * a.clone())
            })
        });
    }
    
    group.finish();
}

fn bench_field_inverse(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let mut group = c.benchmark_group("Field Inversion");
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(50); // Inversion is expensive, use fewer samples
    
    // Benchmark inversion with different input sizes
    let values = vec![
        ("small", 42u32),
        ("medium", u32::MAX),
        ("large", 12345678901u64),
        ("random", 0u32), // placeholder for random
    ];
    
    for (name, val) in values {
        let a = if name == "random" {
            random_fp()
        } else {
            Fp::new(BigUint::from(val), modulus.clone())
        };
        
        group.bench_with_input(BenchmarkId::new("inversion", name), &a, |bench, a| {
            bench.iter(|| {
                black_box(a.inverse().unwrap())
            })
        });
    }
    
    group.finish();
}

fn bench_field_exponentiation(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let mut group = c.benchmark_group("Field Exponentiation");
    group.measurement_time(Duration::from_secs(5));
    group.sample_size(50);
    
    // Benchmark exponentation with different exponents
    let base = Fp::new(BigUint::from(7u32), modulus.clone());
    let exponents = vec![
        ("tiny", 3u64),
        ("small", 42u64),
        ("medium", 65537u64),
        ("large", u32::MAX as u64),
        ("xlarge", (1u64 << 32) + 1),
    ];
    
    for (name, exp) in exponents {
        group.bench_with_input(BenchmarkId::new("exponentiation", name), &(base.clone(), exp), 
            |bench, (base, exp)| {
                bench.iter(|| {
                    black_box(base.pow(*exp))
                })
            }
        );
    }
    
    group.finish();
}

fn bench_montgomery_operations(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let mut group = c.benchmark_group("Montgomery Form");
    group.measurement_time(Duration::from_secs(5));
    
    // Benchmark Montgomery form conversion
    let values = vec![
        ("small", 42u32),
        ("medium", u32::MAX),
        ("large", 12345678901u64),
    ];
    
    for (name, val) in values {
        let a = Fp::new(BigUint::from(val), modulus.clone());
        
        // To Montgomery form
        group.bench_with_input(BenchmarkId::new("to_montgomery", name), &a, |bench, a| {
            bench.iter(|| {
                // Force Montgomery conversion (implementation detail)
                let val = a.clone();
                black_box(val)
            })
        });
        
        // From Montgomery form
        group.bench_with_input(BenchmarkId::new("from_montgomery", name), &a, |bench, a| {
            bench.iter(|| {
                black_box(a.from_montgomery())
            })
        });
    }
    
    group.finish();
}

fn bench_lazy_reduction(c: &mut Criterion) {
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let mut group = c.benchmark_group("Lazy Reduction");
    group.measurement_time(Duration::from_secs(5));
    
    // Generate random values
    let values: Vec<Fp> = (0..10).map(|_| random_fp()).collect();
    
    // Benchmark chain of operations with lazy reduction
    let operations = vec![
        ("10_multiplications", 10),
        ("20_multiplications", 20),
        ("50_multiplications", 50),
    ];
    
    for (name, num_ops) in operations {
        group.bench_function(name, |bench| {
            bench.iter(|| {
                let mut result = values[0].clone();
                for i in 1..num_ops {
                    // This should leverage lazy reduction internally
                    result = result * values[i % values.len()].clone();
                }
                black_box(result)
            })
        });
    }
    
    group.finish();
}

fn bench_gas_estimation(c: &mut Criterion) {
    // This benchmark simulates Arbitrum Stylus gas usage
    let modulus = BigUint::from_str(BN254_MODULUS).unwrap();
    let mut group = c.benchmark_group("Stylus Gas Estimation");
    group.measurement_time(Duration::from_secs(5));
    
    // Create field elements
    let a = random_fp();
    let b = random_fp();
    
    // Benchmark common operations for gas estimation
    group.bench_function("add_gas", |bench| {
        bench.iter(|| {
            black_box(a.clone() + b.clone())
        })
    });
    
    group.bench_function("mul_gas", |bench| {
        bench.iter(|| {
            black_box(a.clone() * b.clone())
        })
    });
    
    group.bench_function("inverse_gas", |bench| {
        bench.iter(|| {
            black_box(a.clone().inverse().unwrap())
        })
    });
    
    // Simulate a common ZK operation pattern: multiple multiplications followed by addition
    group.bench_function("zk_pattern_gas", |bench| {
        bench.iter(|| {
            let t1 = a.clone() * b.clone();
            let t2 = t1.clone() * a.clone();
            let t3 = t2.clone() * b.clone();
            black_box(t3 + a.clone())
        })
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_field_add,
    bench_field_mul,
    bench_field_inverse,
    bench_field_exponentiation,
    bench_montgomery_operations,
    bench_lazy_reduction,
    bench_gas_estimation,
);
criterion_main!(benches); 