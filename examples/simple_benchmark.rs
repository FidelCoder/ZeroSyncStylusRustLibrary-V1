use num_bigint::BigUint;
use std::str::FromStr;
use std::time::{Duration, Instant};

/// A simplified field element implementation
#[derive(Clone, Debug, PartialEq)]
struct Field {
    value: BigUint,
    modulus: BigUint,
}

impl Field {
    fn new(value: BigUint, modulus: BigUint) -> Self {
        let reduced = value % &modulus;
        Self {
            value: reduced,
            modulus,
        }
    }
    
    fn add(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let sum = (&self.value + &other.value) % &self.modulus;
        Self {
            value: sum,
            modulus: self.modulus.clone(),
        }
    }
    
    fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let product = (&self.value * &other.value) % &self.modulus;
        Self {
            value: product,
            modulus: self.modulus.clone(),
        }
    }
    
    fn square(&self) -> Self {
        let squared = (&self.value * &self.value) % &self.modulus;
        Self {
            value: squared,
            modulus: self.modulus.clone(),
        }
    }
}

fn benchmark<F>(name: &str, iterations: u32, f: F) -> Duration
where
    F: Fn(),
{
    // Warm up
    for _ in 0..10 {
        f();
    }
    
    // Actual benchmark
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    let duration = start.elapsed();
    
    let ns_per_op = duration.as_nanos() / iterations as u128;
    println!("{}: {} iterations, {} ns/op", name, iterations, ns_per_op);
    
    duration
}

fn main() {
    println!("ZeroSync Simple Benchmark");
    println!("=========================");
    
    // BN254 field modulus
    let bn254_modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    
    // Create some field elements with various sizes
    let small_a = Field::new(BigUint::from(123u32), bn254_modulus.clone());
    let small_b = Field::new(BigUint::from(456u32), bn254_modulus.clone());
    
    let medium_a = Field::new(BigUint::from_str("123456789012345").unwrap(), bn254_modulus.clone());
    let medium_b = Field::new(BigUint::from_str("987654321098765").unwrap(), bn254_modulus.clone());
    
    let large_a = Field::new(BigUint::from_str("12345678901234567890123456789").unwrap(), bn254_modulus.clone());
    let large_b = Field::new(BigUint::from_str("98765432109876543210987654321").unwrap(), bn254_modulus.clone());
    
    // Benchmark addition with small values
    benchmark("Small Addition", 100_000, || {
        let _ = small_a.add(&small_b);
    });
    
    // Benchmark addition with medium values
    benchmark("Medium Addition", 100_000, || {
        let _ = medium_a.add(&medium_b);
    });
    
    // Benchmark addition with large values
    benchmark("Large Addition", 100_000, || {
        let _ = large_a.add(&large_b);
    });
    
    // Benchmark multiplication with small values
    benchmark("Small Multiplication", 10_000, || {
        let _ = small_a.mul(&small_b);
    });
    
    // Benchmark multiplication with medium values
    benchmark("Medium Multiplication", 10_000, || {
        let _ = medium_a.mul(&medium_b);
    });
    
    // Benchmark multiplication with large values
    benchmark("Large Multiplication", 10_000, || {
        let _ = large_a.mul(&large_b);
    });
    
    // Benchmark squaring
    benchmark("Small Squaring", 10_000, || {
        let _ = small_a.square();
    });
    
    benchmark("Medium Squaring", 10_000, || {
        let _ = medium_a.square();
    });
    
    benchmark("Large Squaring", 10_000, || {
        let _ = large_a.square();
    });
    
    // Benchmark a sequence of operations
    benchmark("Mixed Operations", 10_000, || {
        let t1 = small_a.add(&medium_a);
        let t2 = t1.mul(&large_a);
        let _ = t2.square();
    });
} 