use num_bigint::BigUint;
use num_traits::Num;
use std::str::FromStr;
use zerosync::arithmetic::field::Fp;

fn main() {
    println!("ZeroSync Field Arithmetic Example");
    println!("================================");
    
    // Create a small field with modulus 101 (small prime for demonstration)
    let modulus = BigUint::from(101u32);
    println!("Field modulus: {}", modulus);
    
    // Create some field elements
    let a = Fp::new(BigUint::from(30u32), modulus.clone());
    let b = Fp::new(BigUint::from(50u32), modulus.clone());
    println!("a = 30, b = 50");
    
    // Addition
    let sum = a.clone() + b.clone();
    println!("a + b = {}", sum.from_montgomery());
    
    // Multiplication
    let product = a.clone() * b.clone();
    println!("a * b = {} (expected {})", product.from_montgomery(), (BigUint::from(30u32) * BigUint::from(50u32)) % &modulus);
    
    // Try with BN254 curve parameters
    println!("\nBN254 Curve Demonstration");
    println!("=========================");
    
    // BN254 field modulus
    let bn254_modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    println!("BN254 field modulus: {}", bn254_modulus);
    
    // Create field elements
    let c = Fp::new(BigUint::from(12345u32), bn254_modulus.clone());
    let d = Fp::new(BigUint::from(67890u32), bn254_modulus.clone());
    println!("c = 12345, d = 67890");
    
    // Addition
    let sum2 = c.clone() + d.clone();
    println!("c + d = {}", sum2.from_montgomery());
    
    // Multiplication
    let product2 = c.clone() * d.clone();
    let expected = (BigUint::from(12345u32) * BigUint::from(67890u32)) % &bn254_modulus;
    println!("c * d = {}", product2.from_montgomery());
    println!("Expected: {}", expected);
} 