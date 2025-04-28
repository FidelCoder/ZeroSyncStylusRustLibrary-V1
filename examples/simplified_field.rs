use num_bigint::BigUint;
use std::str::FromStr;
use std::ops::{Add, Sub, Mul, Div, Neg};

/// A simplified field element implementation for ZeroSync demonstration
#[derive(Clone, Debug, PartialEq)]
struct SimpleField {
    value: BigUint,
    modulus: BigUint,
}

impl SimpleField {
    fn new(value: BigUint, modulus: BigUint) -> Self {
        // Make sure value is reduced modulo the modulus
        let reduced_value = value % &modulus;
        Self {
            value: reduced_value,
            modulus,
        }
    }
    
    fn to_string(&self) -> String {
        self.value.to_string()
    }
    
    fn pow(&self, exp: u32) -> Self {
        let result = self.value.modpow(&BigUint::from(exp), &self.modulus);
        Self::new(result, self.modulus.clone())
    }
    
    fn inverse(&self) -> Option<Self> {
        // Extended Euclidean algorithm to find modular inverse
        if self.value == BigUint::from(0u32) {
            return None;
        }
        
        let mut t = BigUint::from(0u32);
        let mut new_t = BigUint::from(1u32);
        let mut r = self.modulus.clone();
        let mut new_r = self.value.clone();
        
        while new_r != BigUint::from(0u32) {
            let quotient = &r / &new_r;
            
            // Update t
            let tmp_t = t.clone();
            t = new_t.clone();
            // t = new_t and new_t = tmp_t - quotient * new_t
            if &quotient * &new_t <= tmp_t {
                new_t = tmp_t - &quotient * &new_t;
            } else {
                // Handle underflow in modular arithmetic
                new_t = &self.modulus - ((&quotient * &new_t - tmp_t) % &self.modulus);
            }
            
            // Update r
            let tmp_r = r.clone();
            r = new_r.clone();
            new_r = tmp_r - &quotient * &new_r;
        }
        
        if r > BigUint::from(1u32) {
            return None; // Not invertible
        }
        
        Some(Self::new(t, self.modulus.clone()))
    }
}

impl Add for SimpleField {
    type Output = Self;
    
    fn add(self, other: Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let sum = (&self.value + &other.value) % &self.modulus;
        Self {
            value: sum,
            modulus: self.modulus,
        }
    }
}

impl Sub for SimpleField {
    type Output = Self;
    
    fn sub(self, other: Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let diff = if self.value >= other.value {
            &self.value - &other.value
        } else {
            // Handle underflow: result = modulus - (other - self)
            &self.modulus - ((&other.value - &self.value) % &self.modulus)
        };
        
        Self {
            value: diff,
            modulus: self.modulus,
        }
    }
}

impl Mul for SimpleField {
    type Output = Self;
    
    fn mul(self, other: Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let product = (&self.value * &other.value) % &self.modulus;
        Self {
            value: product,
            modulus: self.modulus,
        }
    }
}

impl Neg for SimpleField {
    type Output = Self;
    
    fn neg(self) -> Self {
        if self.value == BigUint::from(0u32) {
            return self;
        }
        
        let neg_value = &self.modulus - &self.value;
        Self {
            value: neg_value,
            modulus: self.modulus,
        }
    }
}

impl Div for SimpleField {
    type Output = Self;
    
    fn div(self, other: Self) -> Self {
        assert_eq!(self.modulus, other.modulus, "Field moduli must match");
        
        let inv = other.inverse().expect("Division by zero or non-invertible element");
        self * inv
    }
}

fn main() {
    println!("ZeroSync Simplified Field Arithmetic Example");
    println!("===========================================");
    
    // Create a small field with modulus 101 (small prime for demonstration)
    let modulus = BigUint::from(101u32);
    println!("Field modulus: {}", modulus);
    
    // Create some field elements
    let a = SimpleField::new(BigUint::from(30u32), modulus.clone());
    let b = SimpleField::new(BigUint::from(50u32), modulus.clone());
    println!("a = {}, b = {}", a.to_string(), b.to_string());
    
    // Addition
    let sum = a.clone() + b.clone();
    println!("a + b = {}", sum.to_string());
    
    // Subtraction
    let diff = a.clone() - b.clone();
    println!("a - b = {}", diff.to_string());
    
    // Multiplication
    let product = a.clone() * b.clone();
    println!("a * b = {}", product.to_string());
    
    // Negation
    let neg_a = -a.clone();
    println!("-a = {}", neg_a.to_string());
    
    // Inverse
    if let Some(inv_a) = a.clone().inverse() {
        println!("a^(-1) = {}", inv_a.to_string());
    } else {
        println!("a has no inverse");
    }
    
    // Division
    let quotient = a.clone() / b.clone();
    println!("a / b = {}", quotient.to_string());
    
    // Power
    let power = a.pow(3);
    println!("a^3 = {}", power.to_string());
    
    // Try with BN254 curve parameters
    println!("\nBN254 Curve Demonstration");
    println!("=========================");
    
    // BN254 field modulus
    let bn254_modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    println!("BN254 field modulus: {}", bn254_modulus);
    
    // Create field elements
    let c = SimpleField::new(BigUint::from(12345u32), bn254_modulus.clone());
    let d = SimpleField::new(BigUint::from(67890u32), bn254_modulus.clone());
    println!("c = {}, d = {}", c.to_string(), d.to_string());
    
    // Addition
    let sum2 = c.clone() + d.clone();
    println!("c + d = {}", sum2.to_string());
    
    // Multiplication
    let product2 = c.clone() * d.clone();
    println!("c * d = {}", product2.to_string());
    
    // Verify with direct computation
    let expected = (BigUint::from(12345u32) * BigUint::from(67890u32)) % &bn254_modulus;
    println!("Expected c * d = {}", expected.to_string());
    assert_eq!(product2.value, expected);
} 