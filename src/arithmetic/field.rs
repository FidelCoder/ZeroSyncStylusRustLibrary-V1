use std::ops::{Add, Sub, Mul, Div, Neg};
use std::sync::Arc;
use std::collections::HashMap;
use std::sync::RwLock;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use std::cmp::Ordering;
use num_traits::{Zero, One, ToPrimitive};
use num_integer::Integer;
use crate::arithmetic::{
    traits::{Field, PrimeField},
    montgomery::{MontgomeryConstants, MontgomeryForm, mont_mul, ct_lt, to_limbs, to_bytes},
};
use std::str::FromStr;

/// Word size for field operations
const WORD_SIZE: u32 = 64;
const WORDS_PER_LIMB: usize = 4;

// Cache for commonly used field elements and constants
lazy_static! {
    static ref FIELD_CACHE: RwLock<HashMap<BigUint, Arc<FieldCache>>> = RwLock::new(HashMap::new());
}

struct FieldCache {
    // Common constants in Montgomery form
    zero: MontgomeryForm,
    one: MontgomeryForm,
    // Pre-computed small values
    small_values: Vec<MontgomeryForm>,
    // Montgomery constants
    constants: MontgomeryConstants,
}

impl FieldCache {
    fn new(modulus: &BigUint) -> Self {
        let constants = MontgomeryConstants::new(modulus, 64);
        let zero = MontgomeryForm::new(vec![0; 4], constants.clone());
        let one = {
            let mut value = vec![0; 4];
            value[0] = 1;
            MontgomeryForm::new(value, constants.clone())
        };

        // Pre-compute small values up to 16
        let mut small_values = Vec::with_capacity(16);
        let mut current = one.clone();
        for _ in 0..16 {
            small_values.push(current.clone());
            current = current.mul(&one);
        }

        Self {
            zero,
            one,
            small_values,
            constants,
        }
    }

    fn get_small_value(&self, value: u64) -> Option<MontgomeryForm> {
        if value < self.small_values.len() as u64 {
            Some(self.small_values[value as usize].clone())
        } else {
            None
        }
    }
}

/// Represents an element of a prime field using Montgomery arithmetic
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Fp {
    /// The value in Montgomery form
    mont_form: MontgomeryForm,
}

impl PartialOrd for Fp {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // First, make sure we're working with the same modulus
        if self.mont_form.constants.modulus != other.mont_form.constants.modulus {
            return None;
        }
        
        // Convert both values to standard form for comparison
        let a_val = self.from_montgomery();
        let b_val = other.from_montgomery();
        
        a_val.partial_cmp(&b_val)
    }
}

impl Fp {
    /// Creates a new field element
    pub fn new(value: BigUint, modulus: BigUint) -> Self {
        // Try to get cached constants - handle RwLock errors gracefully
        let cache = match FIELD_CACHE.read() {
            Ok(cache_map) => {
                cache_map.get(&modulus).cloned()
            },
            Err(_) => {
                // If the lock is poisoned, we'll create a new cache instance
                None
            },
        };

        let cache = match cache {
            Some(cache) => cache,
            None => {
                // Create new cache entry - handle RwLock errors gracefully
                let cache = Arc::new(FieldCache::new(&modulus));
                
                // Try to update the cache, but continue even if it fails
                let _ = FIELD_CACHE.write().map(|mut cache_map| {
                    cache_map.insert(modulus.clone(), cache.clone());
                });
                
                cache
            }
        };

        // Check if it's a small value first
        if let Some(small_value) = value.to_u64().and_then(|v| cache.get_small_value(v)) {
            return Self { mont_form: small_value };
        }

        // Convert to Montgomery form
        let reduced_value = value % &cache.constants.modulus;
        let mut mont_form = MontgomeryForm::new(
            to_limbs(&reduced_value, 4),
            cache.constants.clone(),
        );

        // Convert to Montgomery form by multiplying by R^2 mod N
        let r_squared_limbs = to_limbs(&cache.constants.r_squared, 4);
        let modulus_limbs = to_limbs(&cache.constants.modulus, 4);
        let n_prime_limbs = to_limbs(&cache.constants.n_prime, 4);
        
        mont_form.value = mont_mul(
            &mont_form.value,
            &r_squared_limbs,
            &modulus_limbs,
            &n_prime_limbs
        );

        Self { mont_form }
    }

    /// Converts the value from Montgomery form
    pub fn from_montgomery(&self) -> BigUint {
        let mut mont_form = self.mont_form.clone();
        
        // Make sure the value is fully reduced before converting
        mont_form.reduce();
        
        // Convert from Montgomery form by multiplying by 1
        let one_limbs = vec![1u64; 4];
        let modulus_limbs = to_limbs(&mont_form.constants.modulus, 4);
        let n_prime_limbs = to_limbs(&mont_form.constants.n_prime, 4);
        
        let result_limbs = mont_mul(
            &mont_form.value,
            &one_limbs,
            &modulus_limbs,
            &n_prime_limbs
        );
        
        // Convert to BigUint safely
        BigUint::from_bytes_le(&to_bytes(&result_limbs))
    }

    pub fn square(&mut self) -> Self {
        let mut result = self.clone();
        let base_mont = self.mont_form.clone();
        result.mont_form = base_mont.clone().mul(&base_mont);
        result
    }
    
    /// Get the modulus of this field element
    pub fn modulus(&self) -> BigUint {
        self.mont_form.constants.modulus.clone()
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Convert from Montgomery form
        let mut a = self.mont_form.value.clone();
        let modulus_limbs = to_limbs(&self.mont_form.constants.modulus, WORDS_PER_LIMB);
        let one_limbs = vec![1u64; WORDS_PER_LIMB];
        let n_prime_limbs = to_limbs(&self.mont_form.constants.n_prime, WORDS_PER_LIMB);

        // Convert from Montgomery form
        a = mont_mul(&a, &one_limbs, &modulus_limbs, &n_prime_limbs);

        // Convert to BigUint for inverse calculation
        let a_biguint = BigUint::from_bytes_le(&to_bytes(&a));
        let modulus = self.mont_form.constants.modulus.clone();

        // Implementation of modular inverse using extended Euclidean algorithm
        let (g, mut x) = mod_inverse(&a_biguint, &modulus);
        
        // Check if gcd is 1 (meaning inverse exists)
        if g != BigUint::one() {
            return None;
        }
        
        // Ensure x is reduced modulo the modulus
        x = x % &modulus;

        // Convert back to Montgomery form
        let mut result = Self::new(x, modulus);
        result.mont_form.reduce();
        Some(result)
    }
}

impl Field for Fp {
    fn characteristic() -> Vec<u64> {
        // TODO: Implement proper characteristic calculation
        vec![1]
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        // Convert from Montgomery form
        let mut a = self.mont_form.value.clone();
        let modulus_limbs = to_limbs(&self.mont_form.constants.modulus, WORDS_PER_LIMB);
        let one_limbs = vec![1u64; WORDS_PER_LIMB];
        let n_prime_limbs = to_limbs(&self.mont_form.constants.n_prime, WORDS_PER_LIMB);

        // Convert from Montgomery form
        a = mont_mul(&a, &one_limbs, &modulus_limbs, &n_prime_limbs);

        // Convert to BigUint for inverse calculation
        let a_biguint = BigUint::from_bytes_le(&to_bytes(&a));
        let modulus = self.mont_form.constants.modulus.clone();

        // Implementation of modular inverse using extended Euclidean algorithm
        let (g, mut x) = mod_inverse(&a_biguint, &modulus);
        
        // Check if gcd is 1 (meaning inverse exists)
        if g != BigUint::one() {
            return None;
        }
        
        // Ensure x is reduced modulo the modulus
        x = x % &modulus;

        // Convert back to Montgomery form
        let mut result = Self::new(x, modulus);
        result.mont_form.reduce();
        Some(result)
    }

    fn pow(&self, exp: u64) -> Self {
        let mut result = Self::one();
        let mut base = self.clone();
        let mut exp = exp;

        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base.clone();
            }
            base = base.square();
            exp /= 2;
        }

        result
    }
}

impl PrimeField for Fp {
    fn modulus() -> Vec<u64> {
        // TODO: Return proper modulus representation
        vec![1]
    }

    fn to_montgomery(&self) -> Self {
        self.clone()
    }

    fn from_montgomery(&self) -> Self {
        let value = self.from_montgomery();
        let modulus = self.modulus();
        Self::new(value, modulus)
    }
}

impl Add for Fp {
    type Output = Self;

    fn add(mut self, other: Self) -> Self {
        self.mont_form = self.mont_form.add(&other.mont_form);
        self
    }
}

impl Sub for Fp {
    type Output = Self;

    fn sub(mut self, other: Self) -> Self {
        self.mont_form = self.mont_form.sub(&other.mont_form);
        self
    }
}

impl Mul for Fp {
    type Output = Self;

    fn mul(mut self, other: Self) -> Self {
        self.mont_form = self.mont_form.mul(&other.mont_form);
        self
    }
}

impl Zero for Fp {
    fn zero() -> Self {
        let modulus = BigUint::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583"
        ).unwrap();
        Self::new(BigUint::zero(), modulus)
    }

    fn is_zero(&self) -> bool {
        self.mont_form.value.iter().all(|&x| x == 0)
    }
}

impl One for Fp {
    fn one() -> Self {
        let modulus = BigUint::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583"
        ).unwrap();
        Self::new(BigUint::one(), modulus)
    }
}

impl Neg for Fp {
    type Output = Self;

    fn neg(mut self) -> Self {
        let mut zero = Self::zero();
        self.mont_form = zero.mont_form.sub(&self.mont_form);
        self
    }
}

impl Div for Fp {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        let inv = rhs.inverse().expect("Division by zero");
        self * inv
    }
}

// Helper function for modular inverse using extended Euclidean algorithm for BigUint
fn mod_inverse(a: &BigUint, m: &BigUint) -> (BigUint, BigUint) {
    let mut s = BigUint::zero();
    let mut old_s = BigUint::one();
    let mut t = BigUint::one();
    let mut old_t = BigUint::zero();
    let mut r = m.clone();
    let mut old_r = a.clone();

    while !r.is_zero() {
        let quotient = &old_r / &r;
        
        // Update old_r and r
        let temp_r = r.clone();
        r = modular_sub(&old_r, &(quotient.clone() * &r), m);
        old_r = temp_r;
        
        // Update old_s and s
        let temp_s = s.clone();
        let qs = quotient.clone() * &s;
        s = if old_s >= qs {
            old_s.clone() - qs
        } else {
            m - (qs - old_s.clone()) % m
        };
        old_s = temp_s;
        
        // Update old_t and t
        let temp_t = t.clone();
        let qt = quotient * &t;
        t = if old_t >= qt {
            old_t.clone() - qt
        } else {
            m - (qt - old_t.clone()) % m
        };
        old_t = temp_t;
    }

    (old_r, old_s)
}

// Safe modular subtraction that handles the case where b > a
fn modular_sub(a: &BigUint, b: &BigUint, m: &BigUint) -> BigUint {
    if a >= b {
        (a - b) % m
    } else {
        (m - (b - a) % m) % m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;
    
    #[test]
    fn test_field_creation() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        
        // Create field elements
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        
        // Convert to standard form
        assert_eq!(a.from_montgomery(), BigUint::from(5u32));
    }
    
    #[test]
    fn test_field_addition() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        
        // Create field elements
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        let b = Fp::new(BigUint::from(7u32), modulus.clone());
        
        // Addition that doesn't exceed modulus
        let sum = a.clone() + b.clone();
        assert_eq!(sum.from_montgomery(), BigUint::from(12u32));
        
        // Addition that exceeds modulus
        let c = Fp::new(BigUint::from(10u32), modulus.clone());
        let d = Fp::new(BigUint::from(9u32), modulus.clone());
        let sum2 = c + d;
        assert_eq!(sum2.from_montgomery(), BigUint::from(2u32)); // (10 + 9) mod 17 = 19 mod 17 = 2
    }
    
    #[test]
    fn test_field_subtraction() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        
        // Create field elements
        let a = Fp::new(BigUint::from(12u32), modulus.clone());
        let b = Fp::new(BigUint::from(7u32), modulus.clone());
        
        // Subtraction where a > b
        let diff = a.clone() - b.clone();
        assert_eq!(diff.from_montgomery(), BigUint::from(5u32));
        
        // Subtraction where a < b
        let c = Fp::new(BigUint::from(5u32), modulus.clone());
        let diff2 = c - b;
        assert_eq!(diff2.from_montgomery(), BigUint::from(15u32)); // (5 - 7) mod 17 = -2 mod 17 = 15
    }
    
    #[test]
    fn test_field_multiplication() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        
        // Create field elements
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        let b = Fp::new(BigUint::from(7u32), modulus.clone());
        
        // Multiplication
        let prod = a * b;
        assert_eq!(prod.from_montgomery(), BigUint::from(1u32)); // (5 * 7) mod 17 = 35 mod 17 = 1
    }
    
    #[test]
    fn test_field_inverse() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        
        // Create field element
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        
        // Compute inverse
        let a_inv = a.inverse().unwrap();
        
        // Verify a * a^-1 = 1
        let prod = a.clone() * a_inv;
        assert_eq!(prod.from_montgomery(), BigUint::from(1u32));
    }
} 