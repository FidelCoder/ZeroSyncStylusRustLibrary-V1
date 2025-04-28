use std::ops::{Add, Sub, Mul, Div, Neg};
use std::sync::Arc;
use std::collections::HashMap;
use std::sync::RwLock;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use num_traits::{Zero, One, ToPrimitive};
use num_integer::Integer;
use crate::arithmetic::{
    traits::{Field, PrimeField},
    montgomery::{MontgomeryConstants, MontgomeryForm, mont_mul, ct_lt},
};

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
        mont_form.reduce();

        Self { mont_form }
    }

    /// Converts the value from Montgomery form
    pub fn from_montgomery(&self) -> BigUint {
        let mut mont_form = self.mont_form.clone();
        
        // Make sure the value is fully reduced before converting
        mont_form.reduce();
        
        // Convert to BigUint safely
        let mut bytes = Vec::with_capacity(mont_form.value.len() * 8);
        for &limb in &mont_form.value {
            for j in 0..8 {
                bytes.push(((limb >> (j * 8)) & 0xFF) as u8);
            }
        }
        
        BigUint::from_bytes_le(&bytes)
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

        // Extended Binary GCD
        let mut u = a_biguint;
        let mut v = modulus.clone();
        let mut b = BigUint::one();
        let mut c = BigUint::zero();

        while !u.is_zero() {
            while u.is_even() {
                u >>= 1;
                if b.is_even() {
                    b >>= 1;
                } else {
                    b = (b.clone() + modulus.clone()) >> 1;
                }
            }

            while v.is_even() {
                v >>= 1;
                if c.is_even() {
                    c >>= 1;
                } else {
                    c = (c.clone() + modulus.clone()) >> 1;
                }
            }

            if u >= v {
                u = u - v.clone();
                if b >= c.clone() {
                    b = b - c.clone();
                } else {
                    b = modulus.clone() - (c.clone() - b.clone());
                }
            } else {
                v = v - u.clone();
                if c.clone() >= b.clone() {
                    c = c - b.clone();
                } else {
                    c = modulus.clone() - (b.clone() - c.clone());
                }
            }
        }

        if v != BigUint::one() {
            return None;
        }

        Some(Self::new(c, modulus))
    }

    fn pow(&self, exp: u64) -> Self {
        let mut base = self.clone();
        let mut result = Self::new(BigUint::one(), self.mont_form.constants.modulus.clone());
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                let base_mont = base.mont_form.clone();
                result.mont_form = result.mont_form.mul(&base_mont);
            }
            let base_mont = base.mont_form.clone();
            base.mont_form = base_mont.clone().mul(&base_mont);
            e >>= 1;
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
        // TODO: Implement conversion from Montgomery form
        self.clone()
    }
}

impl Add for Fp {
    type Output = Self;

    #[inline(always)]
    fn add(mut self, other: Self) -> Self {
        assert_eq!(self.mont_form.constants.modulus, other.mont_form.constants.modulus);
        
        // Use the built-in add method from MontgomeryForm
        let result = self.mont_form.add(&other.mont_form);
        
        Self {
            mont_form: result,
        }
    }
}

impl Sub for Fp {
    type Output = Self;

    #[inline(always)]
    fn sub(mut self, other: Self) -> Self {
        assert_eq!(self.mont_form.constants.modulus, other.mont_form.constants.modulus);
        
        // Use the built-in sub method from MontgomeryForm
        let result = self.mont_form.sub(&other.mont_form);
        
        Self {
            mont_form: result,
        }
    }
}

impl Mul for Fp {
    type Output = Self;

    #[inline(always)]
    fn mul(self, other: Self) -> Self {
        assert_eq!(self.mont_form.constants.modulus, other.mont_form.constants.modulus);
        
        let mut a = self.mont_form.clone();
        let result = a.mul(&other.mont_form);
        
        Self {
            mont_form: result,
        }
    }
}

impl Zero for Fp {
    fn zero() -> Self {
        // This is a very simplified implementation for demo purposes
        // In a real implementation, we would need to ensure the modulus is correct
        let modulus = BigUint::from(101u32); // Example small prime
        Self::new(BigUint::zero(), modulus)
    }

    fn is_zero(&self) -> bool {
        let mut mont_form = self.mont_form.clone();
        mont_form.reduce();
        mont_form.value.iter().all(|&x| x == 0)
    }
}

impl One for Fp {
    fn one() -> Self {
        // This is a very simplified implementation for demo purposes
        // In a real implementation, we would need to ensure the modulus is correct
        let modulus = BigUint::from(101u32); // Example small prime
        Self::new(BigUint::one(), modulus)
    }
}

impl Neg for Fp {
    type Output = Self;
    
    fn neg(self) -> Self {
        if self.is_zero() {
            return self;
        }
        
        // Get the modulus and implement proper modular negation
        let modulus = self.modulus();
        
        // Create a separate representation to avoid directly using from_montgomery
        // which could lead to additional errors
        let mut result_limbs = vec![0u64; WORDS_PER_LIMB];
        let modulus_limbs = to_limbs(&modulus, WORDS_PER_LIMB);
        let value_limbs = self.mont_form.value.clone();
        
        // If value is non-zero, compute modulus - value
        if !value_limbs.iter().all(|&x| x == 0) {
            let mut borrow = 0i64;
            
            // Subtract: modulus - value
            for i in 0..WORDS_PER_LIMB {
                let diff = (modulus_limbs[i] as i128) - (value_limbs[i] as i128) - (borrow as i128);
                result_limbs[i] = diff as u64;
                borrow = if diff < 0 { 1 } else { 0 };
            }
        }
        
        Self {
            mont_form: MontgomeryForm::new(result_limbs, self.mont_form.constants.clone()),
        }
    }
}

impl Div for Fp {
    type Output = Self;
    
    fn div(self, rhs: Self) -> Self {
        match rhs.inverse() {
            Some(inv) => self * inv,
            None => panic!("Division by zero"),
        }
    }
}

/// Converts a BigUint to a fixed-size array of limbs
#[inline(always)]
fn to_limbs(value: &BigUint, num_limbs: usize) -> Vec<u64> {
    let bytes = value.to_bytes_le();
    let mut limbs = vec![0u64; num_limbs];
    
    // Convert bytes to limbs
    for (i, chunk) in bytes.chunks(8).enumerate() {
        if i >= num_limbs {
            break;
        }
        
        let mut limb = 0u64;
        for (j, &byte) in chunk.iter().enumerate() {
            limb |= (byte as u64) << (j * 8);
        }
        limbs[i] = limb;
    }
    
    limbs
}

/// Converts a slice of u64 limbs to bytes in little-endian order
fn to_bytes(limbs: &[u64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(limbs.len() * 8);
    for &limb in limbs {
        for j in 0..8 {
            bytes.push(((limb >> (j * 8)) & 0xFF) as u8);
        }
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;
    
    #[test]
    fn test_field_arithmetic() {
        // Test with a larger prime modulus to avoid underflow issues
        let modulus = BigUint::from(101u32);
        
        // Create field elements with values less than the modulus
        let a = Fp::new(BigUint::from(30u32), modulus.clone());
        let b = Fp::new(BigUint::from(50u32), modulus.clone());
        
        // Addition that doesn't exceed modulus
        let sum = a.clone() + b.clone();
        let expected_sum = Fp::new(BigUint::from(80u32), modulus.clone());
        assert_eq!(sum.from_montgomery(), expected_sum.from_montgomery());
        
        // Subtraction that needs modular reduction
        let diff = a.clone() - b.clone();
        // 30 - 50 = -20, which in modular arithmetic is 101 - 20 = 81
        let expected_diff = Fp::new(BigUint::from(81u32), modulus.clone());
        assert_eq!(diff.from_montgomery(), expected_diff.from_montgomery());
        
        // Multiplication
        let prod = a.clone() * b.clone();
        // 30 * 50 = 1500, which is 1500 % 101 = 86
        let expected_prod = Fp::new(BigUint::from(86u32), modulus.clone());
        assert_eq!(prod.from_montgomery(), expected_prod.from_montgomery());
    }
    
    #[test]
    fn test_montgomery_form() {
        // Use a small prime for predictable results
        let modulus = BigUint::from(101u32);
        
        let a = Fp::new(BigUint::from(30u32), modulus.clone());
        let b = Fp::new(BigUint::from(40u32), modulus.clone());
        
        // Basic operations
        let c = a.clone() + b.clone();
        let d = a.clone() * b.clone();
        
        // 30 + 40 = 70
        assert_eq!(c.from_montgomery() % &modulus, BigUint::from(70u32));
        
        // 30 * 40 = 1200, which is 1200 % 101 = 88
        assert_eq!(d.from_montgomery() % &modulus, BigUint::from(88u32));
    }
} 