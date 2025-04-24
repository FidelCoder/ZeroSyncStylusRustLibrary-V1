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
        // Try to get cached constants
        let cache = {
            let cache_map = FIELD_CACHE.read().unwrap();
            cache_map.get(&modulus).cloned()
        };

        let cache = match cache {
            Some(cache) => cache,
            None => {
                // Create new cache entry
                let mut cache_map = FIELD_CACHE.write().unwrap();
                let cache = Arc::new(FieldCache::new(&modulus));
                cache_map.insert(modulus.clone(), cache.clone());
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
        mont_form.reduce();
        BigUint::from_bytes_le(&to_bytes(&mont_form.value))
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
                result.mont_form = result.mont_form.mul(&base.mont_form);
            }
            base.mont_form = base.mont_form.mul(&base.mont_form);
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
    fn add(self, other: Self) -> Self {
        assert_eq!(self.mont_form.constants.modulus, other.mont_form.constants.modulus);
        let modulus_limbs = to_limbs(&self.mont_form.constants.modulus, WORDS_PER_LIMB);
        
        let mut sum = vec![0u64; WORDS_PER_LIMB];
        let mut carry = 0u64;
        
        // Constant-time addition with reduction
        for i in 0..WORDS_PER_LIMB {
            let temp = (self.mont_form.value[i] as u128) + (other.mont_form.value[i] as u128) + (carry as u128);
            sum[i] = temp as u64;
            carry = (temp >> 64) as u64;
        }
        
        // Subtract modulus if result is too large
        if carry > 0 || !ct_lt(&sum, &modulus_limbs) {
            let mut borrow = 0i64;
            for i in 0..WORDS_PER_LIMB {
                let diff = (sum[i] as i128) - (modulus_limbs[i] as i128) - (borrow as i128);
                sum[i] = diff as u64;
                borrow = if diff < 0 { 1 } else { 0 };
            }
        }
        
        Self {
            mont_form: MontgomeryForm::new(sum, self.mont_form.constants.clone()),
        }
    }
}

impl Sub for Fp {
    type Output = Self;

    #[inline(always)]
    fn sub(self, other: Self) -> Self {
        assert_eq!(self.mont_form.constants.modulus, other.mont_form.constants.modulus);
        let modulus_limbs = to_limbs(&self.mont_form.constants.modulus, WORDS_PER_LIMB);
        
        let mut diff = vec![0u64; WORDS_PER_LIMB];
        let mut borrow = 0i64;
        
        // Constant-time subtraction
        for i in 0..WORDS_PER_LIMB {
            let temp = (self.mont_form.value[i] as i128) - (other.mont_form.value[i] as i128) - (borrow as i128);
            diff[i] = temp as u64;
            borrow = if temp < 0 { 1 } else { 0 };
        }
        
        // Add modulus if result is negative
        if borrow != 0 {
            let mut carry = 0u64;
            for i in 0..WORDS_PER_LIMB {
                let temp = (diff[i] as u128) + (modulus_limbs[i] as u128) + (carry as u128);
                diff[i] = temp as u64;
                carry = (temp >> 64) as u64;
            }
        }
        
        Self {
            mont_form: MontgomeryForm::new(diff, self.mont_form.constants.clone()),
        }
    }
}

impl Mul for Fp {
    type Output = Self;

    #[inline(always)]
    fn mul(self, other: Self) -> Self {
        let modulus_limbs = to_limbs(&self.mont_form.constants.modulus, WORDS_PER_LIMB);
        let n_prime_limbs = to_limbs(&self.mont_form.constants.n_prime, WORDS_PER_LIMB);
        
        let result = mont_mul(
            &self.mont_form.value,
            &other.mont_form.value,
            &modulus_limbs,
            &n_prime_limbs
        );
        
        Self {
            mont_form: MontgomeryForm::new(result, self.mont_form.constants.clone()),
        }
    }
}

impl Zero for Fp {
    fn zero() -> Self {
        // Use a small prime for testing
        let modulus = BigUint::from(17u64);
        let cache = {
            let cache_map = FIELD_CACHE.read().unwrap();
            cache_map.get(&modulus).cloned()
        };

        match cache {
            Some(cache) => Self { mont_form: cache.zero.clone() },
            None => Self::new(BigUint::zero(), modulus)
        }
    }

    fn is_zero(&self) -> bool {
        // Check if all limbs are zero after reduction
        let mut mont_form = self.mont_form.clone();
        mont_form.reduce();
        mont_form.value.iter().all(|&x| x == 0)
    }
}

impl One for Fp {
    fn one() -> Self {
        // Use a small prime for testing
        let modulus = BigUint::from(17u64);
        let cache = {
            let cache_map = FIELD_CACHE.read().unwrap();
            cache_map.get(&modulus).cloned()
        };

        match cache {
            Some(cache) => Self { mont_form: cache.one.clone() },
            None => Self::new(BigUint::one(), modulus)
        }
    }
}

impl Neg for Fp {
    type Output = Self;

    fn neg(self) -> Self {
        if self.is_zero() {
            return self;
        }
        let modulus = self.mont_form.constants.modulus.clone();
        // Create a field element representing the modulus
        let modulus_elem = Self::new(modulus.clone(), modulus);
        // Subtract self from modulus to get the negation
        modulus_elem - self
    }
}

impl Div for Fp {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        if let Some(inv) = rhs.inverse() {
            self * inv
        } else {
            panic!("Division by zero")
        }
    }
}

/// Converts a BigUint to a fixed-size array of limbs
#[inline(always)]
fn to_limbs(value: &BigUint, num_limbs: usize) -> Vec<u64> {
    let mut limbs = vec![0u64; num_limbs];
    let bytes = value.to_bytes_le();
    
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
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    bytes
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn test_field_arithmetic() {
        let modulus = BigUint::from(17u64);
        let a = Fp::new(BigUint::from(5u64), modulus.clone());
        let b = Fp::new(BigUint::from(3u64), modulus.clone());
        
        let sum = a.clone() + b.clone();
        let product = a.clone() * b.clone();
        let diff = a - b;
        
        // 5 + 3 = 8 mod 17
        assert_eq!(sum, Fp::new(BigUint::from(8u64), modulus.clone()));
        // 5 * 3 = 15 mod 17
        assert_eq!(product, Fp::new(BigUint::from(15u64), modulus.clone()));
        // 5 - 3 = 2 mod 17
        assert_eq!(diff, Fp::new(BigUint::from(2u64), modulus));
    }

    #[test]
    fn test_montgomery_form() {
        let modulus = BigUint::from(17u32);
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        
        // Test that Montgomery multiplication by 1 returns the same value
        let one = Fp::new(BigUint::from(1u32), modulus);
        let result = a.clone() * one;
        
        assert_eq!(result, a);
    }
} 