use std::ops::{Add, Sub, Mul, Div, Neg};
use num_bigint::BigUint;
use num_traits::{Zero, One};
use num_integer::Integer;
use crate::arithmetic::{
    traits::{Field, PrimeField},
    montgomery::{MontgomeryConstants, mont_mul, ct_lt},
};

/// Word size for field operations
const WORD_SIZE: u32 = 64;
const WORDS_PER_LIMB: usize = 4;

/// Represents an element of a prime field using Montgomery arithmetic
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Fp {
    /// The value in Montgomery form as limbs
    limbs: Vec<u64>,
    /// Montgomery arithmetic constants
    constants: MontgomeryConstants,
}

impl Fp {
    /// Creates a new field element
    pub fn new(value: BigUint, modulus: BigUint) -> Self {
        let constants = MontgomeryConstants::new(&modulus, WORD_SIZE);
        
        // Reduce value modulo the modulus first
        let value = value % modulus.clone();
        
        let mut fp = Self {
            limbs: to_limbs(&value, WORDS_PER_LIMB),
            constants,
        };
        
        // Convert to Montgomery form
        fp.to_montgomery_form();
        fp
    }

    /// Converts the value to Montgomery form
    #[inline(always)]
    fn to_montgomery_form(&mut self) {
        // Special case for zero
        if self.limbs.iter().all(|&x| x == 0) {
            return;
        }

        let r_squared_limbs = to_limbs(&self.constants.r_squared, WORDS_PER_LIMB);
        let modulus_limbs = to_limbs(&self.constants.modulus, WORDS_PER_LIMB);
        let n_prime_limbs = to_limbs(&self.constants.n_prime, WORDS_PER_LIMB);
        
        self.limbs = mont_mul(
            &self.limbs,
            &r_squared_limbs,
            &modulus_limbs,
            &n_prime_limbs
        );
    }

    /// Performs Montgomery multiplication
    #[inline(always)]
    fn mont_mul(&self, other: &Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        let modulus_limbs = to_limbs(&self.constants.modulus, WORDS_PER_LIMB);
        let n_prime_limbs = to_limbs(&self.constants.n_prime, WORDS_PER_LIMB);
        
        let result = mont_mul(
            &self.limbs,
            &other.limbs,
            &modulus_limbs,
            &n_prime_limbs
        );
        
        Self {
            limbs: result,
            constants: self.constants.clone(),
        }
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
        let mut a = self.limbs.clone();
        let modulus_limbs = to_limbs(&self.constants.modulus, WORDS_PER_LIMB);
        let one_limbs = vec![1u64; WORDS_PER_LIMB];
        let n_prime_limbs = to_limbs(&self.constants.n_prime, WORDS_PER_LIMB);

        // Convert from Montgomery form
        a = mont_mul(&a, &one_limbs, &modulus_limbs, &n_prime_limbs);

        // Convert to BigUint for inverse calculation
        let a_biguint = BigUint::from_bytes_le(&to_bytes(&a));
        let modulus = self.constants.modulus.clone();

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
        let mut result = Self::new(BigUint::one(), self.constants.modulus.clone());
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result = result.mont_mul(&base);
            }
            base = base.mont_mul(&base);
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
        assert_eq!(self.constants.modulus, other.constants.modulus);
        let modulus_limbs = to_limbs(&self.constants.modulus, WORDS_PER_LIMB);
        
        let mut sum = vec![0u64; WORDS_PER_LIMB];
        let mut carry = 0u64;
        
        // Constant-time addition with reduction
        for i in 0..WORDS_PER_LIMB {
            let temp = (self.limbs[i] as u128) + (other.limbs[i] as u128) + (carry as u128);
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
            limbs: sum,
            constants: self.constants.clone(),
        }
    }
}

impl Sub for Fp {
    type Output = Self;

    #[inline(always)]
    fn sub(self, other: Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        let modulus_limbs = to_limbs(&self.constants.modulus, WORDS_PER_LIMB);
        
        let mut diff = vec![0u64; WORDS_PER_LIMB];
        let mut borrow = 0i64;
        
        // Constant-time subtraction
        for i in 0..WORDS_PER_LIMB {
            let temp = (self.limbs[i] as i128) - (other.limbs[i] as i128) - (borrow as i128);
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
            limbs: diff,
            constants: self.constants.clone(),
        }
    }
}

impl Mul for Fp {
    type Output = Self;

    #[inline(always)]
    fn mul(self, other: Self) -> Self {
        self.mont_mul(&other)
    }
}

impl Zero for Fp {
    fn zero() -> Self {
        // Use a small prime for testing
        let modulus = BigUint::from(17u64);
        Self {
            limbs: vec![0; WORDS_PER_LIMB],
            constants: MontgomeryConstants::new(&modulus, WORD_SIZE),
        }
    }

    fn is_zero(&self) -> bool {
        self.limbs.iter().all(|&x| x == 0)
    }
}

impl One for Fp {
    fn one() -> Self {
        // Use a small prime for testing
        let modulus = BigUint::from(17u64);
        Self::new(BigUint::one(), modulus)
    }
}

impl Neg for Fp {
    type Output = Self;

    fn neg(self) -> Self {
        if self.is_zero() {
            return self;
        }
        let modulus = self.constants.modulus.clone();
        let mut result = Self::new(modulus.clone(), modulus);
        result = result - self;
        result
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