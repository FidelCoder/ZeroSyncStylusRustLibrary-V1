//! Prime field implementation using Montgomery representation
//! Optimized for Stylus environment

use crate::field::{Field, FieldError, FieldResult};
use num_bigint::{BigInt, BigUint, ToBigInt};
use num_traits::{One, Zero};
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::{Rng, thread_rng};

/// Montgomery representation of a field element
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp {
    /// The internal Montgomery representation value
    value: BigUint,
    /// The field parameters
    params: &'static FpParameters,
}

/// Parameters for the prime field
#[derive(Debug)]
pub struct FpParameters {
    /// The modulus of the field
    pub modulus: BigUint,
    /// The number of bits in the modulus
    pub bits: usize,
    /// R = 2^bits mod modulus (Montgomery constant)
    pub r: BigUint,
    /// R^2 mod modulus (used for Montgomery conversion)
    pub r_squared: BigUint,
    /// Inverse of modulus mod 2^64 (used for Montgomery reduction)
    pub inv: u64,
}

impl FpParameters {
    /// Create new field parameters
    pub fn new(modulus: BigUint) -> FieldResult<Self> {
        // Ensure modulus is odd (requirement for Montgomery arithmetic)
        if (&modulus & BigUint::from(1u64)) == BigUint::zero() {
            return Err(FieldError::NonPrimeModulus);
        }
        
        let bits = modulus.bits() as usize;
        let r = (BigUint::from(1u64) << bits) % &modulus;
        let r_squared = (&r * &r) % &modulus;
        
        // Compute -m^-1 mod 2^64 (Montgomery inverse)
        let mut inv = 1u64;
        for _ in 0..63 {
            inv = inv.wrapping_mul(inv);
            inv = inv.wrapping_mul(modulus.to_u64_digits()[0]);
        }
        inv = inv.wrapping_neg();
        
        Ok(Self {
            modulus,
            bits,
            r,
            r_squared,
            inv,
        })
    }
    
    /// Convert a value to Montgomery form
    fn to_montgomery(&self, value: &BigUint) -> BigUint {
        (value * &self.r_squared) % &self.modulus
    }
    
    /// Convert a value from Montgomery form
    fn from_montgomery(&self, value: &BigUint) -> BigUint {
        self.montgomery_reduce(value, &BigUint::from(1u64))
    }
    
    /// Montgomery reduction algorithm
    fn montgomery_reduce(&self, a: &BigUint, b: &BigUint) -> BigUint {
        let product = a * b;
        let u = (&product * BigUint::from(self.inv)) & ((BigUint::from(1u64) << 64) - BigUint::from(1u64));
        let m_u = &self.modulus * u;
        let result = (product + m_u) >> 64;
        
        if result >= self.modulus {
            result - &self.modulus
        } else {
            result
        }
    }
}

/// Precomputed parameters for BN254 scalar field
pub static BN254_FR_PARAMS: FpParameters = {
    // r = 21888242871839275222246405745257275088548364400416034343698204186575808495617
    // Modulus of the BN254 scalar field
    let modulus = BigUint::from_bytes_be(&[
        0x30, 0x64, 0x4e, 0x72, 0xe1, 0x31, 0xa0, 0x29, 0xf2, 0x2b, 0x4d, 0x5b, 
        0x7b, 0x4b, 0x78, 0x78, 0x64, 0xc1, 0x3f, 0x4d, 0x7e, 0xc0, 0x59, 0x29, 
        0x9e, 0x1d, 0b1, 0x21, 0x85, 0xf9, 0x11, 0x3b
    ]);
    
    // Calculate parameters
    let bits = modulus.bits() as usize;
    let r = (BigUint::from(1u64) << bits) % &modulus;
    let r_squared = (&r * &r) % &modulus;
    
    // Precomputed -m^-1 mod 2^64
    let inv = 0x87d20782e4866389;
    
    FpParameters {
        modulus,
        bits,
        r,
        r_squared,
        inv,
    }
};

impl Fp {
    /// Create a new field element from a BigUint value
    pub fn new(value: BigUint, params: &'static FpParameters) -> Self {
        let mont_value = params.to_montgomery(&(value % &params.modulus));
        Self {
            value: mont_value,
            params,
        }
    }
    
    /// Create a new field element from a u64 value
    pub fn from_u64(value: u64, params: &'static FpParameters) -> Self {
        Self::new(BigUint::from(value), params)
    }
    
    /// Get the value in standard (non-Montgomery) form
    pub fn get_value(&self) -> BigUint {
        self.params.from_montgomery(&self.value)
    }
    
    /// Create a field element for the BN254 scalar field
    pub fn bn254_scalar(value: impl Into<BigUint>) -> Self {
        Self::new(value.into(), &BN254_FR_PARAMS)
    }
}

impl Field for Fp {
    fn zero() -> Self {
        Self {
            value: BigUint::zero(),
            params: &BN254_FR_PARAMS, // Default to BN254
        }
    }
    
    fn one() -> Self {
        Self {
            value: BN254_FR_PARAMS.r.clone(), // R in Montgomery form is 1
            params: &BN254_FR_PARAMS,
        }
    }
    
    fn random() -> Self {
        let mut rng = thread_rng();
        let bits = BN254_FR_PARAMS.bits;
        let bytes = (bits + 7) / 8;
        
        let mut buf = vec![0u8; bytes];
        rng.fill(&mut buf[..]);
        
        // Ensure value is less than modulus
        let value = BigUint::from_bytes_be(&buf) % &BN254_FR_PARAMS.modulus;
        Self::new(value, &BN254_FR_PARAMS)
    }
    
    fn is_zero(&self) -> bool {
        self.value.is_zero()
    }
    
    fn is_one(&self) -> bool {
        self.value == self.params.r
    }
    
    fn inverse(&self) -> FieldResult<Self> {
        if self.is_zero() {
            return Err(FieldError::DivisionByZero);
        }
        
        // Extended Euclidean algorithm for modular inverse
        let mut a = self.get_value().to_bigint().unwrap();
        let mut b = self.params.modulus.to_bigint().unwrap();
        let mut x = BigInt::one();
        let mut y = BigInt::zero();
        
        while !b.is_zero() {
            let q = &a / &b;
            let r = &a % &b;
            let nx = &y;
            let ny = &x - &q * &y;
            
            a = b;
            b = r;
            x = nx.clone();
            y = ny;
        }
        
        if a != BigInt::one() {
            return Err(FieldError::OperationError("Inverse does not exist".into()));
        }
        
        let result = if x < BigInt::zero() {
            x + self.params.modulus.to_bigint().unwrap()
        } else {
            x
        };
        
        let value = result.to_biguint().unwrap();
        Ok(Self::new(value, self.params))
    }
    
    fn sqrt(&self) -> FieldResult<Self> {
        // Tonelli-Shanks algorithm for square root
        // Simplified version for fields where p ≡ 3 (mod 4)
        // For BN254, we use the property that p ≡ 3 (mod 4)
        
        if self.is_zero() {
            return Ok(*self);
        }
        
        let modulus = &self.params.modulus;
        let exp = (modulus + BigUint::from(1u64)) / BigUint::from(4u64);
        
        // a^((p+1)/4) mod p
        let value = self.get_value().modpow(&exp, modulus);
        let result = Self::new(value, self.params);
        
        // Verify result
        if result.square() != *self {
            return Err(FieldError::OperationError("No square root exists".into()));
        }
        
        Ok(result)
    }
}

// Implement arithmetic operations using Montgomery representation
impl Add for Fp {
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self::Output {
        assert!(std::ptr::eq(self.params, rhs.params), "Field parameters must match");
        
        let mut result = self.value + rhs.value;
        if result >= self.params.modulus {
            result -= &self.params.modulus;
        }
        
        Self {
            value: result,
            params: self.params,
        }
    }
}

impl AddAssign for Fp {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Fp {
    type Output = Self;
    
    fn sub(self, rhs: Self) -> Self::Output {
        assert!(std::ptr::eq(self.params, rhs.params), "Field parameters must match");
        
        let result = if self.value >= rhs.value {
            self.value - rhs.value
        } else {
            &self.params.modulus + &self.value - rhs.value
        };
        
        Self {
            value: result,
            params: self.params,
        }
    }
}

impl SubAssign for Fp {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for Fp {
    type Output = Self;
    
    fn mul(self, rhs: Self) -> Self::Output {
        assert!(std::ptr::eq(self.params, rhs.params), "Field parameters must match");
        
        let result = self.params.montgomery_reduce(&self.value, &rhs.value);
        
        Self {
            value: result,
            params: self.params,
        }
    }
}

impl MulAssign for Fp {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Div for Fp {
    type Output = Self;
    
    fn div(self, rhs: Self) -> Self::Output {
        let inverse = rhs.inverse().expect("Division by zero");
        self * inverse
    }
}

impl DivAssign for Fp {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl Neg for Fp {
    type Output = Self;
    
    fn neg(self) -> Self::Output {
        if self.is_zero() {
            return self;
        }
        
        Self {
            value: &self.params.modulus - &self.value,
            params: self.params,
        }
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_value())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_field_arithmetic() {
        let a = Fp::bn254_scalar(5u64);
        let b = Fp::bn254_scalar(7u64);
        
        // Addition
        let c = a + b;
        assert_eq!(c.get_value(), BigUint::from(12u64));
        
        // Subtraction
        let d = b - a;
        assert_eq!(d.get_value(), BigUint::from(2u64));
        
        // Multiplication
        let e = a * b;
        assert_eq!(e.get_value(), BigUint::from(35u64));
        
        // Inversion and division
        let a_inv = a.inverse().unwrap();
        let f = b / a;
        let expected = Fp::bn254_scalar(7u64) * a_inv;
        assert_eq!(f, expected);
        
        // Squaring
        let g = a.square();
        assert_eq!(g.get_value(), BigUint::from(25u64));
        
        // Exponentiation
        let h = a.pow(3);
        assert_eq!(h.get_value(), BigUint::from(125u64));
    }
    
    #[test]
    fn test_field_properties() {
        let zero = Fp::zero();
        let one = Fp::one();
        
        // Zero properties
        assert!(zero.is_zero());
        assert!(!zero.is_one());
        
        // One properties
        assert!(!one.is_zero());
        assert!(one.is_one());
        
        // Additive identity
        let a = Fp::bn254_scalar(42u64);
        assert_eq!(a + zero, a);
        assert_eq!(zero + a, a);
        
        // Multiplicative identity
        assert_eq!(a * one, a);
        assert_eq!(one * a, a);
        
        // Additive inverse
        let neg_a = -a;
        assert_eq!(a + neg_a, zero);
        
        // Multiplicative inverse
        let a_inv = a.inverse().unwrap();
        assert_eq!(a * a_inv, one);
    }
    
    #[test]
    fn test_batch_inversion() {
        let a = Fp::bn254_scalar(5u64);
        let b = Fp::bn254_scalar(7u64);
        let c = Fp::bn254_scalar(11u64);
        
        let a_inv = a.inverse().unwrap();
        let b_inv = b.inverse().unwrap();
        let c_inv = c.inverse().unwrap();
        
        let mut elements = vec![a, b, c];
        Fp::batch_invert(&mut elements).unwrap();
        
        assert_eq!(elements[0], a_inv);
        assert_eq!(elements[1], b_inv);
        assert_eq!(elements[2], c_inv);
    }
}