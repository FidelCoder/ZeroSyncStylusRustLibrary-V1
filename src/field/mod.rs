//! Finite field implementations optimized for Stylus
//! 
//! This module provides field arithmetic operations with a focus on 
//! gas efficiency within the Stylus environment.

use std::fmt::{Debug, Display};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub mod fp;
pub mod extension;

/// Error types for field operations
#[derive(Debug, thiserror::Error)]
pub enum FieldError {
    #[error("Division by zero")]
    DivisionByZero,
    
    #[error("Invalid field element")]
    InvalidElement,
    
    #[error("Modulus must be prime")]
    NonPrimeModulus,
    
    #[error("Field operation error: {0}")]
    OperationError(String),
}

/// Result type for field operations
pub type FieldResult<T> = Result<T, FieldError>;

/// Field element trait providing the interface for finite field arithmetic
pub trait Field: 
    Sized 
    + Clone 
    + Copy 
    + Debug 
    + Display
    + PartialEq 
    + Eq
    + Add<Output = Self> 
    + AddAssign
    + Sub<Output = Self> 
    + SubAssign
    + Mul<Output = Self> 
    + MulAssign
    + Div<Output = Self> 
    + DivAssign
    + Neg<Output = Self>
{
    /// Returns the zero element of the field
    fn zero() -> Self;
    
    /// Returns the multiplicative identity (one) of the field
    fn one() -> Self;
    
    /// Generates a random element of the field
    fn random() -> Self;
    
    /// Checks if this element is zero
    fn is_zero(&self) -> bool;
    
    /// Checks if this element is one
    fn is_one(&self) -> bool;
    
    /// Computes the multiplicative inverse of this element if it exists
    fn inverse(&self) -> FieldResult<Self>;
    
    /// Computes the square of this element
    fn square(&self) -> Self {
        *self * *self
    }
    
    /// Computes this element raised to the given power
    fn pow(&self, exp: u64) -> Self {
        if exp == 0 {
            return Self::one();
        }
        
        let mut base = *self;
        let mut result = Self::one();
        let mut exp = exp;
        
        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base.square();
            exp >>= 1;
        }
        
        result
    }
    
    /// Computes the square root of this element if it exists
    fn sqrt(&self) -> FieldResult<Self>;
    
    /// Computes a double of this element
    fn double(&self) -> Self {
        *self + *self
    }
    
    /// Batch inversion of multiple field elements
    /// More efficient than individual inversions
    fn batch_invert(elements: &mut [Self]) -> FieldResult<()> {
        if elements.is_empty() {
            return Ok(());
        }
        
        // Use Montgomery's trick for batch inversion
        let n = elements.len();
        let mut products = Vec::with_capacity(n);
        let mut acc = Self::one();
        
        // Compute products
        for element in elements.iter() {
            if element.is_zero() {
                return Err(FieldError::DivisionByZero);
            }
            products.push(acc);
            acc = acc * *element;
        }
        
        // Compute inverse of the accumulated product
        let mut inv = acc.inverse()?;
        
        // Compute individual inverses
        for i in (0..n).rev() {
            let tmp = elements[i] * inv;
            elements[i] = products[i] * inv;
            inv = tmp;
        }
        
        Ok(())
    }
}