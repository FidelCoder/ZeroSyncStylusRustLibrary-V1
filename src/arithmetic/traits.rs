use std::ops::{Add, Sub, Mul, Div, Neg};
use num_traits::{Zero, One};

/// Trait for field elements with basic arithmetic operations
pub trait Field:
    Sized
    + Clone
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + Zero
    + One
{
    /// Returns the characteristic of the field
    fn characteristic() -> Vec<u64>;

    /// Returns the multiplicative inverse of this element
    fn inverse(&self) -> Option<Self>;

    /// Squares this element
    fn square(&self) -> Self {
        self.clone() * self.clone()
    }

    /// Raises this element to a power
    fn pow(&self, exp: u64) -> Self;
}

/// Trait for prime fields with modular arithmetic
pub trait PrimeField: Field {
    /// The modulus of the field
    fn modulus() -> Vec<u64>;
    
    /// Returns the field element in Montgomery form
    fn to_montgomery(&self) -> Self;
    
    /// Returns the field element from Montgomery form
    fn from_montgomery(&self) -> Self;
} 