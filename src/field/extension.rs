//! Extension field implementations for quadratic and cubic extensions
//! Optimized for gas efficiency in Stylus environment

use crate::field::{Field, FieldError, FieldResult};
use std::fmt;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Quadratic extension field element F_p^2 = F_p[x]/(x^2 + 1)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp2<F: Field> {
    /// c0 + c1 * x representation, where x^2 = -1 (or other quadratic non-residue)
    pub c0: F,
    pub c1: F,
}

impl<F: Field> Fp2<F> {
    /// Create a new Fp2 element from coefficients
    pub fn new(c0: F, c1: F) -> Self {
        Self { c0, c1 }
    }
    
    /// Frobenius map (p-power raising)
    /// For p ≡ 3 (mod 4), Frobenius(a + bi) = a - bi
    pub fn frobenius(&self) -> Self {
        Self::new(self.c0, -self.c1)
    }
    
    /// Multiplication by quadratic non-residue
    /// Assuming x^2 = -1, mul_by_nonresidue(a + bx) = -b + ax
    pub fn mul_by_nonresidue(&self) -> Self {
        Self::new(-self.c1, self.c0)
    }
}

impl<F: Field> Field for Fp2<F> {
    fn zero() -> Self {
        Self::new(F::zero(), F::zero())
    }
    
    fn one() -> Self {
        Self::new(F::one(), F::zero())
    }
    
    fn random() -> Self {
        Self::new(F::random(), F::random())
    }
    
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }
    
    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }
    
    fn inverse(&self) -> FieldResult<Self> {
        if self.is_zero() {
            return Err(FieldError::DivisionByZero);
        }
        
        // (a + bi)^(-1) = (a - bi) / (a^2 + b^2)
        let norm = (self.c0 * self.c0) + (self.c1 * self.c1);
        let inv_norm = norm.inverse()?;
        
        Ok(Self::new(
            self.c0 * inv_norm,
            -self.c1 * inv_norm
        ))
    }
    
    fn sqrt(&self) -> FieldResult<Self> {
        if self.is_zero() {
            return Ok(*self);
        }
        
        // Algorithm for computing sqrt in Fp2
        // For c0 + c1*i, we compute (a + bi) such that (a + bi)^2 = c0 + c1*i
        
        // Check if self is a square
        // For Fp2, every element is a square, but implementation is more complex
        
        // Simplified algorithm for the case where p ≡ 3 (mod 4)
        let norm = (self.c0 * self.c0) + (self.c1 * self.c1);
        let norm_sqrt = norm.sqrt()?;
        
        // Compute a = sqrt((norm + c0) / 2)
        let a_squared = (norm + self.c0) / (F::one() + F::one());
        let a = a_squared.sqrt()?;
        
        // Compute b = c1 / (2a)
        let two_a = a + a;
        if two_a.is_zero() {
            return Err(FieldError::OperationError("Cannot compute square root".into()));
        }
        
        let b = self.c1 / two_a;
        
        // Verify the result
        let result = Self::new(a, b);
        if result.square() != *self && (-result).square() != *self {
            return Err(FieldError::OperationError("Square root computation failed".into()));
        }
        
        Ok(result)
    }
}

// Implement arithmetic operations for Fp2
impl<F: Field> Add for Fp2<F> {
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(
            self.c0 + rhs.c0,
            self.c1 + rhs.c1,
        )
    }
}

impl<F: Field> AddAssign for Fp2<F> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<F: Field> Sub for Fp2<F> {
    type Output = Self;
    
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(
            self.c0 - rhs.c0,
            self.c1 - rhs.c1,
        )
    }
}

impl<F: Field> SubAssign for Fp2<F> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<F: Field> Mul for Fp2<F> {
    type Output = Self;
    
    fn mul(self, rhs: Self) -> Self::Output {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        // Using Karatsuba multiplication to reduce multiplications
        let a = self.c0;
        let b = self.c1;
        let c = rhs.c0;
        let d = rhs.c1;
        
        let ac = a * c;
        let bd = b * d;
        let abcd = (a + b) * (c + d);
        
        Self::new(
            ac - bd,
            abcd - ac - bd,
        )
    }
}

impl<F: Field> MulAssign for Fp2<F> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<F: Field> Div for Fp2<F> {
    type Output = Self;
    
    fn div(self, rhs: Self) -> Self::Output {
        let inverse = rhs.inverse().expect("Division by zero");
        self * inverse
    }
}

impl<F: Field> DivAssign for Fp2<F> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<F: Field> Neg for Fp2<F> {
    type Output = Self;
    
    fn neg(self) -> Self::Output {
        Self::new(-self.c0, -self.c1)
    }
}

impl<F: Field> fmt::Display for Fp2<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({} + {}i)", self.c0, self.c1)
    }
}

/// Fp6 field element as a cubic extension of Fp2
/// F_p^6 = F_p^2[y]/(y^3 - ξ) where ξ is a non-residue in F_p^2
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp6<F: Field> {
    /// c0 + c1 * y + c2 * y^2 representation
    pub c0: Fp2<F>,
    pub c1: Fp2<F>,
    pub c2: Fp2<F>,
}

impl<F: Field> Fp6<F> {
    /// Create a new Fp6 element from coefficients
    pub fn new(c0: Fp2<F>, c1: Fp2<F>, c2: Fp2<F>) -> Self {
        Self { c0, c1, c2 }
    }
    
    /// Frobenius map (p-power raising)
    pub fn frobenius(&self) -> Self {
        // Implementation depends on specific field parameters
        // Simplified version here
        Self::new(
            self.c0.frobenius(),
            self.c1.frobenius().mul_by_nonresidue(),
            self.c2.frobenius().mul_by_nonresidue().mul_by_nonresidue(),
        )
    }
}

impl<F: Field> Field for Fp6<F> {
    fn zero() -> Self {
        Self::new(Fp2::zero(), Fp2::zero(), Fp2::zero())
    }
    
    fn one() -> Self {
        Self::new(Fp2::one(), Fp2::zero(), Fp2::zero())
    }
    
    fn random() -> Self {
        Self::new(Fp2::random(), Fp2::random(), Fp2::random())
    }
    
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }
    
    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero() && self.c2.is_zero()
    }
    
    fn inverse(&self) -> FieldResult<Self> {
        if self.is_zero() {
            return Err(FieldError::DivisionByZero);
        }
        
        // Using the formula for the inverse of a cubic extension
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;
        
        // Precompute some values
        let aa = a * a;
        let bb = b * b;
        let cc = c * c;
        let ac = a * c;
        let ab = a * b