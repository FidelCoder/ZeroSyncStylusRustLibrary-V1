use std::ops::{Add, Mul};
use crate::arithmetic::traits::Field;

/// Represents a univariate polynomial over a field
#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial<F: Field> {
    /// Coefficients of the polynomial in ascending order of degree
    pub coefficients: Vec<F>,
}

impl<F: Field> Polynomial<F> {
    /// Creates a new polynomial from its coefficients
    pub fn new(mut coefficients: Vec<F>) -> Self {
        // Remove trailing zeros
        while coefficients.len() > 1 && coefficients.last().unwrap().is_zero() {
            coefficients.pop();
        }
        // Ensure at least one coefficient (zero polynomial has one zero coefficient)
        if coefficients.is_empty() {
            coefficients.push(F::zero());
        }
        Self { coefficients }
    }

    /// Returns the degree of the polynomial
    pub fn degree(&self) -> usize {
        let mut deg = self.coefficients.len() - 1;
        // Find the highest non-zero coefficient
        while deg > 0 && self.coefficients[deg].is_zero() {
            deg -= 1;
        }
        deg
    }

    /// Returns a reference to the polynomial's coefficients
    pub fn coefficients(&self) -> &[F] {
        &self.coefficients
    }

    /// Creates a zero polynomial
    pub fn zero() -> Self {
        Self::new(vec![F::zero()])
    }

    /// Returns true if this is the zero polynomial
    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.coefficients[0].is_zero()
    }

    /// Computes the formal derivative of the polynomial
    pub fn derivative(&self) -> Self {
        if self.degree() == 0 {
            return Self::zero();
        }

        let mut result = Vec::with_capacity(self.degree());
        for (i, coeff) in self.coefficients.iter().skip(1).enumerate() {
            let mut term = coeff.clone();
            // Multiply by the power (i + 1)
            let mut power = F::one();
            for _ in 0..i+1 {
                power = power + F::one();
            }
            term = term * power;
            result.push(term);
        }
        Self::new(result)
    }

    /// Interpolates a polynomial from a set of points
    pub fn interpolate(points: &[(F, F)]) -> Self {
        let n = points.len();
        if n == 0 {
            return Self::zero();
        }

        // Check for distinct x-coordinates
        for i in 0..n {
            for j in i+1..n {
                if points[i].0 == points[j].0 {
                    panic!("points must have distinct x-coordinates");
                }
            }
        }

        let mut result = Self::zero();
        for i in 0..n {
            let mut term = Self::new(vec![points[i].1.clone()]);
            
            for j in 0..n {
                if i != j {
                    let denom = points[i].0.clone() - points[j].0.clone();
                    let denom_inv = denom.inverse().expect("points must be distinct");
                    
                    let factor = Self::new(vec![
                        -points[j].0.clone() * denom_inv.clone(),
                        denom_inv,
                    ]);
                    term = &term * &factor;
                }
            }
            
            result = &result + &term;
        }
        
        result
    }
}

impl<'a, F: Field> Add for &'a Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, other: Self) -> Self::Output {
        let max_len = self.coefficients.len().max(other.coefficients.len());
        let mut result = vec![F::zero(); max_len];

        for (i, coeff) in self.coefficients.iter().enumerate() {
            result[i] = coeff.clone();
        }

        for (i, coeff) in other.coefficients.iter().enumerate() {
            result[i] = result[i].clone() + coeff.clone();
        }

        Polynomial::new(result)
    }
}

impl<'a, F: Field> Mul for &'a Polynomial<F> {
    type Output = Polynomial<F>;

    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return Polynomial::zero();
        }

        let n = self.coefficients.len();
        let m = other.coefficients.len();
        let mut result = vec![F::zero(); n + m - 1];

        for i in 0..n {
            for j in 0..m {
                let prod = self.coefficients[i].clone() * other.coefficients[j].clone();
                result[i + j] = result[i + j].clone() + prod;
            }
        }

        Polynomial::new(result)
    }
}

/// Evaluates a polynomial at a given point using Horner's method
pub fn evaluate_polynomial<F: Field>(poly: &Polynomial<F>, x: &F) -> F {
    let mut result = F::zero();
    for coeff in poly.coefficients.iter().rev() {
        result = result * x.clone() + coeff.clone();
    }
    result
} 