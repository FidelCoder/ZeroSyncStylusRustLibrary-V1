use crate::arithmetic::field::Fp;
use num_bigint::BigUint;
use std::ops::{Add, Mul, Neg};
use std::str::FromStr;
use crate::arithmetic::traits::Field;
use num_traits::Zero;

/// BN254 elliptic curve implementation
#[derive(Debug, Clone)]
pub struct BN254 {
    /// The prime field modulus
    pub modulus: BigUint,
    /// Coefficient A in curve equation y² = x³ + Ax + B
    pub a: Fp,
    /// Coefficient B in curve equation y² = x³ + Ax + B
    pub b: Fp,
}

/// A point in G1 represented in affine coordinates
#[derive(Debug, Clone, PartialEq)]
pub struct G1Affine {
    pub x: Fp,
    pub y: Fp,
    pub infinity: bool,
}

/// Represents an element in the quadratic extension field Fp2
#[derive(Debug, Clone, PartialEq)]
pub struct Fp2 {
    /// Coefficient of 1 (real part)
    pub c0: Fp,
    /// Coefficient of u (imaginary part)
    pub c1: Fp,
}

/// A point in G2 represented in affine coordinates
#[derive(Debug, Clone, PartialEq)]
pub struct G2Affine {
    pub x: Fp2,
    pub y: Fp2,
    pub infinity: bool,
}

impl BN254 {
    /// Creates a new BN254 curve instance
    pub fn new() -> Self {
        // BN254 parameters
        let modulus = BigUint::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583"
        ).unwrap();
        let a = Fp::new(BigUint::from(0u32), modulus.clone());
        let b = Fp::new(BigUint::from(3u32), modulus.clone());
        
        Self { modulus, a, b }
    }
    
    /// Returns the generator point for G1
    pub fn g1_generator(&self) -> G1Affine {
        let x = Fp::new(
            BigUint::from_str("1").unwrap(),
            self.modulus.clone()
        );
        let y = Fp::new(
            BigUint::from_str("2").unwrap(),
            self.modulus.clone()
        );
        
        G1Affine { x, y, infinity: false }
    }
    
    /// Returns the generator point for G2
    pub fn g2_generator(&self) -> G2Affine {
        // BN254 G2 generator point (simplified values for example)
        let x0 = Fp::new(
            BigUint::from_str("10857046999023057135944570762232829481370756359578518086990519993285655852781").unwrap(),
            self.modulus.clone()
        );
        let x1 = Fp::new(
            BigUint::from_str("11559732032986387107991004021392285783925812861821192530917403151452391805634").unwrap(),
            self.modulus.clone()
        );
        
        let y0 = Fp::new(
            BigUint::from_str("8495653923123431417604973247489272438418190587263600148770280649306958101930").unwrap(),
            self.modulus.clone()
        );
        let y1 = Fp::new(
            BigUint::from_str("4082367875863433681332203403145435568316851327593401208105741076214120093531").unwrap(),
            self.modulus.clone()
        );
        
        G2Affine {
            x: Fp2 { c0: x0, c1: x1 },
            y: Fp2 { c0: y0, c1: y1 },
            infinity: false
        }
    }
    
    /// Checks if a point is on the curve
    pub fn is_on_curve(&self, point: &G1Affine) -> bool {
        if point.infinity {
            return true;
        }
        
        // y² = x³ + ax + b
        let x3 = point.x.clone() * point.x.clone() * point.x.clone();
        let ax = self.a.clone() * point.x.clone();
        let rhs = x3 + ax + self.b.clone();
        let lhs = point.y.clone() * point.y.clone();
        
        lhs == rhs
    }
    
    /// Checks if a G2 point is on the curve
    pub fn is_on_curve_g2(&self, point: &G2Affine) -> bool {
        if point.infinity {
            return true;
        }
        
        // y² = x³ + a*x + b, where a=0 and b=3/(u+9) in Fp2
        let x3 = point.x.mul(&point.x).mul(&point.x);
        
        // Calculate b in Fp2: b = 3/(u+9)
        let b_fp2 = self.twist_constant_b();
        
        // y² = x³ + b
        let rhs = x3.add(&b_fp2);
        let lhs = point.y.mul(&point.y);
        
        lhs == rhs
    }
    
    /// Returns the twisted curve constant b in Fp2
    fn twist_constant_b(&self) -> Fp2 {
        // For BN254, the twisted curve has b' = b / (u+9) = 3 / (u+9)
        // This is represented as an element in Fp2
        // Simplified implementation - actual values would be computed from the curve parameters
        let b_c0 = Fp::new(
            BigUint::from_str("19485874751759354771024239261021720505790618469301721065564631296452457478373").unwrap(),
            self.modulus.clone()
        );
        let b_c1 = Fp::new(
            BigUint::from_str("266929791119991161246907387137283842545076965332900288569378510910307636690").unwrap(),
            self.modulus.clone()
        );
        
        Fp2 { c0: b_c0, c1: b_c1 }
    }
}

impl G1Affine {
    /// Creates the identity point (point at infinity)
    pub fn identity(modulus: &BigUint) -> Self {
        let zero = Fp::new(BigUint::from(0u32), modulus.clone());
        Self {
            x: zero.clone(),
            y: zero,
            infinity: true,
        }
    }
    
    /// Point doubling with lazy reduction
    pub fn double(&self) -> Self {
        if self.infinity {
            return self.clone();
        }
        
        // Formula: λ = (3x²) / (2y)
        let x_squared = self.x.clone() * self.x.clone();
        let three = Fp::new(BigUint::from(3u32), self.modulus());
        let numerator = three * x_squared;
        
        let two = Fp::new(BigUint::from(2u32), self.modulus());
        let denominator = two.clone() * self.y.clone();
        
        // Check if denominator is zero to avoid division by zero
        if denominator.is_zero() {
            return Self::identity(&self.modulus());
        }
        
        let lambda = numerator * denominator.inverse().unwrap();
        
        // x' = λ² - 2x
        let lambda_squared = lambda.clone() * lambda.clone();
        let two_x = two * self.x.clone();
        let x3 = lambda_squared - two_x;
        
        // y' = λ(x - x') - y
        let x_diff = self.x.clone() - x3.clone();
        let lambda_x_diff = lambda * x_diff;
        let y3 = lambda_x_diff - self.y.clone();
        
        Self {
            x: x3,
            y: y3,
            infinity: false,
        }
    }
    
    /// Get the modulus of the field
    pub fn modulus(&self) -> BigUint {
        self.x.modulus()
    }
    
    /// Windowed scalar multiplication using a window size of 4 bits
    pub fn windowed_mul(&self, scalar: &BigUint) -> Self {
        const WINDOW_SIZE: usize = 4;
        const WINDOW_MASK: u64 = (1 << WINDOW_SIZE) - 1;
        
        // Precompute points for each window
        let mut window_points = vec![Self::identity(&self.modulus()); 1 << WINDOW_SIZE];
        window_points[1] = self.clone();
        for i in 2..(1 << WINDOW_SIZE) {
            window_points[i] = window_points[i-1].clone() + self.clone();
        }
        
        let mut result = Self::identity(&self.modulus());
        let mut temp = self.clone();
        
        // Process scalar in windows
        let mut remaining_scalar = scalar.clone();
        while !remaining_scalar.is_zero() {
            let window = (&remaining_scalar & BigUint::from(WINDOW_MASK)).to_u64_digits()[0];
            if window != 0 {
                result = result + window_points[window as usize].clone();
            }
            
            // Shift the window
            for _ in 0..WINDOW_SIZE {
                temp = temp.double();
            }
            remaining_scalar >>= WINDOW_SIZE;
        }
        
        result
    }
}

impl Add for G1Affine {
    type Output = Self;
    
    fn add(self, other: Self) -> Self {
        // Check for identity points
        if self.infinity {
            return other;
        }
        if other.infinity {
            return self;
        }
        
        // Check for point doubling
        if self.x == other.x && self.y == other.y {
            return self.double();
        }
        
        // Check for inverses (P + (-P) = O)
        let neg_y = self.y.clone().neg();
        if self.x == other.x && neg_y == other.y {
            return Self::identity(&self.modulus());
        }
        
        // Addition formula: λ = (y2 - y1) / (x2 - x1)
        let y_diff = other.y - self.y.clone();
        let x_diff = other.x.clone() - self.x.clone();
        
        // Check if x_diff is zero to avoid division by zero
        if x_diff.is_zero() {
            return Self::identity(&self.modulus());
        }
        
        let lambda = y_diff * x_diff.inverse().unwrap();
        
        // x3 = λ² - x1 - x2
        let lambda_squared = lambda.clone() * lambda.clone();
        let x3 = lambda_squared - self.x.clone() - other.x.clone();
        
        // y3 = λ(x1 - x3) - y1
        let x_diff = self.x - x3.clone();
        let y3 = lambda * x_diff - self.y;
        
        Self {
            x: x3,
            y: y3,
            infinity: false,
        }
    }
}

impl Neg for G1Affine {
    type Output = Self;
    
    fn neg(self) -> Self {
        if self.infinity {
            return self;
        }
        
        Self {
            x: self.x,
            y: self.y.neg(),
            infinity: false,
        }
    }
}

impl Mul<u64> for G1Affine {
    type Output = Self;
    
    fn mul(self, scalar: u64) -> Self {
        if scalar == 0 || self.infinity {
            return Self::identity(&self.modulus());
        }
        
        if scalar == 1 {
            return self;
        }
        
        // Double-and-add algorithm with lazy reduction
        let mut result = Self::identity(&self.modulus());
        let mut temp = self.clone();
        let mut n = scalar;
        
        while n > 0 {
            if n & 1 == 1 {
                result = result + temp.clone();
            }
            temp = temp.double();
            n >>= 1;
        }
        
        result
    }
}

impl G2Affine {
    /// Creates the identity point (point at infinity)
    pub fn identity(modulus: &BigUint) -> Self {
        let zero = Fp2::zero(modulus);
        Self {
            x: zero.clone(),
            y: zero,
            infinity: true,
        }
    }
    
    /// Point doubling with lazy reduction
    pub fn double(&self) -> Self {
        if self.infinity {
            return self.clone();
        }
        
        // Formula: λ = (3x²) / (2y)
        let x_squared = self.x.square();
        
        // Create field element for 3
        let three_c0 = Fp::new(BigUint::from(3u32), self.x.c0.modulus());
        let three_c1 = Fp::new(BigUint::from(0u32), self.x.c0.modulus()); 
        let three = Fp2::new(three_c0, three_c1);
        
        let numerator = three.mul(&x_squared);
        
        // Create field element for 2
        let two_c0 = Fp::new(BigUint::from(2u32), self.x.c0.modulus());
        let two_c1 = Fp::new(BigUint::from(0u32), self.x.c0.modulus());
        let two = Fp2::new(two_c0, two_c1);
        
        let denominator = two.mul(&self.y);
        
        // Check if denominator is zero to avoid division by zero
        if denominator.c0.is_zero() && denominator.c1.is_zero() {
            return Self::identity(&self.x.c0.modulus());
        }
        
        let lambda = match denominator.inverse() {
            Some(inv) => numerator.mul(&inv),
            None => return Self::identity(&self.x.c0.modulus())
        };
        
        // x' = λ² - 2x
        let lambda_squared = lambda.square();
        let two_x = two.mul(&self.x);
        let x3 = lambda_squared.sub(&two_x);
        
        // y' = λ(x - x') - y
        let x_diff = self.x.sub(&x3);
        let lambda_x_diff = lambda.mul(&x_diff);
        let y3 = lambda_x_diff.sub(&self.y);
        
        Self {
            x: x3,
            y: y3,
            infinity: false,
        }
    }
    
    /// Windowed scalar multiplication using a window size of 4 bits
    pub fn windowed_mul(&self, scalar: &BigUint) -> Self {
        const WINDOW_SIZE: usize = 4;
        const WINDOW_MASK: u64 = (1 << WINDOW_SIZE) - 1;
        
        // Precompute points for each window
        let mut window_points = vec![Self::identity(&self.x.c0.modulus()); 1 << WINDOW_SIZE];
        window_points[1] = self.clone();
        for i in 2..(1 << WINDOW_SIZE) {
            window_points[i] = window_points[i-1].clone() + self.clone();
        }
        
        let mut result = Self::identity(&self.x.c0.modulus());
        let mut temp = self.clone();
        
        // Process scalar in windows
        let mut remaining_scalar = scalar.clone();
        while !remaining_scalar.is_zero() {
            let window = (&remaining_scalar & BigUint::from(WINDOW_MASK)).to_u64_digits()[0];
            if window != 0 {
                result = result + window_points[window as usize].clone();
            }
            
            // Shift the window
            for _ in 0..WINDOW_SIZE {
                temp = temp.double();
            }
            remaining_scalar >>= WINDOW_SIZE;
        }
        
        result
    }
}

impl Add for G2Affine {
    type Output = Self;
    
    fn add(self, other: Self) -> Self {
        // Check for identity points
        if self.infinity {
            return other;
        }
        if other.infinity {
            return self;
        }
        
        // Check for point doubling
        if self.x == other.x && self.y == other.y {
            return self.double();
        }
        
        // Check for inverses (P + (-P) = O)
        let neg_y = self.y.neg();
        if self.x == other.x && neg_y == other.y {
            return Self::identity(&self.x.c0.modulus());
        }
        
        // Addition formula: λ = (y2 - y1) / (x2 - x1)
        let y_diff = other.y.sub(&self.y);
        let x_diff = other.x.sub(&self.x);
        
        // Check if x_diff is zero to avoid division by zero
        if x_diff.c0.is_zero() && x_diff.c1.is_zero() {
            return Self::identity(&self.x.c0.modulus());
        }
        
        let lambda = match x_diff.inverse() {
            Some(inv) => y_diff.mul(&inv),
            None => return Self::identity(&self.x.c0.modulus())
        };
        
        // x3 = λ² - x1 - x2
        let lambda_squared = lambda.square();
        let x3 = lambda_squared.sub(&self.x).sub(&other.x);
        
        // y3 = λ(x1 - x3) - y1
        let x_diff = self.x.sub(&x3);
        let y3 = lambda.mul(&x_diff).sub(&self.y);
        
        Self {
            x: x3,
            y: y3,
            infinity: false,
        }
    }
}

impl Neg for G2Affine {
    type Output = Self;
    
    fn neg(self) -> Self {
        if self.infinity {
            return self;
        }
        
        Self {
            x: self.x,
            y: self.y.neg(),
            infinity: self.infinity,
        }
    }
}

impl Mul<u64> for G2Affine {
    type Output = Self;
    
    fn mul(self, scalar: u64) -> Self {
        let mut result = Self::identity(&self.x.c0.modulus());
        let mut temp = self;
        let mut s = scalar;
        
        while s > 0 {
            if s & 1 == 1 {
                result = result + temp.clone();
            }
            temp = temp.double();
            s >>= 1;
        }
        
        result
    }
}

impl Fp2 {
    /// Create a new Fp2 element
    pub fn new(c0: Fp, c1: Fp) -> Self {
        // Ideally, we'd validate that both elements use the same modulus
        // But since we can't access the private field directly, we're omitting this check
        // in a production environment, we would need to design a better solution
        Self { c0, c1 }
    }
    
    /// Create zero in Fp2
    pub fn zero(modulus: &BigUint) -> Self {
        let zero = Fp::new(BigUint::from(0u32), modulus.clone());
        Self { c0: zero.clone(), c1: zero }
    }
    
    /// Create one in Fp2
    pub fn one(modulus: &BigUint) -> Self {
        let zero = Fp::new(BigUint::from(0u32), modulus.clone());
        let one = Fp::new(BigUint::from(1u32), modulus.clone());
        Self { c0: one, c1: zero }
    }
    
    /// Multiply two Fp2 elements
    pub fn mul(&self, other: &Self) -> Self {
        // (a + bu) * (c + du) = (ac - bd) + (ad + bc)u
        let a = &self.c0;
        let b = &self.c1;
        let c = &other.c0;
        let d = &other.c1;
        
        let ac = a.clone() * c.clone();
        let bd = b.clone() * d.clone();
        
        // Calculate (a+b)(c+d)
        let a_plus_b = a.clone() + b.clone();
        let c_plus_d = c.clone() + d.clone();
        let abcd = a_plus_b * c_plus_d;
        
        // For c0 = ac - bd, we need to handle modular subtraction
        // We'll use the modular subtraction built into the field implementation
        let c0 = ac.clone() - bd.clone();
        
        // For c1 = abcd - ac - bd, we'll use field operations
        let c1 = abcd - ac - bd;
        
        Self { c0, c1 }
    }
    
    /// Add two Fp2 elements
    pub fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.clone() + other.c0.clone(),
            c1: self.c1.clone() + other.c1.clone(),
        }
    }
    
    /// Subtract two Fp2 elements
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.clone() - other.c0.clone(),
            c1: self.c1.clone() - other.c1.clone(),
        }
    }
    
    /// Negate an Fp2 element
    pub fn neg(&self) -> Self {
        Self {
            c0: self.c0.clone().neg(),
            c1: self.c1.clone().neg(),
        }
    }
    
    /// Square an Fp2 element (optimized multiplication by self)
    pub fn square(&self) -> Self {
        // (a + bu)² = (a² - b²) + 2abu
        let a = &self.c0;
        let b = &self.c1;
        
        let a2 = a.clone() * a.clone();
        let b2 = b.clone() * b.clone();
        let ab = a.clone() * b.clone();
        
        // For c0 = a² - b², we'll use field operations
        let c0 = a2 - b2;
        
        // For c1 = 2ab, we'll use field operations
        let two = Fp::new(BigUint::from(2u32), a.modulus());
        let c1 = two * ab;
        
        Self { c0, c1 }
    }
    
    /// Compute the inverse of an Fp2 element
    pub fn inverse(&self) -> Option<Self> {
        // (a + bu)⁻¹ = (a - bu) / (a² + b²)
        let a = &self.c0;
        let b = &self.c1;
        
        let a2 = a.clone() * a.clone();
        let b2 = b.clone() * b.clone();
        let denom = a2 + b2;
        
        match denom.inverse() {
            Some(inv) => Some(Self {
                c0: a.clone() * inv.clone(),
                c1: (b.clone() * inv).neg(),
            }),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;
    
    #[test]
    fn test_point_on_curve() {
        let curve = BN254::new();
        let g1 = curve.g1_generator();
        
        assert!(curve.is_on_curve(&g1));
    }
    
    #[test]
    fn test_point_addition() {
        let curve = BN254::new();
        let g1 = curve.g1_generator();
        
        let p2 = g1.clone() + g1.clone();
        assert!(curve.is_on_curve(&p2));
        
        let p3 = p2 + g1.clone();
        assert!(curve.is_on_curve(&p3));
    }
    
    #[test]
    fn test_scalar_multiplication() {
        let curve = BN254::new();
        let g1 = curve.g1_generator();
        
        let p2 = g1.clone() * 2;
        let p2_add = g1.clone() + g1.clone();
        assert_eq!(p2, p2_add);
        
        let p3 = g1.clone() * 3;
        let p3_add = g1.clone() + g1.clone() + g1.clone();
        assert_eq!(p3, p3_add);
    }
    
    #[test]
    fn test_fp2_arithmetic() {
        let curve = BN254::new();
        let modulus = curve.modulus.clone();
        
        // Create some Fp2 elements
        let a = Fp::new(BigUint::from(2u32), modulus.clone());
        let b = Fp::new(BigUint::from(3u32), modulus.clone());
        let c = Fp::new(BigUint::from(5u32), modulus.clone());
        let d = Fp::new(BigUint::from(7u32), modulus.clone());
        
        let fp2_a = Fp2::new(a, b);
        let fp2_b = Fp2::new(c, d);
        
        // Test multiplication
        let prod = fp2_a.mul(&fp2_b);
        
        // (2 + 3u) * (5 + 7u) = (2*5 - 3*7) + (2*7 + 3*5)u = (10 - 21) + (14 + 15)u = -11 + 29u
        let expected_c0 = Fp::new(
            (BigUint::from(10u32) + &modulus - BigUint::from(21u32)) % &modulus,
            modulus.clone()
        );
        let expected_c1 = Fp::new(BigUint::from(29u32), modulus.clone());
        let expected = Fp2::new(expected_c0, expected_c1);
        
        assert_eq!(prod, expected);
        
        // Test squaring
        let sq = fp2_a.square();
        
        // (2 + 3u)² = (2² - 3²) + 2*2*3*u = (4 - 9) + 12u = -5 + 12u
        let expected_c0 = Fp::new(
            (BigUint::from(4u32) + &modulus - BigUint::from(9u32)) % &modulus,
            modulus.clone()
        );
        let expected_c1 = Fp::new(BigUint::from(12u32), modulus);
        let expected = Fp2::new(expected_c0, expected_c1);
        
        assert_eq!(sq, expected);
    }
    
    #[test]
    fn test_g2_generator() {
        let curve = BN254::new();
        let g2 = curve.g2_generator();
        
        // Due to the simplified implementation, we're just validating that this runs
        // A real implementation would verify that g2 is on the twisted curve
        assert!(!g2.infinity);
    }
    
    #[test]
    fn test_g2_point_addition() {
        let curve = BN254::new();
        let g = curve.g2_generator();
        
        // Test point addition
        let p2 = g.clone() + g.clone();
        assert!(curve.is_on_curve_g2(&p2));
        
        // Test point doubling
        let p2_double = g.clone().double();
        assert!(curve.is_on_curve_g2(&p2_double));
        assert_eq!(p2, p2_double);
        
        // Test identity point
        let identity = G2Affine::identity(&curve.modulus);
        let p_plus_identity = g.clone() + identity.clone();
        assert_eq!(p_plus_identity, g);
        
        // Test inverse points
        let neg_g = g.clone().neg();
        let p_plus_neg = g.clone() + neg_g;
        assert!(p_plus_neg.infinity);
    }
    
    #[test]
    fn test_g2_scalar_multiplication() {
        let curve = BN254::new();
        let g = curve.g2_generator();
        
        // Test scalar multiplication
        let p2 = g.clone() * 2;
        let p3 = g.clone() * 3;
        let p4 = g.clone() * 4;
        
        assert!(curve.is_on_curve_g2(&p2));
        assert!(curve.is_on_curve_g2(&p3));
        assert!(curve.is_on_curve_g2(&p4));
        
        // Verify scalar multiplication properties
        let p2_plus_p2 = p2.clone() + p2.clone();
        assert_eq!(p2_plus_p2, p4);
        
        let p2_plus_g = p2 + g;
        assert_eq!(p2_plus_g, p3);
    }
    
    #[test]
    fn test_g2_generator_correctness() {
        let curve = BN254::new();
        let g = curve.g2_generator();
        
        // Verify generator point is on the curve
        assert!(curve.is_on_curve_g2(&g));
        
        // Verify generator point is not the identity
        assert!(!g.infinity);
        
        // Verify generator point coordinates are non-zero
        assert!(!g.x.c0.is_zero());
        assert!(!g.x.c1.is_zero());
        assert!(!g.y.c0.is_zero());
        assert!(!g.y.c1.is_zero());
    }
    
    #[test]
    fn test_windowed_scalar_multiplication_g1() {
        let curve = BN254::new();
        let g = curve.g1_generator();
        
        // Test with various scalar values
        let scalars = vec![
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from(10u32),
            BigUint::from(100u32),
            BigUint::from_str("123456789").unwrap(),
        ];
        
        for scalar in scalars {
            let p1 = g.windowed_mul(&scalar);
            let p2 = g.clone() * scalar.to_u64_digits()[0];
            assert_eq!(p1, p2);
            assert!(curve.is_on_curve(&p1));
        }
    }
    
    #[test]
    fn test_windowed_scalar_multiplication_g2() {
        let curve = BN254::new();
        let g = curve.g2_generator();
        
        // Test with various scalar values
        let scalars = vec![
            BigUint::from(1u32),
            BigUint::from(2u32),
            BigUint::from(10u32),
            BigUint::from(100u32),
            BigUint::from_str("123456789").unwrap(),
        ];
        
        for scalar in scalars {
            let p1 = g.windowed_mul(&scalar);
            let p2 = g.clone() * scalar.to_u64_digits()[0];
            assert_eq!(p1, p2);
            assert!(curve.is_on_curve_g2(&p1));
        }
    }
} 