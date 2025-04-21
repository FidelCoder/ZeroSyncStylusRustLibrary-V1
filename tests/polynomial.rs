use zerosync::polynomial::{Polynomial, evaluate_polynomial};
use num_bigint::BigUint;

#[test]
fn test_polynomial_creation() {
    let coefficients = vec![1u64, 2, 3];  // x^2 + 2x + 3
    let poly = Polynomial::new(coefficients.clone());
    assert_eq!(poly.degree(), 2);
    assert_eq!(poly.coefficients(), &coefficients);
}

#[test]
fn test_polynomial_addition() {
    let poly1 = Polynomial::new(vec![1u64, 2, 3]);  // x^2 + 2x + 3
    let poly2 = Polynomial::new(vec![4u64, 5, 6]);  // x^2 + 5x + 6
    
    let result = &poly1 + &poly2;
    assert_eq!(result.coefficients(), &[5u64, 7, 9]);  // (x^2 + 2x + 3) + (x^2 + 5x + 6) = 2x^2 + 7x + 9
}

#[test]
fn test_polynomial_multiplication() {
    let poly1 = Polynomial::new(vec![1u64, 2]);     // 2x + 1
    let poly2 = Polynomial::new(vec![3u64, 4]);     // 4x + 3
    
    let result = &poly1 * &poly2;
    // (2x + 1)(4x + 3) = 8x^2 + 10x + 3
    assert_eq!(result.coefficients(), &[3u64, 10, 8]);
}

#[test]
fn test_polynomial_evaluation() {
    let poly = Polynomial::new(vec![1u64, 2, 3]);  // 3x^2 + 2x + 1
    let x = BigUint::from(2u32);
    
    let result = evaluate_polynomial(&poly, &x);
    // At x = 2: 3(2^2) + 2(2) + 1 = 12 + 4 + 1 = 17
    assert_eq!(result, BigUint::from(17u32));
}

#[test]
fn test_polynomial_degree() {
    let poly = Polynomial::new(vec![1u64, 0, 0, 4]);  // 4x^3 + 1
    assert_eq!(poly.degree(), 3);
    
    let zero_poly = Polynomial::new(vec![0u64]);
    assert_eq!(zero_poly.degree(), 0);
}

#[test]
fn test_polynomial_zero() {
    let zero = Polynomial::zero();
    assert_eq!(zero.coefficients(), &[0u64]);
    assert_eq!(zero.degree(), 0);
    
    let poly = Polynomial::new(vec![1u64, 2, 3]);
    assert_eq!(&poly + &zero, poly);
}

#[test]
fn test_polynomial_scalar_multiplication() {
    let poly = Polynomial::new(vec![1u64, 2, 3]);  // 3x^2 + 2x + 1
    let scalar = 2u64;
    
    let result = poly.scalar_mul(scalar);
    assert_eq!(result.coefficients(), &[2u64, 4, 6]);  // 6x^2 + 4x + 2
}

#[test]
fn test_polynomial_composition() {
    let outer = Polynomial::new(vec![1u64, 1]);     // x + 1
    let inner = Polynomial::new(vec![2u64, 1]);     // x + 2
    
    let result = outer.compose(&inner);
    // (x + 2) + 1 = x + 3
    assert_eq!(result.coefficients(), &[3u64, 1]);
}

#[test]
fn test_polynomial_derivative() {
    let poly = Polynomial::new(vec![1u64, 2, 3]);  // 3x^2 + 2x + 1
    let derivative = poly.derivative();
    // Derivative: 6x + 2
    assert_eq!(derivative.coefficients(), &[2u64, 6]);
}

#[test]
fn test_polynomial_interpolation() {
    let points = vec![(0u64, 1u64), (1u64, 2u64), (2u64, 4u64)];
    let poly = Polynomial::interpolate(&points);
    
    // Verify the polynomial passes through all points
    for (x, y) in points {
        let result = evaluate_polynomial(&poly, &BigUint::from(x));
        assert_eq!(result, BigUint::from(y));
    }
}

#[test]
fn test_polynomial_properties() {
    let p1 = Polynomial::new(vec![1u64, 2]);  // 2x + 1
    let p2 = Polynomial::new(vec![3u64, 4]);  // 4x + 3
    let p3 = Polynomial::new(vec![5u64, 6]);  // 6x + 5
    
    // Test associativity: (p1 + p2) + p3 = p1 + (p2 + p3)
    let left = (&(&p1 + &p2) + &p3).coefficients().to_vec();
    let right = (&p1 + &(&p2 + &p3)).coefficients().to_vec();
    assert_eq!(left, right);
    
    // Test distributivity: p1 * (p2 + p3) = (p1 * p2) + (p1 * p3)
    let left = (&p1 * &(&p2 + &p3)).coefficients().to_vec();
    let right = (&(&p1 * &p2) + &(&p1 * &p3)).coefficients().to_vec();
    assert_eq!(left, right);
} 