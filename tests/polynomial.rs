use num_bigint::BigUint;
use num_traits::Num;
use zerosync::arithmetic::field::Fp;
use zerosync::polynomial::{Polynomial, evaluate_polynomial};

const TEST_MODULUS: &str = "17";  // Small prime for testing

#[test]
fn test_polynomial_creation() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let coefficients = vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(0u64), modulus.clone()),
        Fp::new(BigUint::from(0u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ];  // 4x^3 + 1
    
    let poly = Polynomial::new(coefficients);
    assert_eq!(poly.degree(), 3);
}

#[test]
fn test_polynomial_zero() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let zero_poly = Polynomial::new(vec![
        Fp::new(BigUint::from(0u64), modulus.clone())
    ]);
    assert_eq!(zero_poly.degree(), 0);
}

#[test]
fn test_polynomial_evaluation() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let poly = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(3u64), modulus.clone()),
    ]);  // 3x^2 + 2x + 1
    
    let x = Fp::new(BigUint::from(2u64), modulus.clone());
    let result = evaluate_polynomial(&poly, &x);
    
    // At x = 2: 3(2^2) + 2(2) + 1 = 12 + 4 + 1 = 17 ≡ 0 (mod 17)
    let expected = Fp::new(BigUint::from(0u64), modulus);
    assert_eq!(result, expected);
}

#[test]
fn test_polynomial_addition() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let p1 = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
    ]);  // 2x + 1
    
    let p2 = Polynomial::new(vec![
        Fp::new(BigUint::from(3u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ]);  // 4x + 3
    
    let result = &p1 + &p2;
    let expected = Polynomial::new(vec![
        Fp::new(BigUint::from(4u64), modulus.clone()),  // (1 + 3) mod 17
        Fp::new(BigUint::from(6u64), modulus.clone()),  // (2 + 4) mod 17
    ]);
    assert_eq!(result, expected);
}

#[test]
fn test_polynomial_multiplication() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let p1 = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
    ]);  // 2x + 1
    
    let p2 = Polynomial::new(vec![
        Fp::new(BigUint::from(3u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ]);  // 4x + 3
    
    let result = &p1 * &p2;
    // (2x + 1)(4x + 3) = 8x^2 + 10x + 3
    let expected = Polynomial::new(vec![
        Fp::new(BigUint::from(3u64), modulus.clone()),
        Fp::new(BigUint::from(10u64), modulus.clone()),
        Fp::new(BigUint::from(8u64), modulus.clone()),
    ]);
    assert_eq!(result, expected);
}

#[test]
fn test_polynomial_associativity() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let p1 = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
    ]);  // 2x + 1
    
    let p2 = Polynomial::new(vec![
        Fp::new(BigUint::from(3u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ]);  // 4x + 3
    
    let p3 = Polynomial::new(vec![
        Fp::new(BigUint::from(5u64), modulus.clone()),
        Fp::new(BigUint::from(6u64), modulus.clone()),
    ]);  // 6x + 5
    
    // Test associativity of addition: (p1 + p2) + p3 = p1 + (p2 + p3)
    let left = &(&p1 + &p2) + &p3;
    let right = &p1 + &(&p2 + &p3);
    assert_eq!(left, right);
    
    // Test distributivity: p1 * (p2 + p3) = (p1 * p2) + (p1 * p3)
    let left = &p1 * &(&p2 + &p3);
    let right = &(&p1 * &p2) + &(&p1 * &p3);
    assert_eq!(left, right);
}

#[test]
fn test_polynomial_degree() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let coefficients = vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(0u64), modulus.clone()),
        Fp::new(BigUint::from(0u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ];  // 4x^3 + 1
    
    let poly = Polynomial::new(coefficients);
    assert_eq!(poly.degree(), 3);
    
    let zero_poly = Polynomial::new(vec![
        Fp::new(BigUint::from(0u64), modulus)
    ]);
    assert_eq!(zero_poly.degree(), 0);
}

#[test]
fn test_polynomial_scalar_multiplication() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let poly = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(3u64), modulus.clone()),
    ]);  // 3x^2 + 2x + 1
    
    let scalar = Fp::new(BigUint::from(2u64), modulus.clone());
    let result = &poly * &Polynomial::new(vec![scalar]);
    
    let expected = Polynomial::new(vec![
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
        Fp::new(BigUint::from(6u64), modulus.clone()),
    ]);  // 6x^2 + 4x + 2
    
    assert_eq!(result, expected);
}

#[test]
fn test_polynomial_composition() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let outer = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(1u64), modulus.clone()),
    ]);  // x + 1
    
    let inner = Polynomial::new(vec![
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(1u64), modulus.clone()),
    ]);  // x + 2
    
    let x = Fp::new(BigUint::from(3u64), modulus.clone());
    let inner_result = evaluate_polynomial(&inner, &x);  // g(3)
    
    // Verify that f(g(x)) = f(g(3))
    let composed_result = evaluate_polynomial(&outer, &inner_result);  // f(g(3))
    assert_eq!(composed_result, Fp::new(BigUint::from(6u64), modulus));  // (3 + 2) + 1 = 6
}

#[test]
fn test_polynomial_derivative() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let poly = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(3u64), modulus.clone()),
    ]);  // 3x^2 + 2x + 1
    
    let derivative = poly.derivative();
    // Derivative: 6x + 2
    let expected = Polynomial::new(vec![
        Fp::new(BigUint::from(2u64), modulus.clone()),
        Fp::new(BigUint::from(6u64), modulus),
    ]);
    assert_eq!(derivative, expected);
}

#[test]
fn test_polynomial_interpolation() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let points = vec![
        (
            Fp::new(BigUint::from(0u64), modulus.clone()),
            Fp::new(BigUint::from(1u64), modulus.clone())
        ),
        (
            Fp::new(BigUint::from(1u64), modulus.clone()),
            Fp::new(BigUint::from(2u64), modulus.clone())
        ),
        (
            Fp::new(BigUint::from(2u64), modulus.clone()),
            Fp::new(BigUint::from(4u64), modulus.clone())
        ),
    ];
    
    let poly = Polynomial::interpolate(&points);
    
    // Verify the polynomial passes through all points
    for (x, y) in points {
        let result = evaluate_polynomial(&poly, &x);
        assert_eq!(result, y);
    }
}

#[test]
fn test_polynomial_properties() {
    let modulus = BigUint::from_str_radix(TEST_MODULUS, 10).unwrap();
    let p1 = Polynomial::new(vec![
        Fp::new(BigUint::from(1u64), modulus.clone()),
        Fp::new(BigUint::from(2u64), modulus.clone()),
    ]);  // 2x + 1
    
    let p2 = Polynomial::new(vec![
        Fp::new(BigUint::from(3u64), modulus.clone()),
        Fp::new(BigUint::from(4u64), modulus.clone()),
    ]);  // 4x + 3
    
    let p3 = Polynomial::new(vec![
        Fp::new(BigUint::from(5u64), modulus.clone()),
        Fp::new(BigUint::from(6u64), modulus.clone()),
    ]);  // 6x + 5
    
    // Test associativity: (p1 + p2) + p3 = p1 + (p2 + p3)
    let left = &(&p1 + &p2) + &p3;
    let right = &p1 + &(&p2 + &p3);
    assert_eq!(left, right);
    
    // Test distributivity: p1 * (p2 + p3) = (p1 * p2) + (p1 * p3)
    let left = &p1 * &(&p2 + &p3);
    let right = &(&p1 * &p2) + &(&p1 * &p3);
    assert_eq!(left, right);
} 