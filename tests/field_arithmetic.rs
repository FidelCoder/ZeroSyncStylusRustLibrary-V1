use num_bigint::BigUint;
use num_traits::{One, Zero};
use zerosync::arithmetic::{
    field::Fp,
    montgomery::{mont_mul, ct_lt},
};
use proptest::prelude::*;

// BN254 base field modulus
const BN254_MODULUS: &str = "21888242871839275222246405745257275088696311157297823662689037894645226208583";

#[test]
fn test_field_basic_arithmetic() {
    let modulus = BigUint::from_str_radix(BN254_MODULUS, 10).unwrap();
    let a = Fp::new(BigUint::from(5u32), modulus.clone());
    let b = Fp::new(BigUint::from(3u32), modulus.clone());
    let c = Fp::new(BigUint::from(7u32), modulus.clone());

    // Test addition
    let sum = a.clone() + b.clone();
    assert_eq!(sum, Fp::new(BigUint::from(8u32), modulus.clone()));

    // Test subtraction
    let diff = a.clone() - b.clone();
    assert_eq!(diff, Fp::new(BigUint::from(2u32), modulus.clone()));

    // Test multiplication
    let prod = a.clone() * b.clone();
    assert_eq!(prod, Fp::new(BigUint::from(15u32), modulus.clone()));

    // Test associativity: (a + b) + c = a + (b + c)
    let left = (a.clone() + b.clone()) + c.clone();
    let right = a.clone() + (b.clone() + c.clone());
    assert_eq!(left, right);

    // Test commutativity: a * b = b * a
    let prod1 = a.clone() * b.clone();
    let prod2 = b.clone() * a.clone();
    assert_eq!(prod1, prod2);
}

#[test]
fn test_field_edge_cases() {
    let modulus = BigUint::from_str_radix(BN254_MODULUS, 10).unwrap();
    
    // Test zero
    let zero = Fp::new(BigUint::zero(), modulus.clone());
    let a = Fp::new(BigUint::from(5u32), modulus.clone());
    
    assert_eq!(a.clone() + zero.clone(), a.clone());
    assert_eq!(a.clone() * zero.clone(), zero.clone());
    
    // Test one
    let one = Fp::new(BigUint::one(), modulus.clone());
    assert_eq!(a.clone() * one.clone(), a.clone());
    
    // Test modulus reduction
    let large = Fp::new(modulus.clone(), modulus.clone());
    assert_eq!(large, zero);
}

#[test]
fn test_montgomery_form() {
    let modulus = BigUint::from_str_radix(BN254_MODULUS, 10).unwrap();
    let a = Fp::new(BigUint::from(5u32), modulus.clone());
    
    // Test Montgomery multiplication identity
    let one = Fp::new(BigUint::one(), modulus.clone());
    let result = a.clone() * one;
    assert_eq!(result, a);
}

#[test]
fn test_field_exponentiation() {
    let modulus = BigUint::from_str_radix(BN254_MODULUS, 10).unwrap();
    let base = Fp::new(BigUint::from(2u32), modulus.clone());
    
    // Test small exponents
    assert_eq!(base.pow(0), Fp::new(BigUint::one(), modulus.clone()));
    assert_eq!(base.pow(1), base.clone());
    assert_eq!(base.pow(2), base.clone() * base.clone());
    
    // Test larger exponent
    let expected = Fp::new(BigUint::from(16u32), modulus.clone());
    assert_eq!(base.pow(4), expected);
}

proptest! {
    #[test]
    fn test_field_properties(
        a in 0u32..1000u32,
        b in 0u32..1000u32,
        c in 0u32..1000u32
    ) {
        let modulus = BigUint::from_str_radix(BN254_MODULUS, 10).unwrap();
        let fa = Fp::new(BigUint::from(a), modulus.clone());
        let fb = Fp::new(BigUint::from(b), modulus.clone());
        let fc = Fp::new(BigUint::from(c), modulus.clone());

        // Associativity
        prop_assert_eq!((fa.clone() + fb.clone()) + fc.clone(), fa.clone() + (fb.clone() + fc.clone()));
        prop_assert_eq!((fa.clone() * fb.clone()) * fc.clone(), fa.clone() * (fb.clone() * fc.clone()));

        // Commutativity
        prop_assert_eq!(fa.clone() + fb.clone(), fb.clone() + fa.clone());
        prop_assert_eq!(fa.clone() * fb.clone(), fb.clone() * fa.clone());

        // Distributivity
        prop_assert_eq!(
            fa.clone() * (fb.clone() + fc.clone()),
            (fa.clone() * fb.clone()) + (fa.clone() * fc.clone())
        );
    }
}

#[test]
fn test_constant_time_comparison() {
    let a = [1u64, 2, 3, 4];
    let b = [1u64, 2, 3, 5];
    let c = [1u64, 2, 3, 4];

    assert!(ct_lt(&a, &b));
    assert!(!ct_lt(&b, &a));
    assert!(!ct_lt(&a, &c));
}

#[test]
fn test_montgomery_multiplication() {
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];
    let a = [1u64, 0, 0, 0];
    let b = [2u64, 0, 0, 0];
    let n_prime = [1u64, 0, 0, 0];  // Simplified for test

    let result = mont_mul(&a, &b, &modulus, &n_prime);
    assert_eq!(result.len(), 4);
    // Add more specific assertions based on expected results
} 