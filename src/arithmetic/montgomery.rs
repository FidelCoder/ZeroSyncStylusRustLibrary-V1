use std::sync::OnceLock;
use num_bigint::BigUint;
use std::ops::Shl;
use crate::arithmetic::traits::Field;

/// Thread-local storage for Montgomery constants
thread_local! {
    static MONTGOMERY_CONSTANTS: OnceLock<MontgomeryConstants> = OnceLock::new();
}

#[derive(Clone, Debug)]
pub(crate) struct MontgomeryConstants {
    r: BigUint,       // R = 2^word_size mod N
    r_squared: BigUint, // R^2 mod N
    n_prime: BigUint,   // -N^(-1) mod R
}

impl MontgomeryConstants {
    pub fn new(modulus: &BigUint, word_size: u32) -> Self {
        let r = BigUint::from(1u64).shl(word_size);
        let r_squared = (&r * &r) % modulus;
        
        // Calculate n_prime using extended GCD
        let n_prime = Self::compute_n_prime(modulus, &r);
        
        Self {
            r,
            r_squared,
            n_prime,
        }
    }

    /// Compute -N^(-1) mod R using extended GCD
    fn compute_n_prime(n: &BigUint, r: &BigUint) -> BigUint {
        let (mut t, mut new_t) = (BigUint::from(0u32), BigUint::from(1u32));
        let (mut r_copy, mut new_r) = (r.clone(), n.clone());

        while !new_r.is_zero() {
            let quotient = &r_copy / &new_r;
            (t, new_t) = (new_t.clone(), t - quotient.clone() * new_t);
            (r_copy, new_r) = (new_r.clone(), r_copy - quotient * new_r);
        }

        if t < BigUint::from(0u32) {
            t = t + r;
        }
        r - t
    }
}

/// Optimized Montgomery multiplication
#[inline(always)]
pub fn mont_mul(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    let len = a.len();
    let mut t = vec![0u64; len * 2];
    
    // Step 1: Compute t = a * b
    for i in 0..len {
        let mut carry = 0u64;
        for j in 0..len {
            let prod = (t[i + j] as u128)
                    + (a[i] as u128 * b[j] as u128)
                    + (carry as u128);
            t[i + j] = prod as u64;
            carry = (prod >> 64) as u64;
        }
        t[i + len] = carry;
    }
    
    // Step 2: Compute m = t * n_prime mod 2^64
    let mut m = vec![0u64; len];
    for i in 0..len {
        let mut carry = 0u64;
        let mu = (t[i] as u128 * n_prime[0] as u128) as u64;
        for j in 0..len {
            let prod = (t[i + j] as u128)
                    + (mu as u128 * n[j] as u128)
                    + (carry as u128);
            t[i + j] = prod as u64;
            carry = (prod >> 64) as u64;
        }
        t[i + len] += carry;
    }
    
    // Step 3: Divide by R (shift right by word_size)
    let result: Vec<u64> = t[len..].to_vec();
    
    // Step 4: Final reduction
    if result.iter().zip(n.iter()).rev()
        .find(|(&a, &b)| a != b)
        .map_or(false, |(a, b)| a >= b) {
        let mut borrow = 0i64;
        let mut final_result = vec![0u64; len];
        for i in 0..len {
            let diff = (result[i] as i128) - (n[i] as i128) - (borrow as i128);
            final_result[i] = diff as u64;
            borrow = if diff < 0 { 1 } else { 0 };
        }
        final_result
    } else {
        result
    }
}

/// Constant-time comparison
#[inline(always)]
pub fn ct_lt(a: &[u64], b: &[u64]) -> bool {
    let mut result = 0u8;
    let mut equal = 1u8;
    
    for (x, y) in a.iter().zip(b.iter()).rev() {
        let lt = ((x < y) as u8) & equal;
        let gt = ((x > y) as u8) & equal;
        result |= lt;
        equal &= !lt & !gt;
    }
    
    result == 1
} 