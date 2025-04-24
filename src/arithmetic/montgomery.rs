use num_bigint::BigUint;
use std::cmp::{PartialEq, Eq};

/// Constants used for Montgomery arithmetic
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MontgomeryConstants {
    /// The modulus
    pub modulus: BigUint,
    /// R^2 mod N where R = 2^(word_size * num_words)
    pub r_squared: BigUint,
    /// -N^(-1) mod R where R = 2^word_size
    pub n_prime: BigUint,
    /// Word size in bits
    pub word_size: u32,
}

/// Represents a number in Montgomery form with lazy reduction
#[derive(Clone, Debug)]
pub struct MontgomeryForm {
    /// The value in Montgomery representation
    pub value: Vec<u64>,
    /// Number of extra bits of precision (for lazy reduction)
    extra_precision: u32,
    /// Constants for Montgomery arithmetic
    pub constants: MontgomeryConstants,
}

impl PartialEq for MontgomeryForm {
    fn eq(&self, other: &Self) -> bool {
        if self.constants.modulus != other.constants.modulus {
            return false;
        }

        // Create clones and reduce them before comparison
        let mut a = self.clone();
        let mut b = other.clone();
        a.reduce();
        b.reduce();

        // Compare the reduced values
        a.value == b.value
    }
}

impl Eq for MontgomeryForm {}

impl MontgomeryConstants {
    /// Creates new Montgomery constants for the given modulus
    pub fn new(modulus: &BigUint, word_size: u32) -> Self {
        // R = 2^(word_size * 4) for 4 limbs
        let total_bits = word_size * 4;
        let r = BigUint::from(1u64) << total_bits;
        
        // R^2 mod N
        let r_squared = (&r * &r) % modulus;
        
        // -N^(-1) mod R
        let n_prime = calculate_n_prime(modulus, word_size);
        
        Self {
            modulus: modulus.clone(),
            r_squared,
            n_prime,
            word_size,
        }
    }
}

impl MontgomeryForm {
    /// Creates a new value in Montgomery form
    pub fn new(value: Vec<u64>, constants: MontgomeryConstants) -> Self {
        Self {
            value,
            extra_precision: 0,
            constants,
        }
    }

    /// Performs lazy reduction if extra precision is too high
    pub fn reduce(&mut self) {
        // Always reduce if there's any extra precision
        if self.extra_precision > 0 {
            let n = to_limbs(&self.constants.modulus, self.value.len());
            let n_prime = to_limbs(&self.constants.n_prime, self.value.len());
            
            // Convert to Montgomery form with R = 2^(64*4)
            let r_squared = to_limbs(&self.constants.r_squared, self.value.len());
            let mut result = mont_mul(&self.value, &r_squared, &n, &n_prime);
            
            // Final reduction step
            if !ct_lt(&result, &n) {
                let mut borrow = 0i64;
                for i in 0..self.value.len() {
                    let diff = (result[i] as i128) - (n[i] as i128) - (borrow as i128);
                    result[i] = diff as u64;
                    borrow = if diff < 0 { 1 } else { 0 };
                }
            }
            
            self.value = result;
            self.extra_precision = 0;
        }
    }

    /// Multiplies two Montgomery values with lazy reduction
    pub fn mul(&self, other: &Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        let n = to_limbs(&self.constants.modulus, self.value.len());
        let n_prime = to_limbs(&self.constants.n_prime, self.value.len());
        
        let mut result = Self::new(
            mont_mul(&self.value, &other.value, &n, &n_prime),
            self.constants.clone(),
        );
        result.extra_precision = self.extra_precision + other.extra_precision + 1;
        
        // Reduce if precision is too high
        if result.extra_precision >= result.constants.word_size {
            result.reduce();
        }
        
        result
    }
}

/// Performs Montgomery multiplication with lazy reduction
#[inline(always)]
pub fn mont_mul(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    assert_eq!(a.len(), n.len());
    assert_eq!(a.len(), n_prime.len());
    
    let len = a.len();
    let mut t = vec![0u64; len * 2];
    
    // Compute t = a * b with lazy reduction
    for i in 0..len {
        let mut carry = 0u64;
        for j in 0..len {
            let temp = (t[i + j] as u128) + (a[i] as u128 * b[j] as u128) + (carry as u128);
            t[i + j] = temp as u64;
            carry = (temp >> 64) as u64;
        }
        t[i + len] = carry;
    }
    
    // Compute m = t * n' mod R with lazy reduction
    let mut m = vec![0u64; len];
    for i in 0..len {
        let mut carry = 0u64;
        for j in 0..len {
            let temp = (m[j] as u128) + (t[i] as u128 * n_prime[j] as u128) + (carry as u128);
            m[j] = temp as u64;
            carry = (temp >> 64) as u64;
        }
    }
    
    // Compute t = (t + m * n) / R with lazy reduction
    let mut carry = 0u64;
    for i in 0..len {
        let mut sum = carry;
        for j in 0..len {
            let temp = (sum as u128) + (m[j] as u128 * n[i] as u128);
            sum = temp as u64;
            carry = (temp >> 64) as u64;
        }
        t[i] = sum;
    }
    
    // Return result with potential extra bits
    let mut result = t[len..2*len].to_vec();
    
    // Final reduction step if needed
    if !ct_lt(&result, n) {
        let mut borrow = 0i64;
        for i in 0..len {
            let diff = (result[i] as i128) - (n[i] as i128) - (borrow as i128);
            result[i] = diff as u64;
            borrow = if diff < 0 { 1 } else { 0 };
        }
    }
    
    result
}

/// Converts a BigUint to a fixed-size array of limbs
#[inline(always)]
fn to_limbs(value: &BigUint, num_limbs: usize) -> Vec<u64> {
    let mut limbs = vec![0u64; num_limbs];
    let bytes = value.to_bytes_le();
    
    for (i, chunk) in bytes.chunks(8).enumerate() {
        if i >= num_limbs {
            break;
        }
        let mut limb = 0u64;
        for (j, &byte) in chunk.iter().enumerate() {
            limb |= (byte as u64) << (j * 8);
        }
        limbs[i] = limb;
    }
    limbs
}

/// Constant-time less-than comparison
#[inline(always)]
pub fn ct_lt(a: &[u64], b: &[u64]) -> bool {
    let mut lt = false;
    let mut eq = true;
    
    for i in (0..a.len()).rev() {
        lt |= eq && (a[i] < b[i]);
        eq &= a[i] == b[i];
    }
    
    lt
}

/// Calculates -N^(-1) mod R where R = 2^word_size
fn calculate_n_prime(n: &BigUint, word_size: u32) -> BigUint {
    let r = BigUint::from(1u64) << word_size;
    let mut t = BigUint::from(1u64);
    
    // Find n_prime such that: n * n_prime ≡ -1 (mod R)
    // This means: n * n_prime + k * R = -1 for some k
    for _ in 0..word_size {
        if (&t * n) % &r == (r.clone() - BigUint::from(1u64)) {
            return t;
        }
        t = t << 1;
    }
    
    panic!("Failed to find Montgomery inverse");
} 