use num_bigint::BigUint;
use std::cmp::{PartialEq, Eq};
use num_traits::{One, Zero};

#[cfg(feature = "simd")]
use super::simd_avx512;

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
    /// Number of extra bits available for lazy reduction
    pub extra_bits: u32,
}

/// Represents a number in Montgomery form with lazy reduction
#[derive(Clone, Debug)]
pub struct MontgomeryForm {
    /// The value in Montgomery representation
    pub value: Vec<u64>,
    /// Number of extra bits of precision (for lazy reduction)
    pub extra_precision: u32,
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
            extra_bits: 64, // Default to 64 extra bits for lazy reduction
        }
    }
    
    /// Creates new Montgomery constants with specified extra bits for lazy reduction
    pub fn new_with_lazy_reduction(modulus: &BigUint, word_size: u32, extra_bits: u32) -> Self {
        let mut constants = Self::new(modulus, word_size);
        constants.extra_bits = extra_bits;
        constants
    }
}

impl MontgomeryForm {
    /// Creates a new value in Montgomery form
    pub fn new(mut value: Vec<u64>, constants: MontgomeryConstants) -> Self {
        // Ensure the value is properly reduced before conversion
        let modulus_limbs = to_limbs(&constants.modulus, value.len());
        if !ct_lt(&value, &modulus_limbs) {
            let mut borrow = 0i64;
            for i in 0..value.len() {
                let diff = (value[i] as i128) - (modulus_limbs[i] as i128) - (borrow as i128);
                value[i] = diff as u64;
                borrow = if diff < 0 { 1 } else { 0 };
            }
        }
        
        // Convert to Montgomery form by multiplying by R^2 mod N
        let r_squared_limbs = to_limbs(&constants.r_squared, value.len());
        let n_prime_limbs = to_limbs(&constants.n_prime, value.len());
        
        let result = mont_mul(
            &value,
            &r_squared_limbs,
            &modulus_limbs,
            &n_prime_limbs
        );
        
        Self {
            value: result,
            extra_precision: 0,
            constants,
        }
    }

    /// Performs lazy reduction if extra precision is too high
    pub fn reduce(&mut self) {
        // Only reduce if the extra precision is beyond our threshold
        if self.extra_precision > self.constants.extra_bits {
            let n = to_limbs(&self.constants.modulus, self.value.len());
            let n_prime = to_limbs(&self.constants.n_prime, self.value.len());
            
            // Use optimized SIMD implementation when available
            #[cfg(all(feature = "simd", target_arch = "x86_64"))]
            if self.value.len() == 8 && simd_avx512::has_avx512f() {
                unsafe {
                    self.value = simd_avx512::mont_reduce_avx512(&self.value, &n, &n_prime);
                    self.extra_precision = 0;
                    return;
                }
            }
            
            // Fallback to standard Montgomery reduction
            self.value = mont_reduce(&self.value, &n, &n_prime);
            self.extra_precision = 0;
        }
    }

    /// Multiplies two Montgomery values with lazy reduction
    pub fn mul(&mut self, other: &Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        
        #[cfg(all(feature = "simd", target_arch = "x86_64"))]
        if self.value.len() == 8 && simd_avx512::has_avx512f() {
            unsafe {
                let result = simd_avx512::mont_mul_avx512(
                    &self.value,
                    &other.value,
                    &to_limbs(&self.constants.modulus, self.value.len()),
                    &to_limbs(&self.constants.n_prime, self.value.len())
                );
                
                // Track extra precision for lazy reduction
                return Self {
                    value: result,
                    constants: self.constants.clone(),
                    extra_precision: self.extra_precision + other.extra_precision + 1
                };
            }
        }
        
        // Fallback to non-SIMD implementation with lazy reduction
        let result = mont_mul_lazy(
            &self.value,
            &other.value,
            &to_limbs(&self.constants.modulus, self.value.len()),
            &to_limbs(&self.constants.n_prime, self.value.len())
        );
        
        Self {
            value: result,
            constants: self.constants.clone(),
            extra_precision: self.extra_precision + other.extra_precision + 1
        }
    }

    pub fn add(&mut self, other: &Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        
        #[cfg(all(feature = "simd", target_arch = "x86_64"))]
        if self.value.len() == 8 && simd_avx512::has_avx512f() {
            unsafe {
                let result = simd_avx512::field_add_avx512(
                    &self.value,
                    &other.value,
                    &to_limbs(&self.constants.modulus, self.value.len())
                );
                
                // Track extra precision
                return Self {
                    value: result,
                    constants: self.constants.clone(),
                    extra_precision: self.extra_precision.max(other.extra_precision)
                };
            }
        }
        
        // Perform addition with lazy reduction
        let mut result = self.value.clone();
        let modulus_limbs = to_limbs(&self.constants.modulus, self.value.len());
        
        let mut carry = 0u64;
        for i in 0..result.len() {
            let sum = (result[i] as u128) + (other.value[i] as u128) + (carry as u128);
            result[i] = sum as u64;
            carry = (sum >> 64) as u64;
        }
        
        // Only reduce if value exceeds modulus by a threshold
        let needs_reduction = carry > 0 || !ct_lt(&result, &modulus_limbs);
        if needs_reduction {
            let mut borrow = 0i64;
            for i in 0..result.len() {
                let diff = (result[i] as i128) - (modulus_limbs[i] as i128) - (borrow as i128);
                result[i] = diff as u64;
                borrow = if diff < 0 { 1 } else { 0 };
            }
        }
        
        // Track extra precision - addition shouldn't increase precision much
        let new_precision = self.extra_precision.max(other.extra_precision);
        
        Self {
            value: result,
            constants: self.constants.clone(),
            extra_precision: new_precision
        }
    }
    
    /// Subtracts another Montgomery value with lazy reduction
    pub fn sub(&mut self, other: &Self) -> Self {
        assert_eq!(self.constants.modulus, other.constants.modulus);
        
        #[cfg(all(feature = "simd", target_arch = "x86_64"))]
        if self.value.len() == 8 && simd_avx512::has_avx512f() {
            unsafe {
                let result = simd_avx512::field_sub_avx512(
                    &self.value,
                    &other.value,
                    &to_limbs(&self.constants.modulus, self.value.len())
                );
                
                // Track extra precision
                return Self {
                    value: result,
                    constants: self.constants.clone(),
                    extra_precision: self.extra_precision.max(other.extra_precision)
                };
            }
        }
        
        // Perform subtraction with lazy reduction
        let mut result = self.value.clone();
        let modulus_limbs = to_limbs(&self.constants.modulus, self.value.len());
        
        // If self < other, add the modulus first
        if ct_lt(&self.value, &other.value) {
            let mut carry = 0u64;
            for i in 0..result.len() {
                let sum = (result[i] as u128) + (modulus_limbs[i] as u128) + (carry as u128);
                result[i] = sum as u64;
                carry = (sum >> 64) as u64;
            }
        }
        
        // Now perform the subtraction
        let mut borrow = 0i64;
        for i in 0..result.len() {
            let diff = (result[i] as i128) - (other.value[i] as i128) - (borrow as i128);
            result[i] = diff as u64;
            borrow = if diff < 0 { 1 } else { 0 };
        }
        
        // If we still have a borrow, add the modulus again
        if borrow != 0 {
            let mut carry = 0u64;
            for i in 0..result.len() {
                let sum = (result[i] as u128) + (modulus_limbs[i] as u128) + (carry as u128);
                result[i] = sum as u64;
                carry = (sum >> 64) as u64;
            }
        }
        
        Self {
            value: result,
            constants: self.constants.clone(),
            extra_precision: self.extra_precision.max(other.extra_precision)
        }
    }
    
    /// Squares the value with lazy reduction (optimized version of mul)
    pub fn square(&mut self) -> Self {
        #[cfg(all(feature = "simd", target_arch = "x86_64"))]
        if self.value.len() == 8 && simd_avx512::has_avx512f() {
            // For squaring, we can use the same mul routine
            return self.mul(self);
        }
        
        // Specialized squaring algorithm with lazy reduction
        let modulus_limbs = to_limbs(&self.constants.modulus, self.value.len());
        let n_prime_limbs = to_limbs(&self.constants.n_prime, self.value.len());
        
        // Compute efficient squaring with fewer multiplications
        let mut t = vec![0u64; self.value.len() * 2];
        
        // Diagonal terms
        for i in 0..self.value.len() {
            let sq = (self.value[i] as u128) * (self.value[i] as u128);
            t[i * 2] += sq as u64;
            t[i * 2 + 1] += (sq >> 64) as u64;
        }
        
        // Off-diagonal terms (multiplied by 2)
        for i in 0..self.value.len() {
            for j in i+1..self.value.len() {
                let product = (self.value[i] as u128) * (self.value[j] as u128);
                let lo = (product as u64).wrapping_mul(2);
                let hi = ((product >> 64) as u64).wrapping_mul(2) 
                       + (if lo < (product as u64) { 1 } else { 0 });
                
                // Add to result with carrying
                let mut pos = i + j;
                let mut carry = 0u64;
                
                let sum = (t[pos] as u128) + (lo as u128) + (carry as u128);
                t[pos] = sum as u64;
                carry = (sum >> 64) as u64;
                
                pos += 1;
                let sum = (t[pos] as u128) + (hi as u128) + (carry as u128);
                t[pos] = sum as u64;
                carry = (sum >> 64) as u64;
                
                // Propagate the carry
                pos += 1;
                while carry > 0 && pos < t.len() {
                    let sum = (t[pos] as u128) + (carry as u128);
                    t[pos] = sum as u64;
                    carry = (sum >> 64) as u64;
                    pos += 1;
                }
            }
        }
        
        // Now perform Montgomery reduction (can be lazy)
        let result = mont_reduce(&t, &modulus_limbs, &n_prime_limbs);
        
        Self {
            value: result,
            constants: self.constants.clone(),
            extra_precision: self.extra_precision * 2 + 1
        }
    }
}

/// Montgomery multiplication with optimized implementation
pub fn mont_mul(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    // Compute t = a * b
    let mut t = vec![0u64; a.len() * 2];
    for i in 0..a.len() {
        let mut carry = 0u64;
        for j in 0..b.len() {
            let product = (a[i] as u128) * (b[j] as u128) + (t[i + j] as u128) + (carry as u128);
            t[i + j] = product as u64;
            carry = (product >> 64) as u64;
        }
        t[i + b.len()] = carry;
    }
    
    // Compute m = (t mod R) * n' mod R
    let mut m = vec![0u64; a.len()];
    for i in 0..a.len() {
        let mu = ((t[i] as u128) * (n_prime[0] as u128)) as u64;
        m[i] = mu;
    }
    
    // Compute t = (t + m*n) / R
    let mut carry = 0u64;
    for i in 0..a.len() {
        let mut carry2 = 0u64;
        for j in 0..n.len() {
            let product = (m[i] as u128) * (n[j] as u128) + (t[i + j] as u128) + (carry2 as u128);
            t[i + j] = product as u64;
            carry2 = (product >> 64) as u64;
        }
        let sum = (t[i + n.len()] as u128) + (carry2 as u128) + (carry as u128);
        t[i + n.len()] = sum as u64;
        carry = (sum >> 64) as u64;
    }
    
    // Extract higher limbs as the result
    let mut result = vec![0u64; a.len()];
    for i in 0..a.len() {
        result[i] = t[i + a.len()];
    }
    
    // Final reduction step
    if !ct_lt(&result, n) {
        let mut borrow = 0i64;
        for i in 0..result.len() {
            let diff = (result[i] as i128) - (n[i] as i128) - (borrow as i128);
            result[i] = diff as u64;
            borrow = if diff < 0 { 1 } else { 0 };
        }
    }
    
    result
}

/// Montgomery multiplication with lazy reduction for improved performance
pub fn mont_mul_lazy(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    // Specialized version that skips the final reduction step
    let mut t = vec![0u64; a.len() * 2];
    
    // Step 1: Compute t = a * b
    for i in 0..a.len() {
        let mut carry = 0u64;
        for j in 0..b.len() {
            let product = (a[i] as u128) * (b[j] as u128) + (t[i + j] as u128) + (carry as u128);
            t[i + j] = product as u64;
            carry = (product >> 64) as u64;
        }
        t[i + b.len()] = carry;
    }
    
    // Step 2: Compute m = (t mod R) * n' mod R
    let mut m = vec![0u64; a.len()];
    for i in 0..a.len() {
        let mu = ((t[i] as u128) * (n_prime[0] as u128)) as u64;
        m[i] = mu;
    }
    
    // Step 3: Compute t = (t + m*n) / R
    let mut carry = 0u64;
    for i in 0..a.len() {
        let mut carry2 = 0u64;
        for j in 0..n.len() {
            let product = (m[i] as u128) * (n[j] as u128) + (t[i + j] as u128) + (carry2 as u128);
            t[i + j] = product as u64;
            carry2 = (product >> 64) as u64;
        }
        let sum = (t[i + n.len()] as u128) + (carry2 as u128) + (carry as u128);
        t[i + n.len()] = sum as u64;
        carry = (sum >> 64) as u64;
    }
    
    // Extract higher limbs as the result
    let mut result = vec![0u64; a.len()];
    for i in 0..a.len() {
        result[i] = t[i + a.len()];
    }
    
    // Skip final reduction step for lazy reduction
    result
}

/// Montgomery reduction from double-precision to single-precision
pub fn mont_reduce(t: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    let num_limbs = n.len();
    let mut result = vec![0u64; num_limbs];
    
    // Extract lower half of t
    let mut t_lo = vec![0u64; num_limbs];
    for i in 0..num_limbs {
        t_lo[i] = t[i];
    }
    
    // Compute m = t_lo * n' mod R
    let mut m = vec![0u64; num_limbs];
    for i in 0..num_limbs {
        let mu = ((t_lo[i] as u128) * (n_prime[0] as u128)) as u64;
        m[i] = mu;
    }
    
    // Compute t = (t + m*n) / R
    let mut carry = 0u64;
    for i in 0..num_limbs {
        let mut carry2 = 0u64;
        for j in 0..num_limbs {
            let idx = i + j;
            if idx < t.len() {
                let product = (m[i] as u128) * (n[j] as u128) + (t[idx] as u128) + (carry2 as u128);
                t_lo[idx % num_limbs] = product as u64;
                carry2 = (product >> 64) as u64;
            }
        }
        carry += carry2;
    }
    
    // Extract result from higher half of t
    for i in 0..num_limbs {
        let idx = i + num_limbs;
        result[i] = if idx < t.len() { t[idx] } else { 0 };
    }
    
    // Add carry from the lower half computation
    let mut c = 0u64;
    for i in 0..num_limbs {
        let sum = (result[i] as u128) + (c as u128);
        result[i] = sum as u64;
        c = (sum >> 64) as u64;
    }
    
    // Final reduction step
    if !ct_lt(&result, n) {
        let mut borrow = 0i64;
        for i in 0..result.len() {
            let diff = (result[i] as i128) - (n[i] as i128) - (borrow as i128);
            result[i] = diff as u64;
            borrow = if diff < 0 { 1 } else { 0 };
        }
    }
    
    result
}

/// Converts a BigUint to a fixed-size array of limbs
pub fn to_limbs(value: &BigUint, num_limbs: usize) -> Vec<u64> {
    let bytes = value.to_bytes_le();
    let mut limbs = vec![0u64; num_limbs];
    
    for i in 0..num_limbs {
        let start = i * 8;
        if start >= bytes.len() {
            break;
        }
        
        let end = std::cmp::min(start + 8, bytes.len());
        let mut limb = 0u64;
        
        for j in start..end {
            limb |= (bytes[j] as u64) << (8 * (j - start));
        }
        
        limbs[i] = limb;
    }
    
    limbs
}

/// Converts a slice of u64 limbs to bytes in little-endian order
pub fn to_bytes(limbs: &[u64]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(limbs.len() * 8);
    for &limb in limbs {
        for j in 0..8 {
            bytes.push(((limb >> (j * 8)) & 0xFF) as u8);
        }
    }
    bytes
}

/// Constant-time less-than comparison
pub fn ct_lt(a: &[u64], b: &[u64]) -> bool {
    assert_eq!(a.len(), b.len());
    
    let mut result = 0u8;
    let mut cmp = 1u8;
    
    for i in (0..a.len()).rev() {
        let eq = (a[i] == b[i]) as u8;
        let lt = (a[i] < b[i]) as u8;
        
        result |= cmp & lt;
        cmp &= eq;
    }
    
    result == 1
}

/// Calculates n' such that n * n' ≡ -1 (mod 2^64)
fn calculate_n_prime(n: &BigUint, word_size: u32) -> BigUint {
    let r = BigUint::from(1u64) << word_size;
    let one = BigUint::from(1u64);
    
    // n' = -n^(-1) mod 2^64
    // First, calculate n mod 2^64
    let n_mod_r = n % &r;
    
    // Extended Euclidean algorithm to find inverse of n modulo 2^64
    let mut t = BigUint::from(0u64);
    let mut new_t = BigUint::from(1u64);
    let mut r_val = r.clone();
    let mut new_r = n_mod_r.clone();
    
    while !new_r.is_zero() {
        let quotient = &r_val / &new_r;
        let temp_t = t.clone();
        t = new_t.clone();
        new_t = temp_t - &quotient * &new_t;
        
        let temp_r = r_val.clone();
        r_val = new_r.clone();
        new_r = temp_r - &quotient * &new_r;
    }
    
    // If r_val > 1, n is not invertible modulo 2^64
    assert_eq!(r_val, one, "Modulus is not invertible modulo 2^64");
    
    // Ensure t is positive
    while t < BigUint::from(0u64) {
        t = t + &r;
    }
    
    // Compute n' = -t mod 2^64
    (&r - &t) % &r
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use std::str::FromStr;
    
    #[test]
    fn test_basic_conversion() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        let constants = MontgomeryConstants::new(&modulus, 64);
        
        // Convert to Montgomery form and back
        let value = BigUint::from(5u32);
        let mont_form = MontgomeryForm::new(to_limbs(&value, 4), constants.clone());
        
        // Convert back to regular form
        let one_limbs = vec![1u64; 4];
        let n_limbs = to_limbs(&modulus, 4);
        let n_prime_limbs = to_limbs(&constants.n_prime, 4);
        
        let result_limbs = mont_mul(&mont_form.value, &one_limbs, &n_limbs, &n_prime_limbs);
        let result = BigUint::from_bytes_le(&to_bytes(&result_limbs));
        
        // Should get back original value
        assert_eq!(result, value);
    }
    
    #[test]
    fn test_addition() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        let constants = MontgomeryConstants::new(&modulus, 64);
        
        // Create Montgomery forms
        let a_val = BigUint::from(5u32);
        let b_val = BigUint::from(7u32);
        let expected = (a_val.clone() + b_val.clone()) % &modulus; // 5 + 7 = 12
        
        let mut a_mont = MontgomeryForm::new(to_limbs(&a_val, 4), constants.clone());
        let b_mont = MontgomeryForm::new(to_limbs(&b_val, 4), constants.clone());
        
        // Add in Montgomery form
        let result = a_mont.add(&b_mont);
        
        // Convert back to regular form
        let one_limbs = vec![1u64; 4];
        let n_limbs = to_limbs(&modulus, 4);
        let n_prime_limbs = to_limbs(&constants.n_prime, 4);
        
        let result_limbs = mont_mul(&result.value, &one_limbs, &n_limbs, &n_prime_limbs);
        let result_val = BigUint::from_bytes_le(&to_bytes(&result_limbs));
        
        // Verify result
        assert_eq!(result_val, expected);
    }
    
    #[test]
    fn test_subtraction() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        let constants = MontgomeryConstants::new(&modulus, 64);
        
        // Test case 1: a > b
        let a_val = BigUint::from(12u32);
        let b_val = BigUint::from(7u32);
        let expected = (a_val.clone() - b_val.clone()) % &modulus; // 12 - 7 = 5
        
        let mut a_mont = MontgomeryForm::new(to_limbs(&a_val, 4), constants.clone());
        let b_mont = MontgomeryForm::new(to_limbs(&b_val, 4), constants.clone());
        
        // Subtract in Montgomery form
        let result = a_mont.sub(&b_mont);
        
        // Convert back to regular form
        let one_limbs = vec![1u64; 4];
        let n_limbs = to_limbs(&modulus, 4);
        let n_prime_limbs = to_limbs(&constants.n_prime, 4);
        
        let result_limbs = mont_mul(&result.value, &one_limbs, &n_limbs, &n_prime_limbs);
        let result_val = BigUint::from_bytes_le(&to_bytes(&result_limbs));
        
        // Verify result
        assert_eq!(result_val, expected);
        
        // Test case 2: a < b (modular subtraction)
        let a_val = BigUint::from(5u32);
        let b_val = BigUint::from(10u32);
        let expected = ((&modulus - &b_val) + &a_val) % &modulus; // 5 - 10 = -5 = 12 mod 17
        
        let mut a_mont = MontgomeryForm::new(to_limbs(&a_val, 4), constants.clone());
        let b_mont = MontgomeryForm::new(to_limbs(&b_val, 4), constants.clone());
        
        // Subtract in Montgomery form
        let result = a_mont.sub(&b_mont);
        
        // Convert back to regular form
        let result_limbs = mont_mul(&result.value, &one_limbs, &n_limbs, &n_prime_limbs);
        let result_val = BigUint::from_bytes_le(&to_bytes(&result_limbs));
        
        // Verify result
        assert_eq!(result_val, expected);
    }
    
    #[test]
    fn test_multiplication() {
        // Use a small prime field for testing
        let modulus = BigUint::from(17u32);
        let constants = MontgomeryConstants::new(&modulus, 64);
        
        // Create Montgomery forms
        let a_val = BigUint::from(5u32);
        let b_val = BigUint::from(7u32);
        let expected = (a_val.clone() * b_val.clone()) % &modulus; // 5 * 7 = 35 = 1 mod 17
        
        let mut a_mont = MontgomeryForm::new(to_limbs(&a_val, 4), constants.clone());
        let b_mont = MontgomeryForm::new(to_limbs(&b_val, 4), constants.clone());
        
        // Multiply in Montgomery form
        let result = a_mont.mul(&b_mont);
        
        // Convert back to regular form
        let one_limbs = vec![1u64; 4];
        let n_limbs = to_limbs(&modulus, 4);
        let n_prime_limbs = to_limbs(&constants.n_prime, 4);
        
        let result_limbs = mont_mul(&result.value, &one_limbs, &n_limbs, &n_prime_limbs);
        let result_val = BigUint::from_bytes_le(&to_bytes(&result_limbs));
        
        // Verify result
        assert_eq!(result_val, expected);
    }
} 