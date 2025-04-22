use num_bigint::BigUint;

/// Constants used for Montgomery arithmetic
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MontgomeryConstants {
    /// The modulus
    pub modulus: BigUint,
    /// R^2 mod N where R = 2^(word_size * num_words)
    pub r_squared: BigUint,
    /// -N^(-1) mod R where R = 2^word_size
    pub n_prime: BigUint,
}

impl MontgomeryConstants {
    /// Creates new Montgomery constants for the given modulus
    pub fn new(modulus: &BigUint, word_size: u32) -> Self {
        let r = BigUint::from(1u64) << word_size;
        let r_squared = (&r * &r) % modulus;
        let n_prime = calculate_n_prime(modulus, word_size);
        
        Self {
            modulus: modulus.clone(),
            r_squared,
            n_prime,
        }
    }
}

/// Performs Montgomery multiplication: (a * b * R^(-1)) mod N
#[inline(always)]
pub fn mont_mul(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    let len = a.len();
    let mut t = vec![0u64; len * 2];
    
    // Compute t = a * b
    for i in 0..len {
        let mut carry = 0u64;
        for j in 0..len {
            let temp = (t[i + j] as u128) + (a[i] as u128 * b[j] as u128) + (carry as u128);
            t[i + j] = temp as u64;
            carry = (temp >> 64) as u64;
        }
        t[i + len] = carry;
    }
    
    // Compute m = t * n' mod R
    let mut m = vec![0u64; len];
    for i in 0..len {
        let mut carry = 0u64;
        for j in 0..len {
            let temp = (m[j] as u128) + (t[i] as u128 * n_prime[j] as u128) + (carry as u128);
            m[j] = temp as u64;
            carry = (temp >> 64) as u64;
        }
    }
    
    // Compute t = (t + m * n) / R
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
    
    // Final reduction
    let mut result = t[len..2*len].to_vec();
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
    let mut n_prime = BigUint::from(0u64);
    
    for _ in 0..word_size {
        if (&t * n) % &r == BigUint::from(1u64) {
            n_prime = t.clone();
            break;
        }
        t = t << 1;
    }
    
    n_prime
} 