#![cfg(feature = "simd")]
#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::arithmetic::field::Fp;
use std::mem;
use num_bigint::BigUint;

/// Returns true if AVX2 instructions are available
#[inline]
pub fn has_avx2() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        is_x86_feature_detected!("avx2")
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

/// Performs field multiplication using AVX2 instructions
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub unsafe fn field_mul_avx2(a: &[u64], b: &[u64], modulus: &[u64]) -> Vec<u64> {
    assert_eq!(a.len(), 4);
    assert_eq!(b.len(), 4);
    assert_eq!(modulus.len(), 4);

    let mut result = vec![0u64; 4];

    // Load vectors
    let a_vec = _mm256_loadu_si256(a.as_ptr() as *const __m256i);
    let b_vec = _mm256_loadu_si256(b.as_ptr() as *const __m256i);
    let m_vec = _mm256_loadu_si256(modulus.as_ptr() as *const __m256i);

    // Multiply and reduce modulo m
    let mut temp = _mm256_mul_epu32(a_vec, b_vec);
    temp = _mm256_rem_epu32(temp, m_vec);

    // Store result
    _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, temp);
    result
}

/// Performs field addition using AVX2 instructions
#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
pub unsafe fn field_add_avx2(a: &[u64], b: &[u64], modulus: &[u64]) -> Vec<u64> {
    assert_eq!(a.len(), 4);
    assert_eq!(b.len(), 4);
    assert_eq!(modulus.len(), 4);

    let mut result = vec![0u64; 4];

    // Load vectors
    let a_vec = _mm256_loadu_si256(a.as_ptr() as *const __m256i);
    let b_vec = _mm256_loadu_si256(b.as_ptr() as *const __m256i);
    let m_vec = _mm256_loadu_si256(modulus.as_ptr() as *const __m256i);

    // Add and reduce modulo m
    let mut temp = _mm256_add_epi64(a_vec, b_vec);
    temp = _mm256_rem_epu32(temp, m_vec);

    // Store result
    _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, temp);
    result
}

// Fallback implementations for non-x86_64 or non-AVX2 platforms
#[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
pub fn field_mul_avx2(_a: &[u64], _b: &[u64], _modulus: &[u64]) -> Vec<u64> {
    panic!("AVX2 not available");
}

#[cfg(not(all(target_arch = "x86_64", target_feature = "avx2")))]
pub fn field_add_avx2(_a: &[u64], _b: &[u64], _modulus: &[u64]) -> Vec<u64> {
    panic!("AVX2 not available");
}

/// Performs 4-way parallel Montgomery reduction using AVX2
#[target_feature(enable = "avx2")]
pub unsafe fn mont_reduce_avx2(t: &[u64; 8], modulus: &[u64; 4], n_prime: u64) -> [u64; 4] {
    let mut result = [0u64; 4];
    
    // Load lower and upper halves
    let t_lo = _mm256_loadu_si256(t.as_ptr() as *const __m256i);
    let t_hi = _mm256_loadu_si256(t[4..].as_ptr() as *const __m256i);
    let modulus_vec = _mm256_loadu_si256(modulus.as_ptr() as *const __m256i);
    let n_prime_vec = _mm256_set1_epi64x(n_prime as i64);
    
    // Compute m = t * n_prime mod 2^64
    let m = _mm256_mul_epu32(t_lo, n_prime_vec);
    
    // Compute t + m * N
    let mn = _mm256_mul_epu32(m, modulus_vec);
    let sum = _mm256_add_epi64(t_hi, mn);
    
    // Final reduction if needed
    let mask = _mm256_cmpgt_epi64(sum, modulus_vec);
    let reduced = _mm256_sub_epi64(sum, _mm256_and_si256(mask, modulus_vec));
    
    _mm256_storeu_si256(result.as_mut_ptr() as *mut __m256i, reduced);
    result
}

/// SIMD-optimized field arithmetic operations
pub struct SimdFieldOps {
    /// Field modulus
    modulus: BigUint,
    /// Precomputed Montgomery constants
    r_squared: BigUint,
    /// Precomputed inverse of modulus
    inv_modulus: BigUint,
}

impl SimdFieldOps {
    /// Creates a new instance with SIMD-optimized operations
    pub fn new(modulus: BigUint) -> Self {
        // Precompute constants
        let r_squared = Self::compute_r_squared(&modulus);
        let inv_modulus = Self::compute_inv_modulus(&modulus);
        
        Self {
            modulus,
            r_squared,
            inv_modulus,
        }
    }
    
    /// Computes R² mod N where R = 2^256
    fn compute_r_squared(modulus: &BigUint) -> BigUint {
        let r = BigUint::from(1u32) << 256;
        (r * r) % modulus
    }
    
    /// Computes the inverse of modulus for Montgomery reduction
    fn compute_inv_modulus(modulus: &BigUint) -> BigUint {
        // Using extended Euclidean algorithm
        let mut a = modulus.clone();
        let mut b = BigUint::from(1u32) << 256;
        let mut x0 = BigUint::from(0u32);
        let mut x1 = BigUint::from(1u32);
        
        while !a.is_zero() {
            let q = &b / &a;
            let r = &b % &a;
            let x2 = x0 - &q * &x1;
            
            b = a;
            a = r;
            x0 = x1;
            x1 = x2;
        }
        
        x0
    }
    
    /// SIMD-optimized field addition
    pub unsafe fn add(&self, a: &Fp, b: &Fp) -> Fp {
        if is_x86_feature_detected!("avx2") {
            self.add_avx2(a, b)
        } else {
            a + b
        }
    }
    
    /// SIMD-optimized field multiplication
    pub unsafe fn mul(&self, a: &Fp, b: &Fp) -> Fp {
        if is_x86_feature_detected!("avx2") {
            self.mul_avx2(a, b)
        } else {
            a * b
        }
    }
    
    /// AVX2-optimized field addition
    unsafe fn add_avx2(&self, a: &Fp, b: &Fp) -> Fp {
        let a_limbs = a.to_limbs();
        let b_limbs = b.to_limbs();
        
        // Load limbs into AVX2 registers
        let a_vec = _mm256_loadu_si256(a_limbs.as_ptr() as *const __m256i);
        let b_vec = _mm256_loadu_si256(b_limbs.as_ptr() as *const __m256i);
        
        // Perform addition
        let sum_vec = _mm256_add_epi64(a_vec, b_vec);
        
        // Store result
        let mut result_limbs = [0u64; 4];
        _mm256_storeu_si256(result_limbs.as_mut_ptr() as *mut __m256i, sum_vec);
        
        // Reduce if necessary
        let mut result = Fp::from_limbs(&result_limbs, &self.modulus);
        if result >= Fp::new(self.modulus.clone(), self.modulus.clone()) {
            result = result - Fp::new(self.modulus.clone(), self.modulus.clone());
        }
        
        result
    }
    
    /// AVX2-optimized field multiplication
    unsafe fn mul_avx2(&self, a: &Fp, b: &Fp) -> Fp {
        let a_limbs = a.to_limbs();
        let b_limbs = b.to_limbs();
        
        // Load limbs into AVX2 registers
        let a_vec = _mm256_loadu_si256(a_limbs.as_ptr() as *const __m256i);
        let b_vec = _mm256_loadu_si256(b_limbs.as_ptr() as *const __m256i);
        
        // Perform multiplication
        let mut result_limbs = [0u64; 8];
        for i in 0..4 {
            let a_limb = _mm256_set1_epi64x(a_limbs[i] as i64);
            let b_limb = _mm256_set1_epi64x(b_limbs[i] as i64);
            let product = _mm256_mul_epu32(a_limb, b_limb);
            
            // Store partial products
            let mut partial = [0u64; 4];
            _mm256_storeu_si256(partial.as_mut_ptr() as *mut __m256i, product);
            for j in 0..4 {
                result_limbs[i + j] += partial[j];
            }
        }
        
        // Reduce the result
        let mut result = Fp::from_limbs(&result_limbs[..4], &self.modulus);
        result = result % Fp::new(self.modulus.clone(), self.modulus.clone());
        
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use std::str::FromStr;

    #[test]
    fn test_simd_operations() {
        if !has_avx2() {
            println!("Skipping SIMD tests - AVX2 not available");
            return;
        }

        unsafe {
            let mut rng = rand::thread_rng();
            let a: [u64; 4] = std::array::from_fn(|_| rng.gen());
            let b: [u64; 4] = std::array::from_fn(|_| rng.gen());
            let modulus: [u64; 4] = [
                0xFFFFFFFFFFFFFFFF,
                0xFFFFFFFFFFFFFFFF,
                0xFFFFFFFFFFFFFFFF,
                0x0FFFFFFFFFFFFFFF,
            ];

            let result = field_mul_avx2(&a, &b, &modulus);
            assert!(result.iter().all(|&x| x < modulus[3]));

            let result = field_add_avx2(&a, &b, &modulus);
            assert!(result.iter().all(|&x| x < modulus[3]));
        }
    }

    #[test]
    fn test_simd_field_ops() {
        let modulus = BigUint::from_str(
            "21888242871839275222246405745257275088696311157297823662689037894645226208583"
        ).unwrap();
        
        let simd_ops = SimdFieldOps::new(modulus.clone());
        
        let a = Fp::new(BigUint::from(5u32), modulus.clone());
        let b = Fp::new(BigUint::from(3u32), modulus.clone());
        
        unsafe {
            // Test addition
            let sum = simd_ops.add(&a, &b);
            assert_eq!(sum, a + b);
            
            // Test multiplication
            let product = simd_ops.mul(&a, &b);
            assert_eq!(product, a * b);
        }
    }
} 