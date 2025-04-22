#![cfg(feature = "simd")]
#![cfg(any(target_arch = "x86", target_arch = "x86_64"))]

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::arithmetic::field::Fp;

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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

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
} 