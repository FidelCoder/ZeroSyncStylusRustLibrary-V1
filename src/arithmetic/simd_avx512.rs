#[cfg(all(feature = "simd", target_arch = "x86_64"))]
use std::arch::x86_64::*;

/// Checks if AVX-512F is available
#[inline]
pub fn has_avx512f() -> bool {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        std::is_x86_feature_detected!("avx512f")
    }
    #[cfg(not(target_arch = "x86_64"))]
    false
}

/// Performs field multiplication using AVX-512
#[cfg(target_feature = "avx512f")]
#[target_feature(enable = "avx512f")]
pub unsafe fn field_mul_avx512(a: &[u64], b: &[u64], modulus: &[u64]) -> Vec<u64> {
    assert!(a.len() == 8 && b.len() == 8 && modulus.len() == 8);
    
    let mut result = vec![0u64; 8];
    
    // Load values into AVX-512 registers
    let a_vec = _mm512_loadu_epi64(a.as_ptr() as *const i64);
    let b_vec = _mm512_loadu_epi64(b.as_ptr() as *const i64);
    let m_vec = _mm512_loadu_epi64(modulus.as_ptr() as *const i64);
    
    // Multiply using AVX-512
    let mut t = _mm512_mul_epu32(a_vec, b_vec);
    
    // Reduce modulo m
    let mask = _mm512_cmpgt_epu64_mask(t, m_vec);
    if mask != 0 {
        t = _mm512_sub_epi64(t, m_vec);
    }
    
    // Store result
    _mm512_storeu_epi64(result.as_mut_ptr() as *mut i64, t);
    
    result
}

/// Performs field addition using AVX-512
#[cfg(target_feature = "avx512f")]
#[target_feature(enable = "avx512f")]
pub unsafe fn field_add_avx512(a: &[u64], b: &[u64], modulus: &[u64]) -> Vec<u64> {
    assert!(a.len() == 8 && b.len() == 8 && modulus.len() == 8);
    
    let mut result = vec![0u64; 8];
    
    // Load values into AVX-512 registers
    let a_vec = _mm512_loadu_epi64(a.as_ptr() as *const i64);
    let b_vec = _mm512_loadu_epi64(b.as_ptr() as *const i64);
    let m_vec = _mm512_loadu_epi64(modulus.as_ptr() as *const i64);
    
    // Add using AVX-512
    let mut sum = _mm512_add_epi64(a_vec, b_vec);
    
    // Reduce modulo m
    let mask = _mm512_cmpgt_epu64_mask(sum, m_vec);
    if mask != 0 {
        sum = _mm512_sub_epi64(sum, m_vec);
    }
    
    // Store result
    _mm512_storeu_epi64(result.as_mut_ptr() as *mut i64, sum);
    
    result
}

/// Performs vectorized Montgomery multiplication using AVX-512
#[cfg(target_feature = "avx512f")]
#[target_feature(enable = "avx512f")]
pub unsafe fn mont_mul_avx512(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64> {
    assert!(a.len() == 8 && b.len() == 8 && n.len() == 8 && n_prime.len() == 8);
    
    let mut result = vec![0u64; 8];
    
    // Load values into AVX-512 registers
    let a_vec = _mm512_loadu_epi64(a.as_ptr() as *const i64);
    let b_vec = _mm512_loadu_epi64(b.as_ptr() as *const i64);
    let n_vec = _mm512_loadu_epi64(n.as_ptr() as *const i64);
    let np_vec = _mm512_loadu_epi64(n_prime.as_ptr() as *const i64);
    
    // Compute t = a * b
    let t = _mm512_mul_epu32(a_vec, b_vec);
    
    // Compute m = t * n' mod R
    let m = _mm512_mul_epu32(t, np_vec);
    
    // Compute t = (t + m * n) / R
    let mn = _mm512_mul_epu32(m, n_vec);
    let sum = _mm512_add_epi64(t, mn);
    let r = _mm512_srli_epi64(sum, 32);
    
    // Final reduction
    let mask = _mm512_cmpgt_epu64_mask(r, n_vec);
    let final_result = if mask != 0 {
        _mm512_sub_epi64(r, n_vec)
    } else {
        r
    };
    
    // Store result
    _mm512_storeu_epi64(result.as_mut_ptr() as *mut i64, final_result);
    
    result
} 