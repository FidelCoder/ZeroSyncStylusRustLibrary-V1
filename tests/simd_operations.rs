#[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
use zerosync::arithmetic::simd::{field_mul_avx2, field_add_avx2, has_avx2};

#[cfg(all(feature = "simd", any(target_arch = "x86", target_arch = "x86_64")))]
#[test]
fn test_avx2_detection() {
    // Just verify that the function runs without panicking
    let _ = has_avx2();
}

#[test]
#[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
fn test_simd_operations() {
    if !has_avx2() {
        println!("Skipping SIMD tests - AVX2 not available");
        return;
    }

    // Test vectors
    let a = vec![1u64, 2, 3, 4];
    let b = vec![5u64, 6, 7, 8];
    let modulus = vec![17u64, 17, 17, 17];

    unsafe {
        // Test field multiplication
        let result = field_mul_avx2(&a, &b, &modulus);
        assert_eq!(result.len(), 4);
        for i in 0..4 {
            assert!(result[i] < modulus[i]);
        }

        // Test field addition
        let result = field_add_avx2(&a, &b, &modulus);
        assert_eq!(result.len(), 4);
        for i in 0..4 {
            assert!(result[i] < modulus[i]);
        }
    }
}

#[test]
#[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
fn test_simd_field_addition() {
    if !has_avx2() {
        println!("Skipping SIMD tests - AVX2 not available");
        return;
    }

    let a = vec![1u64, 2, 3, 4];
    let b = vec![5u64, 6, 7, 8];
    let modulus = vec![17u64, 17, 17, 17];

    unsafe {
        let result = field_add_avx2(&a, &b, &modulus);
        
        // Test basic addition
        assert_eq!(result[0], (1 + 5) % 17);  // 6
        assert_eq!(result[1], (2 + 6) % 17);  // 8
        assert_eq!(result[2], (3 + 7) % 17);  // 10
        assert_eq!(result[3], (4 + 8) % 17);  // 12
    }
}

#[test]
#[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
fn test_simd_field_multiplication() {
    if !has_avx2() {
        println!("Skipping SIMD tests - AVX2 not available");
        return;
    }

    let a = vec![1u64, 2, 3, 4];
    let b = vec![2u64, 3, 4, 5];
    let modulus = vec![17u64, 17, 17, 17];

    unsafe {
        let result = field_mul_avx2(&a, &b, &modulus);
        
        // Test basic multiplication
        assert_eq!(result[0], (1 * 2) % 17);  // 2
        assert_eq!(result[1], (2 * 3) % 17);  // 6
        assert_eq!(result[2], (3 * 4) % 17);  // 12
        assert_eq!(result[3], (4 * 5) % 17);  // 3
    }
}

#[test]
#[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
fn test_simd_edge_cases() {
    if !has_avx2() {
        println!("Skipping SIMD tests - AVX2 not available");
        return;
    }

    let modulus = vec![17u64, 17, 17, 17];

    unsafe {
        // Test multiplication by zero
        let a = vec![1u64, 2, 3, 4];
        let zero = vec![0u64, 0, 0, 0];
        let result = field_mul_avx2(&a, &zero, &modulus);
        assert_eq!(result, vec![0u64, 0, 0, 0]);

        // Test multiplication by one
        let one = vec![1u64, 1, 1, 1];
        let result = field_mul_avx2(&a, &one, &modulus);
        for i in 0..4 {
            assert_eq!(result[i], a[i] % 17);
        }

        // Test addition with zero
        let result = field_add_avx2(&a, &zero, &modulus);
        for i in 0..4 {
            assert_eq!(result[i], a[i] % 17);
        }
    }
}

#[test]
#[cfg(target_feature = "avx2")]
fn test_simd_associativity() {
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];
    let a = [1u64, 2, 3, 4];
    let b = [2u64, 3, 4, 5];
    let c = [3u64, 4, 5, 6];

    unsafe {
        // Test (a + b) + c = a + (b + c)
        let ab = field_add_avx2(&a, &b, &modulus);
        let abc1 = field_add_avx2(&ab, &c, &modulus);

        let bc = field_add_avx2(&b, &c, &modulus);
        let abc2 = field_add_avx2(&a, &bc, &modulus);

        assert_eq!(abc1, abc2);

        // Test (a * b) * c = a * (b * c)
        let ab = field_mul_avx2(&a, &b, &modulus);
        let abc1 = field_mul_avx2(&ab, &c, &modulus);

        let bc = field_mul_avx2(&b, &c, &modulus);
        let abc2 = field_mul_avx2(&a, &bc, &modulus);

        assert_eq!(abc1, abc2);
    }
}

#[test]
#[cfg(target_feature = "avx2")]
fn test_simd_distributivity() {
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];
    let a = [1u64, 2, 3, 4];
    let b = [2u64, 3, 4, 5];
    let c = [3u64, 4, 5, 6];

    unsafe {
        // Test a * (b + c) = (a * b) + (a * c)
        let b_plus_c = field_add_avx2(&b, &c, &modulus);
        let left = field_mul_avx2(&a, &b_plus_c, &modulus);

        let a_times_b = field_mul_avx2(&a, &b, &modulus);
        let a_times_c = field_mul_avx2(&a, &c, &modulus);
        let right = field_add_avx2(&a_times_b, &a_times_c, &modulus);

        assert_eq!(left, right);
    }
} 