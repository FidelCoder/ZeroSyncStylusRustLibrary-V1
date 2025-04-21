use zerosync::arithmetic::simd::{field_mul_avx2, field_add_avx2, has_avx2};

#[test]
fn test_avx2_detection() {
    // Just verify that the function runs without panicking
    let _ = has_avx2();
}

#[test]
#[cfg(target_feature = "avx2")]
fn test_simd_field_addition() {
    let a = [1u64, 2, 3, 4];
    let b = [5u64, 6, 7, 8];
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];

    unsafe {
        let result = field_add_avx2(&a, &b, &modulus);
        
        // Test basic addition
        assert_eq!(result[0], 6);  // 1 + 5
        assert_eq!(result[1], 8);  // 2 + 6
        assert_eq!(result[2], 10); // 3 + 7
        assert_eq!(result[3], 12); // 4 + 8
        
        // Test modular reduction
        let large_a = [
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0x0FFFFFFFFFFFFFFF,
        ];
        let result = field_add_avx2(&large_a, &b, &modulus);
        // Verify results are properly reduced modulo p
        for i in 0..4 {
            assert!(result[i] < modulus[i]);
        }
    }
}

#[test]
#[cfg(target_feature = "avx2")]
fn test_simd_field_multiplication() {
    let a = [1u64, 2, 3, 4];
    let b = [2u64, 3, 4, 5];
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];

    unsafe {
        let result = field_mul_avx2(&a, &b, &modulus);
        
        // Test basic multiplication
        assert_eq!(result[0], 2);  // 1 * 2
        assert_eq!(result[1], 6);  // 2 * 3
        assert_eq!(result[2], 12); // 3 * 4
        assert_eq!(result[3], 20); // 4 * 5
        
        // Test modular reduction
        let large_a = [
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0x0FFFFFFFFFFFFFFF,
        ];
        let result = field_mul_avx2(&large_a, &b, &modulus);
        // Verify results are properly reduced modulo p
        for i in 0..4 {
            assert!(result[i] < modulus[i]);
        }
    }
}

#[test]
#[cfg(target_feature = "avx2")]
fn test_simd_edge_cases() {
    let modulus = [
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0xFFFFFFFFFFFFFFFF,
        0x0FFFFFFFFFFFFFFF,
    ];

    unsafe {
        // Test multiplication by zero
        let a = [1u64, 2, 3, 4];
        let zero = [0u64; 4];
        let result = field_mul_avx2(&a, &zero, &modulus);
        assert_eq!(result, [0u64; 4]);

        // Test multiplication by one
        let one = [1u64, 1, 1, 1];
        let result = field_mul_avx2(&a, &one, &modulus);
        assert_eq!(result, a);

        // Test addition with zero
        let result = field_add_avx2(&a, &zero, &modulus);
        assert_eq!(result, a);
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