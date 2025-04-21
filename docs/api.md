# ZeroSync API Documentation

## Field Arithmetic

### Prime Field Operations (`arithmetic::field::Fp`)

The `Fp` struct represents elements in a prime field using Montgomery arithmetic for efficient operations.

```rust
pub struct Fp {
    limbs: Vec<u64>,      // Value in Montgomery form
    constants: MontgomeryConstants,
}
```

#### Constructor

```rust
/// Creates a new field element
pub fn new(value: BigUint, modulus: BigUint) -> Self
```

Example:
```rust
let modulus = BigUint::from(17u32);
let element = Fp::new(BigUint::from(5u32), modulus);
```

#### Field Operations

All operations are constant-time and automatically handle Montgomery form conversion:

```rust
// Addition
let sum = a + b;  // Modular addition

// Subtraction
let diff = a - b;  // Modular subtraction

// Multiplication
let product = a * b;  // Montgomery multiplication

// Exponentiation
let power = a.pow(exponent);  // Efficient square-and-multiply
```

### SIMD Optimizations (`arithmetic::simd`)

AVX2-accelerated field operations for 4-way parallel computation:

```rust
#[target_feature(enable = "avx2")]
pub unsafe fn field_mul_avx2(a: &[u64; 4], b: &[u64; 4], modulus: &[u64; 4]) -> [u64; 4]

#[target_feature(enable = "avx2")]
pub unsafe fn field_add_avx2(a: &[u64; 4], b: &[u64; 4], modulus: &[u64; 4]) -> [u64; 4]
```

Feature detection:
```rust
if has_avx2() {
    // Use SIMD operations
} else {
    // Use standard operations
}
```

### Montgomery Arithmetic (`arithmetic::montgomery`)

Low-level Montgomery arithmetic implementation:

```rust
/// Optimized Montgomery multiplication
pub fn mont_mul(a: &[u64], b: &[u64], n: &[u64], n_prime: &[u64]) -> Vec<u64>

/// Constant-time comparison
pub fn ct_lt(a: &[u64], b: &[u64]) -> bool
```

## Performance Considerations

### Memory Management

- Use stack allocation for small field elements
- Avoid unnecessary cloning of field elements
- Reuse temporary variables in loops

Example:
```rust
// Efficient batch operation
let mut acc = Fp::new(BigUint::from(1u32), modulus.clone());
for element in elements {
    acc = acc * element;  // Reuses temporary variables
}
```

### SIMD Optimization

- Align data to 32-byte boundaries for AVX2
- Process multiple elements in parallel
- Use runtime feature detection

Example:
```rust
if has_avx2() {
    unsafe {
        let result = field_mul_avx2(&a, &b, &modulus);
    }
} else {
    let result = standard_multiplication(&a, &b, &modulus);
}
```

## Error Handling

Operations that may fail return `Option` or `Result`:

```rust
// Returns None if element has no inverse
pub fn inverse(&self) -> Option<Self>

// Returns Err for invalid inputs
pub fn from_bytes(bytes: &[u8]) -> Result<Self, Error>
```

## Thread Safety

- All field elements are `Send` and `Sync`
- Constants are cached thread-locally
- No global mutable state

## Gas Optimization

Field operations are optimized for Arbitrum Stylus:

| Operation    | Gas Cost | Notes                    |
|-------------|----------|--------------------------|
| Addition    | 500      | Constant gas cost        |
| Multiply    | 1500     | Uses Montgomery form     |
| Inverse     | 5000     | Uses binary GCD         |
| Batch Ops   | Varies   | Amortized cost per op   |

## Security Notes

1. All operations are constant-time
2. No secret-dependent branches
3. Memory is cleared after use
4. Regular security audits

## Examples

### Basic Field Arithmetic

```rust
use zerosync::arithmetic::field::Fp;
use num_bigint::BigUint;

// Create field elements
let modulus = BigUint::from(17u32);
let a = Fp::new(BigUint::from(5u32), modulus.clone());
let b = Fp::new(BigUint::from(3u32), modulus);

// Field operations
let sum = a.clone() + b.clone();
let product = a * b;
let power = sum.pow(5);
```

### SIMD Operations

```rust
use zerosync::arithmetic::simd::{field_mul_avx2, has_avx2};

if has_avx2() {
    unsafe {
        let result = field_mul_avx2(&values1, &values2, &modulus);
        // Process four field elements in parallel
    }
}
```

### Batch Operations

```rust
// Efficient batch multiplication
let product = elements.iter().fold(
    Fp::new(BigUint::from(1u32), modulus.clone()),
    |acc, x| acc * x.clone()
);
``` 