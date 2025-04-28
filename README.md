# ZeroSync: Zero-Knowledge Proofs for Arbitrum Stylus

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/zerosync/zerosync/workflows/ZeroSync%20CI/badge.svg)](https://github.com/zerosync/zerosync/actions)

A high-performance Rust library for zero-knowledge proofs optimized for Arbitrum Stylus, featuring efficient field arithmetic and curve operations.

## Features

- 🚀 Optimized field arithmetic with Montgomery form
- 🔒 Constant-time operations for cryptographic security
- ⚡ Lazy reduction for fast modular operations
- 📊 Comprehensive benchmarking suite
- 🛠️ BN254 curve implementation with G1 and G2 support

## Prerequisites

- Rust 1.70 or later
- Cargo (Rust's package manager)
- Basic understanding of elliptic curve cryptography

## Installation

1. Add ZeroSync to your `Cargo.toml`:

```toml
[dependencies]
zerosync = "0.1.0"
num-bigint = "0.4"  # Required for BigUint operations
```

2. Build the library:

```bash
cargo build --release
```

## Quick Start

### Field Arithmetic

```rust
use zerosync::arithmetic::field::Fp;
use num_bigint::BigUint;
use std::str::FromStr;

fn main() {
    // Create field elements
    let modulus = BigUint::from_str(
        "21888242871839275222246405745257275088696311157297823662689037894645226208583"
    ).unwrap();
    
    let a = Fp::new(BigUint::from(5u32), modulus.clone());
    let b = Fp::new(BigUint::from(3u32), modulus);

    // Perform field operations
    let sum = a.clone() + b.clone();
    let product = a * b;
    
    println!("Sum: {}", sum.from_montgomery());
    println!("Product: {}", product.from_montgomery());
}
```

### Curve Operations

```rust
use zerosync::curve::bn254::{BN254, G1Affine};

fn main() {
    // Create a BN254 curve instance
    let curve = BN254::new();
    
    // Get the generator point for G1
    let g = curve.g1_generator();
    
    // Perform point operations
    let p2 = g.clone() + g.clone(); // Point addition
    let p3 = g * 3;                // Scalar multiplication
    
    // Verify points are on the curve
    assert!(curve.is_on_curve(&p2));
    assert!(curve.is_on_curve(&p3));
}
```

## Examples

The library includes several example programs to help you get started:

1. **Basic Field Operations** (`examples/basic_field_ops.rs`):
   ```bash
   cargo run --example basic_field_ops
   ```

2. **Simplified Field Implementation** (`examples/simplified_field.rs`):
   ```bash
   cargo run --example simplified_field
   ```

3. **Simplified Curve Operations** (`examples/simplified_curve.rs`):
   ```bash
   cargo run --example simplified_curve
   ```

4. **Performance Benchmarking** (`examples/simple_benchmark.rs`):
   ```bash
   cargo run --release --example simple_benchmark
   ```

## Performance

The library includes optimized implementations of field arithmetic and curve operations. Here are some typical performance metrics:

| Operation | Performance |
|-----------|-------------|
| Field Addition | ~150ns |
| Field Multiplication | ~120ns |
| Field Squaring | ~130ns |
| Point Addition | ~500ns |
| Scalar Multiplication | ~100μs |

Note: Performance may vary based on your hardware and input sizes.

## Architecture

The library is structured in layers:

1. **Field Arithmetic Layer** (`src/arithmetic/`)
   - Prime field operations with Montgomery form
   - Lazy reduction techniques
   - Efficient modular arithmetic

2. **Curve Operations Layer** (`src/curve/`)
   - BN254 curve implementation
   - Point arithmetic with lazy reduction
   - Extension field (Fp2) operations

## Testing

Run the test suite:

```bash
cargo test
```

## Documentation

- [API Documentation](docs/api.md)
- [Architecture Guide](docs/architecture.md)
- [Performance Guide](docs/performance.md)
- [Security Guide](docs/security.md)
- [Technical Specification](docs/milestone1_spec.md)

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
