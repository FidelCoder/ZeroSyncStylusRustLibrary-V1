//! ZeroSync: Zero-knowledge proof library optimized for Arbitrum Stylus
//! 
//! This library provides efficient implementations of ZK proving systems
//! with a focus on gas optimization and developer experience.

pub mod arithmetic;
pub mod curve;
pub mod util;
pub mod polynomial;

// Re-export commonly used types
pub use arithmetic::traits::{Field, PrimeField};
pub use curve::bn254::{BN254, G1Affine, G2Affine};

/// Feature flags
#[cfg(feature = "parallel")]
pub use rayon;

// Version information
pub const VERSION: &str = env!("CARGO_PKG_VERSION"); 