use crate::arithmetic::field::Fp;

/// BN254 elliptic curve implementation
#[allow(dead_code)]
pub struct BN254;

/// A point in G1 represented in affine coordinates
#[allow(dead_code)]
pub struct G1Affine {
    x: Fp,
    y: Fp,
    infinity: bool,
}

/// A point in G2 represented in affine coordinates
#[allow(dead_code)]
pub struct G2Affine {
    x: Fp,
    y: Fp,
    infinity: bool,
} 