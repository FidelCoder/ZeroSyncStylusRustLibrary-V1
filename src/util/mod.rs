/// Utility functions for the ZeroSync library

/// Converts a byte slice to a fixed-size array
pub fn to_fixed_bytes<const N: usize>(bytes: &[u8]) -> [u8; N] {
    let mut result = [0u8; N];
    result[..bytes.len()].copy_from_slice(bytes);
    result
}

/// Converts a slice of bytes to a vector of u64 limbs
pub fn bytes_to_limbs(bytes: &[u8]) -> Vec<u64> {
    bytes
        .chunks(8)
        .map(|chunk| {
            let mut limb = 0u64;
            for (i, &byte) in chunk.iter().enumerate() {
                limb |= (byte as u64) << (i * 8);
            }
            limb
        })
        .collect()
} 