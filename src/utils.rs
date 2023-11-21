use ff::PrimeField;

use crate::{Fr, FrRepr};

/// U128 representation of 21888242871839275222246405745257275088548364400416034343698204186575808495617
const MODULUS_AS_U64s: [u64; 4] = [
    0x30644E72E131A029,
    0xB85045B68181585D, 
    0x2833E84879B97091,
    0x43E1F593F0000001
];

/// Used for rejection sampling
pub fn u64s_overflow_field(x: &[u64; 4]) -> bool {
    (x[0] > MODULUS_AS_U64s[0])
 || (x[0] == MODULUS_AS_U64s[0] && x[1] > MODULUS_AS_U64s[1])
 || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] > MODULUS_AS_U64s[2]) 
 || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] == MODULUS_AS_U64s[2] && x[3] >= MODULUS_AS_U64s[3])
}

/// Efficient way of making an Fr from a big-endian u64 array representing that number
pub fn fr_from_be_u64_slice(from: &[u64; 4]) -> Fr {
    let mut repr = [0u8; 32];
    repr[0..8].copy_from_slice(&from[0].to_be_bytes());
    repr[8..16].copy_from_slice(&from[1].to_be_bytes());
    repr[16..24].copy_from_slice(&from[2].to_be_bytes());
    repr[24..32].copy_from_slice(&from[3].to_be_bytes());
    Fr::from_repr_vartime(FrRepr(repr)).unwrap()
}
/// Used for rejection sampling.
/// Converts a [u8; 32] to a [u64; 4] and removes the (big-endian) first two bits of the first u64 so it *might* fit inside the modulus
pub fn truncate_u8_32_to_254_bit_u64s_be(x: &[u8; 32]) -> [u64; 4] {
    let mut u64s = [
        u64::from_be_bytes(x[0..8].try_into().unwrap()), 
        u64::from_be_bytes(x[8..16].try_into().unwrap()), 
        u64::from_be_bytes(x[16..24].try_into().unwrap()), 
        u64::from_be_bytes(x[24..32].try_into().unwrap())
    ];
    u64s[0] &= 0x3FFFFFFFFFFFFFFF;
    u64s
}

/// Rejection samples to fit a 32 byte number into a field
/// Is geared for moderate performance but can be improved in some cases
pub fn rejection_sample_u8s(input: &[u8; 32]) -> Fr {
    let mut attempt = input.clone();
    let mut chopped_attempt = truncate_u8_32_to_254_bit_u64s_be(&attempt);
    while u64s_overflow_field(&chopped_attempt) {
        attempt = *blake3::hash(&attempt).as_bytes();
        chopped_attempt = truncate_u8_32_to_254_bit_u64s_be(&attempt);
    }
    fr_from_be_u64_slice(&chopped_attempt)
}