use std::{str::FromStr, ops::MulAssign};

use ff::Field;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use sha2::{Digest,Sha512};

#[macro_use]
extern crate ff;

use crate::ff::PrimeField;

/// TODO: Check whether ff or ff_ce has this method already and if not, contribute it.
fn fr_from_256bit(from: &[u64; 4]) -> Fr {
    Fr::from(from[3]) +
    Fr::from(from[2]) * *TWO_64 +
    Fr::from(from[1]) * *TWO_128 +
    Fr::from(from[0]) * *TWO_192
}


#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
// Important this matches the endianness of MODULUS_AS_U128s
#[PrimeFieldReprEndianness = "big"]
pub struct Fr([u64; 4]);

lazy_static! {
    pub static ref FIELD_MODULUS: BigUint = BigUint::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495617").unwrap();
    /// 2^64
    pub static ref TWO_64: Fr = Fr::from_str_vartime("18446744073709551616").unwrap();
    /// 2^128
    pub static ref TWO_128: Fr = Fr::from_str_vartime("340282366920938463463374607431768211456").unwrap();
    /// 2^192
    pub static ref TWO_192: Fr = Fr::from_str_vartime("6277101735386680763835789423207666416102355444464034512896").unwrap();
}
/// U128 representation of 21888242871839275222246405745257275088548364400416034343698204186575808495617
const MODULUS_AS_U64s: [u64; 4] = [
    0x30644E72E131A029,
    0xB85045B68181585D, 
    0x2833E84879B97091,
    0x43E1F593F0000001
];

// /// Used for rejection sampling
// fn u128s_overflow_field(x: &[u128; 2]) -> bool {
//     x[0] > MODULUS_AS_U128s[0] || (x[0] == MODULUS_AS_U128s[0] && x[1] >= MODULUS_AS_U128s[1])
// }

fn u64s_overflow_field(x: &[u64; 4]) -> bool {
       (x[0] > MODULUS_AS_U64s[0])
    || (x[0] == MODULUS_AS_U64s[0] && x[1] > MODULUS_AS_U64s[1])
    || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] > MODULUS_AS_U64s[2]) 
    || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] == MODULUS_AS_U64s[2] && x[3] >= MODULUS_AS_U64s[3])
}
// /// Generates 2^log_length Fr elements. 2^long_length means it can be put into a binary merkle tree
// pub fn rand_fr_vec(log_len: u32) -> Vec<Fr> {
//     let len = 2usize.pow(log_len);
//     let mut rng = rand::thread_rng();
//     let mut ret = Vec::with_capacity(len);
//     for _i in 0..len {
//         ret.push(Fr::random(&mut rng))
//     }
//     ret
// }


/// Takes a seed and uses it to create a Vec of Fr elements. Uses rejection sampling so isn't constant time!
/// Returns 2^log_length Fr elements, represented as big-endian bytes
/// Don't be tempted to return 2x-1 as many elements by including each recent_lvl in the output. That would make some resulting elements predictable from others
/// TODO: speed up by not converting to/from string
/// TODO: possible speed up in browser by using crypto sha256 function in the browser rather than wasm
pub fn expand_seed_to_Fr_vec(seed: Vec<u8>, log_length: u32) -> Vec<Fr> {
    let mut recent_lvl = vec![seed];
    for i in 0..log_length {
        let mut tmp = Vec::with_capacity(2usize.pow(i));
        for j in 0..recent_lvl.len() {
            // tmp.push(*blake3::hash(recent_lvl[j].as_ref()).as_bytes());
            let mut hasher = Sha512::new();
            hasher.update(&recent_lvl[j]);
            let digest = hasher.finalize();
            // let d2 = digest.as_slice();
            tmp.push(digest[0..32].to_vec());
            tmp.push(digest[32..64].to_vec());
        }
        recent_lvl = tmp;
    }

    // Collect the children as a Vec<Fr> using rejection sampling
    recent_lvl.iter().map(|x| {
        let mut b = x.clone();
        // Prime modulus is only 254 bits; & the largest (0th) byte with 0x3F to ensure as close as possible to the prime modulus
        b[0] = b[0] & 0x3F;
        // let mut as_2u128s = [
        //     u128::from_be_bytes(x[0..16].try_into().unwrap()),
        //     u128::from_be_bytes(x[16..32].try_into().unwrap())
        // ];
        // b
        // while u128s_overflow_field(&as_2u128s) {
        let mut as_4u64s = [
            u64::from_be_bytes(b[0..8].try_into().unwrap()),
            u64::from_be_bytes(b[8..16].try_into().unwrap()),
            u64::from_be_bytes(b[16..24].try_into().unwrap()),
            u64::from_be_bytes(b[24..32].try_into().unwrap())
        ];
        while u64s_overflow_field(&as_4u64s) {
            let mut hasher = Sha512::new();
            hasher.update(b);
            b = hasher.finalize().to_vec();
            // Make it close to the 254-bit modulus
            b[0] = b[0] & 0x3F;
            // as_2u128s = [
            //     u128::from_be_bytes(b[0..16].try_into().unwrap()),
            //     u128::from_be_bytes(b[16..32].try_into().unwrap())
            // ];
            
            as_4u64s = [
                u64::from_be_bytes(b[0..8].try_into().unwrap()),
                u64::from_be_bytes(b[8..16].try_into().unwrap()),
                u64::from_be_bytes(b[16..24].try_into().unwrap()),
                u64::from_be_bytes(b[24..32].try_into().unwrap())
            ];
        }
        fr_from_256bit(&as_4u64s)
    }).collect()
    
}

/// Instead of long vectors in most VOLE protocols, we're just doing a "vector" commitment to three values,
/// This means lambda for our SoftSpokenVOLE instantiation is 3, i.e. ∆ has just under two bits of entropy.
/// Since we have to open and transmit all but one of the seeds, a larger lambda for SoftSpokenVOLE doesn't save significant communication and wastes computation. But anything smaller than 3 is insecure because with 2, there is a 50% chance adversary learns V
pub fn commit_seeds(seed1: &Vec<u8>, seed2: &Vec<u8>) -> Vec<u8> {
    let mut hasher1 = Sha512::new();
    hasher1.update(seed1);
    let d1 = hasher1.finalize();

    let mut hasher2 = Sha512::new();
    hasher2.update(seed2);
    let d2 = hasher2.finalize();

    let mut hasher = Sha512::new();
    hasher.update(d1);
    hasher.update(d2);
    hasher.finalize().to_vec()
}
/// Just open one seed and hide the other since only two were committed :P. The proof an element is just the hash of the other hidden element
pub fn proof_for_revealed_seed(other_seed: &Vec<u8>) -> Vec<u8> {
    let mut hasher = Sha512::new();
    hasher.update(other_seed);
    hasher.finalize().to_vec()
}

/// Verifies a proof for a committed seed
pub fn verify_proof_of_revealed_seed(
    commitment: &Vec<u8>,
    revealed_seed: &Vec<u8>, 
    revealed_seed_idx: bool, 
    proof: &Vec<u8>
) -> bool {
    let mut hasher = Sha512::new();
    if revealed_seed_idx {
        hasher.update(proof);
        hasher.update(revealed_seed);
    } else {
        hasher.update(revealed_seed);
        hasher.update(proof);
    }
    &hasher.finalize().to_vec() == commitment
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_seed_expansion_len() {
        let seed = [0u8; 32];
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 0).len(), 1);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 1).len(), 2);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 2).len(), 4);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 10).len(), 1024);
    }

    #[test]
    // TODO: add more edge cases
    fn test_fr_from_512bits() {
        // Try with modulus + 1
        let mut x = [0u64; 4];
        x[0] = 0x30644E72E131A029;
        x[1] = 0xB85045B68181585D;
        x[2] = 0x2833E84879B97091;
        x[3] = 0x43E1F593F0000002;

        assert_eq!(fr_from_256bit(&x), Fr::from_str_vartime("1").unwrap());

        x[0] = 0x30644E72E131A029;
        x[1] = 0xB85045B68181585D;
        x[2] = 0x2833E84879B97091;
        x[3] = 0x43E1F593F0000001;

        assert_eq!(fr_from_256bit(&x), Fr::from_str_vartime("0").unwrap());

        x[0] = 0x0;
        x[1] = 0x0;
        x[2] = 0x0;
        x[3] = 0x03;

        assert_eq!(fr_from_256bit(&x), Fr::from_str_vartime("3").unwrap());
    }
}