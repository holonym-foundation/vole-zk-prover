use std::{str::FromStr};
use rand::prelude::*;
use ff::PrimeField;
use lazy_static::lazy_static;
use num_bigint::BigUint;

use crate::Fr;

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

/// Used for rejection sampling
fn u64s_overflow_field(x: &[u64; 4]) -> bool {
    (x[0] > MODULUS_AS_U64s[0])
 || (x[0] == MODULUS_AS_U64s[0] && x[1] > MODULUS_AS_U64s[1])
 || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] > MODULUS_AS_U64s[2]) 
 || (x[0] == MODULUS_AS_U64s[0] && x[1] == MODULUS_AS_U64s[1] && x[2] == MODULUS_AS_U64s[2] && x[3] >= MODULUS_AS_U64s[3])
}

// /// Used for rejection sampling
// fn u128s_overflow_field(x: &[u128; 2]) -> bool {
//     x[0] > MODULUS_AS_U128s[0] || (x[0] == MODULUS_AS_U128s[0] && x[1] >= MODULUS_AS_U128s[1])
// }

/// Highly efficient way of making an Fr from a big-endian u64 array representing that number
/// Does not check for overflow
/// TODO: Check whether ff or ff_ce has this method already and if not, contributing a more generic version.
pub fn fr_from_be_u64slice_unchecked(from: &[u64; 4]) -> Fr {
    Fr::from(from[3]) +
    Fr::from(from[2]) * *TWO_64 +
    Fr::from(from[1]) * *TWO_128 +
    Fr::from(from[0]) * *TWO_192
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


/// Newer method much faster: use a CSPRNG
/// Returns N Frs
/// As long as the adversary doesn't learn the seed (for a couple reasons throughout the protocol, they shouldn't), they can't predict any of the outputs
pub fn expand_seed_to_Fr_vec_faster(seed: &[u8; 32], num_outputs: usize) -> Vec<Fr> {
    let mut seed: <StdRng as SeedableRng>::Seed = *seed;
    let mut r = StdRng::from_seed(seed);
    let mut out: Vec<Fr> = Vec::with_capacity(num_outputs);
    // let mut output = r.gen::<Vec<[u64; 4]>>()
    // let mut output = [0u64; N];
    // let mut output = Vec::<u64>::with_capacity(N);

    for i in 0..num_outputs {
        // Rejection sample a random 254 bit big-endian value
        // Generate 4 u64s and try to put them into a Fr
        // AND the two most significant bits with 00, producing a random 254 bit big-endian value
        
        let mut candidate = [
            r.next_u64() & 0x3F,
            r.next_u64(),
            r.next_u64(),
            r.next_u64(),
        ];

        while u64s_overflow_field(&candidate) {
            candidate = [
                r.next_u64() & 0x3F,
                r.next_u64(),
                r.next_u64(),
                r.next_u64(),
            ];
        }
        out.push(fr_from_be_u64slice_unchecked(&candidate));
    }
    out
}

/// New method over 2x faster: generate the 0th output o_0 by doing sha256(seed), then o_1 by sha256(seed || o_0), then o_2 by sha256(seed || o_1)
/// Note: be careful of length extension attacks if modifying this
/// As long as the adversary doesn't learn the seed (for a couple reasons throughout the protocol, they shouldn't), they can't predict any of the outputs
/// TODO: see if native browser sha256 improves performance even compared to blake3
pub fn expand_seed_to_Fr_vec(seed: &[u8; 32], num_outputs: usize) -> Vec<Fr> {
    let mut outputs_bytes: Vec<[u8; 32]> = Vec::with_capacity(num_outputs);
    outputs_bytes.push(
        *blake3::hash(seed).as_bytes()
    );
    for i in 0..num_outputs-1 {
        outputs_bytes.push(
            *blake3::hash(
                &[
                    seed,
                    outputs_bytes[i].as_ref()
                ].concat()
            ).as_bytes()
        )
    }

    outputs_bytes.iter().map(|x|{
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
            b = *blake3::hash(b.as_ref()).as_bytes();
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
        fr_from_be_u64slice_unchecked(&as_4u64s)
    }).collect()
}

// /// Old method: make a tree of sha512 hashes and use two 256-bit halves each level as inputs to the next level:
// /// Takes a seed and uses it to create a Vec of Fr elements. Uses rejection sampling so isn't constant time!
// /// Returns 2^log_length Fr elements, represented as big-endian bytes
// /// Don't be tempted to return 2x-1 as many elements by including each recent_lvl in the output. That would make some resulting elements predictable from others
// /// TODO: speed up by not converting to/from string
// /// TODO: possible speed up in browser by using crypto sha256 function in the browser rather than wasm
// pub fn expand_seed_to_Fr_vec(seed: Vec<u8>, length: usize/*log_length: u32*/) -> Vec<Fr> {


//     let mut recent_lvl = vec![seed];
//     for i in 0..log_length {
//         let mut tmp = Vec::with_capacity(2usize.pow(i));
//         for j in 0..recent_lvl.len() {
//             // tmp.push(*blake3::hash(recent_lvl[j].as_ref()).as_bytes());
//             let mut hasher = Sha512::new();
//             hasher.update(&recent_lvl[j]);
//             let digest = hasher.finalize();
//             // let d2 = digest.as_slice();
//             tmp.push(digest[0..32].to_vec());
//             tmp.push(digest[32..64].to_vec());
//         }
//         recent_lvl = tmp;
//     }

//     // Collect the children as a Vec<Fr> using rejection sampling
//     recent_lvl.iter().map(|x| {
//         let mut b = x.clone();
//         // Prime modulus is only 254 bits; & the largest (0th) byte with 0x3F to ensure as close as possible to the prime modulus
//         b[0] = b[0] & 0x3F;
//         // let mut as_2u128s = [
//         //     u128::from_be_bytes(x[0..16].try_into().unwrap()),
//         //     u128::from_be_bytes(x[16..32].try_into().unwrap())
//         // ];
//         // b
//         // while u128s_overflow_field(&as_2u128s) {
//         let mut as_4u64s = [
//             u64::from_be_bytes(b[0..8].try_into().unwrap()),
//             u64::from_be_bytes(b[8..16].try_into().unwrap()),
//             u64::from_be_bytes(b[16..24].try_into().unwrap()),
//             u64::from_be_bytes(b[24..32].try_into().unwrap())
//         ];
//         while u64s_overflow_field(&as_4u64s) {
//             let mut hasher = Sha512::new();
//             hasher.update(b);
//             b = hasher.finalize().to_vec();
//             // Make it close to the 254-bit modulus
//             b[0] = b[0] & 0x3F;
//             // as_2u128s = [
//             //     u128::from_be_bytes(b[0..16].try_into().unwrap()),
//             //     u128::from_be_bytes(b[16..32].try_into().unwrap())
//             // ];
            
//             as_4u64s = [
//                 u64::from_be_bytes(b[0..8].try_into().unwrap()),
//                 u64::from_be_bytes(b[8..16].try_into().unwrap()),
//                 u64::from_be_bytes(b[16..24].try_into().unwrap()),
//                 u64::from_be_bytes(b[24..32].try_into().unwrap())
//             ];
//         }
//         fr_from_256bit(&as_4u64s)
//     }).collect()
    
// }

/// Instead of long vectors in most VOLE protocols, we're just doing a "vector" commitment to three values,
/// This means k for our SoftSpokenVOLE instantiation is 2, i.e. ∆ has just two bits of entropy.
/// Since we have to open and transmit all but one of the seeds, using a larger k for SoftSpokenVOLE doesn't save significant communication and solely wastes computation. 
pub fn commit_seeds(seed0: &Vec<u8>, seed1: &Vec<u8>) -> Vec<u8> {
    blake3::hash(
        &[
        *blake3::hash(seed0).as_bytes(),
        *blake3::hash(seed1).as_bytes(),
        ].concat()
    ).as_bytes().to_vec()
}

/// Just open one seed and hide the other since only two were committed :P. The proof an element is just the hash of the other hidden element
pub fn proof_for_revealed_seed(other_seed: &Vec<u8>) -> [u8; 32] {
    *blake3::hash(other_seed).as_bytes()
}

/// Verifies a proof for a committed seed
pub fn verify_proof_of_revealed_seed(
    commitment: &Vec<u8>,
    revealed_seed: &Vec<u8>, 
    revealed_seed_idx: bool, 
    proof: &[u8; 32]
) -> bool {
    let digest_of_revealed = *blake3::hash(revealed_seed).as_bytes();
    let preimage = 
    if revealed_seed_idx {
        [proof.clone(), digest_of_revealed].concat()
    } else {
        [digest_of_revealed, proof.clone()].concat()
    };
    blake3::hash(&preimage).as_bytes().to_vec() == *commitment
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_seed_expansion_len() {
        let seed = [0u8; 32];
        assert_eq!(super::expand_seed_to_Fr_vec(&seed, 1).len(), 1);
        assert_eq!(super::expand_seed_to_Fr_vec(&seed, 2).len(), 2);
        assert_eq!(super::expand_seed_to_Fr_vec(&seed, 4).len(), 4);
        assert_eq!(super::expand_seed_to_Fr_vec(&seed, 1000).len(), 1000);
    }

    // #[test]
    // fn test_seed_expansion_len() {
    //     let seed = [0u8; 32];
    //     assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 0).len(), 1);
    //     assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 1).len(), 2);
    //     assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 2).len(), 4);
    //     assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 10).len(), 1024);
    // }

    #[test]
    // TODO: add more edge cases
    fn test_fr_from_512bits() {
        // Try with modulus + 1
        let mut x = [0u64; 4];
        x[0] = 0x30644E72E131A029;
        x[1] = 0xB85045B68181585D;
        x[2] = 0x2833E84879B97091;
        x[3] = 0x43E1F593F0000002;

        assert_eq!(fr_from_be_u64slice_unchecked(&x), Fr::from_str_vartime("1").unwrap());

        x[0] = 0x30644E72E131A029;
        x[1] = 0xB85045B68181585D;
        x[2] = 0x2833E84879B97091;
        x[3] = 0x43E1F593F0000001;

        assert_eq!(fr_from_be_u64slice_unchecked(&x), Fr::from_str_vartime("0").unwrap());

        x[0] = 0x0;
        x[1] = 0x0;
        x[2] = 0x0;
        x[3] = 0x03;

        assert_eq!(fr_from_be_u64slice_unchecked(&x), Fr::from_str_vartime("3").unwrap());
    }
    
    #[test]
    fn test_seed_commit_prove() {
        let seed0 = [5u8; 32].to_vec();
        let seed1 = [6u8; 32].to_vec();
        let commitment = commit_seeds(&seed0, &seed1);

        let proof0 = proof_for_revealed_seed(&seed1);
        let proof1 = proof_for_revealed_seed(&seed0);

        assert!(verify_proof_of_revealed_seed(&commitment, &seed0, false, &proof0));
        assert!(!verify_proof_of_revealed_seed(&commitment, &seed0, true, &proof0));

        assert!(verify_proof_of_revealed_seed(&commitment, &seed1, true, &proof1));
        assert!(!verify_proof_of_revealed_seed(&commitment, &seed1, false, &proof1));

        assert!(!verify_proof_of_revealed_seed(&commitment, &seed0, true, &proof1));
        assert!(!verify_proof_of_revealed_seed(&commitment, &seed0, false, &proof1));

    }
}