use std::{str::FromStr};
use rand::prelude::*;
use ff::{PrimeField, Field};
use lazy_static::lazy_static;
use num_bigint::BigUint;
use rand_chacha::ChaCha12Rng;

use crate::{Fr, FrVec};

lazy_static! {
    pub static ref FIELD_MODULUS: BigUint = BigUint::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495617").unwrap();
    /// 2^64
    pub static ref TWO_64: Fr = Fr::from_str_vartime("18446744073709551616").unwrap();
    /// 2^128
    pub static ref TWO_128: Fr = Fr::from_str_vartime("340282366920938463463374607431768211456").unwrap();
    /// 2^192
    pub static ref TWO_192: Fr = Fr::from_str_vartime("6277101735386680763835789423207666416102355444464034512896").unwrap();
}


/// Newer method much faster: use a CSPRNG
/// Returns N Frs
/// As long as the adversary doesn't learn the seed (for a couple reasons throughout the protocol, they shouldn't), they can't predict any of the outputs
pub fn expand_seed_to_Fr_vec(seed: [u8; 32], num_outputs: usize) -> FrVec {
    let mut r = ChaCha12Rng::from_seed(seed);
    let mut out: Vec<Fr> = Vec::with_capacity(num_outputs);

    for _i in 0..num_outputs {
        out.push(Fr::random(&mut r));
    }

    FrVec(out)
}

/// Instead of long vectors in most VOLE protocols, we're just doing a "vector" commitment to two values,
/// This means k for our SoftSpokenVOLE instantiation is 2, i.e. ∆ has just two bits of entropy.
/// Since we have to open and transmit all but one of the seeds, using a larger k for SoftSpokenVOLE doesn't save significant communication and solely wastes computation. 
pub fn commit_seeds<T: AsRef<[u8]>>(seed0: &T, seed1: &T) -> [u8; 32] {
    *blake3::hash(
        &[
        *blake3::hash(seed0.as_ref()).as_bytes(),
        *blake3::hash(seed1.as_ref()).as_bytes(),
        ].concat()
    ).as_bytes()
}
/// Makes one hash of many seed commitments
pub fn commit_seed_commitments<T: AsRef<[u8]>>(comms: &Vec<T>) -> [u8; 32] {
    let mut hasher = blake3::Hasher::new();
    comms.iter().for_each(|c|{hasher.update(c.as_ref());});
    *hasher.finalize().as_bytes()
}

/// Just open one seed and hide the other since only two were committed :P. The proof an element is just the hash of the other hidden element
pub fn proof_for_revealed_seed(other_seed: &[u8; 32]) -> [u8; 32] {
    *blake3::hash(other_seed).as_bytes()
}

/// Verifies a proof for a committed seed
pub fn verify_proof_of_revealed_seed(
    commitment: &[u8; 32],
    revealed_seed: &[u8; 32], 
    revealed_seed_idx: bool, 
    proof: &[u8; 32]
) -> bool {
    &reconstruct_commitment(revealed_seed, revealed_seed_idx, proof) == commitment
}
/// Reconstructs a commitment to a seed given a known seed and a proof for the other seed. If this commitment checks out the proof is valid
pub fn reconstruct_commitment(
    revealed_seed: &[u8; 32], 
    revealed_seed_idx: bool, 
    proof: &[u8; 32]
) -> [u8; 32] {
    let digest_of_revealed = *blake3::hash(revealed_seed).as_bytes();
    let preimage = 
    if revealed_seed_idx {
        [proof.clone(), digest_of_revealed].concat()
    } else {
        [digest_of_revealed, proof.clone()].concat()
    };
    *blake3::hash(&preimage).as_bytes()
}

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    fn test_seed_expansion_len() {
        let seed = [0u8; 32];
        assert_eq!(super::expand_seed_to_Fr_vec(seed.clone(), 1).0.len(), 1);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.clone(), 2).0.len(), 2);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.clone(), 4).0.len(), 4);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.clone(), 1000).0.len(), 1000);
    }

    #[test]
    fn test_seed_commit_prove() {
        let seed0 = [5u8; 32];
        let seed1 = [6u8; 32];
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