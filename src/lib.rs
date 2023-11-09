use std::str::FromStr;

use ff::Field;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use sha2::{Digest,Sha512};

#[macro_use]
extern crate ff;

use crate::ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fr([u64; 4]);

lazy_static! {
    pub static ref FIELD_MODULUS: BigUint = BigUint::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495617").unwrap();
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
/// Returns 2^log_length Fr elements
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
    // Collect the children as a Vec<Fr>
    recent_lvl.iter().map(|x| {
        let mut b = BigUint::from_bytes_le(x.as_ref());
        while b > *FIELD_MODULUS {
            b = BigUint::from_bytes_le(blake3::hash(&b.to_bytes_le()).as_bytes());
        }
        // rejection sampling
        Fr::from_str_vartime(&b.to_string()).unwrap()
    }).collect()
}
#[cfg(test)]
mod test {
    #[test]
    fn test_seed_expansion_len() {
        let seed = [0u8; 32];
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 0).len(), 1);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 1).len(), 2);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 2).len(), 4);
        assert_eq!(super::expand_seed_to_Fr_vec(seed.to_vec(), 10).len(), 1024);
    }
}