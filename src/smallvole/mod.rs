//! VOLE with only two options for delta

use ff::{PrimeField, Field};
use lazy_static::lazy_static;
use num_bigint::BigInt;
use poseidon_rs::Poseidon;

use crate::{Fr, vecccom::fr_from_be_u64slice};

lazy_static! {
    // Commented out original logic to generate these delta choices, so it can still be seen / verified
    pub static ref DELTA_CHOICES: [Fr; 2] = {
        let mut first_digest = *blake3::hash("FIRST ∆".as_bytes()).as_bytes();
        let mut second_digest = *blake3::hash("SECOND ∆".as_bytes()).as_bytes();
        // Truncate to 254 bits:
        first_digest[0] &= 0b0011_1111;
        second_digest[0] &= 0b0011_1111;
        println!("First and second digest: {:?} {:?}", first_digest, second_digest);
        let first_as_u64s = [
            u64::from_le_bytes(first_digest[0..8].try_into().unwrap()),
            u64::from_le_bytes(first_digest[8..16].try_into().unwrap()),
            u64::from_le_bytes(first_digest[16..24].try_into().unwrap()),
            u64::from_le_bytes(first_digest[24..32].try_into().unwrap()),
        ];
        let second_as_u64s = [
            u64::from_le_bytes(second_digest[0..8].try_into().unwrap()),
            u64::from_le_bytes(second_digest[8..16].try_into().unwrap()),
            u64::from_le_bytes(second_digest[16..24].try_into().unwrap()),
            u64::from_le_bytes(second_digest[24..32].try_into().unwrap()),
        ];
        [fr_from_be_u64slice(&first_as_u64s), fr_from_be_u64slice(&second_as_u64s)]
    };
}

#[cfg(test)]
mod test {
    use super::DELTA_CHOICES;

    #[test]
    fn test_vole_works() {
        let x = &*DELTA_CHOICES;
        println!("DELTA_CHOICES: {:?}", x)
    }
}
