//! VOLE with only two options for delta

use ff::{PrimeField, Field};
use lazy_static::lazy_static;

use crate::{Fr, vecccom::{fr_from_be_u64slice_unchecked, expand_seed_to_Fr_vec}};

lazy_static! {
    // Commented out original logic to generate these delta choices, so it can still be seen / verified
    pub static ref DELTA_CHOICES: [Fr; 2] = {
        let mut first_digest = *blake3::hash("FIRST ∆".as_bytes()).as_bytes();
        let mut second_digest = *blake3::hash("SECOND ∆".as_bytes()).as_bytes();
        // Truncate to 254 bits:
        first_digest[0] &= 0b0011_1111;
        second_digest[0] &= 0b0011_1111;

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
        [fr_from_be_u64slice_unchecked(&first_as_u64s), fr_from_be_u64slice_unchecked(&second_as_u64s)]
    };
}
pub struct ProverSmallVOLEOutputs { 
    u: Vec<Fr>,
    v: Vec<Fr>,
}
pub struct VerifierSmallVOLEOutputs {
    delta: Fr,
    q: Vec<Fr>,
}

pub struct VOLE;
impl VOLE {
    /// Creates a small VOLE from two seeds and two Deltas
    pub fn prover_outputs(seed1: &[u8; 32], seed2: &[u8; 32], vole_length: usize) -> ProverSmallVOLEOutputs {
        let out1 = expand_seed_to_Fr_vec(seed1.clone(), vole_length);
        let out2 = expand_seed_to_Fr_vec(seed2.clone(), vole_length);
        let zipped = out1.iter().zip(out2.iter());
        let u = zipped.clone().map(|(o1, o2)| *o1 + o2).collect();
        let v = zipped.map(|(o1, o2)| Fr::ZERO - (*o1 * DELTA_CHOICES[0] + *o2 * DELTA_CHOICES[1]) ).collect();
        ProverSmallVOLEOutputs { u, v }
    }
    /// Verifier shuold call this after (get) to receive their small VOLE output
    pub fn verifier_outputs(seed_i_know: &[u8; 32], idx_i_dont_know: bool, vole_length: usize) -> VerifierSmallVOLEOutputs {
        let out = expand_seed_to_Fr_vec(seed_i_know.clone(), vole_length);
        let (delta, delta_minus_other_delta) = if idx_i_dont_know { 
            (DELTA_CHOICES[1], DELTA_CHOICES[1] - DELTA_CHOICES[0])
        } else { 
            (DELTA_CHOICES[0], DELTA_CHOICES[0] - DELTA_CHOICES[1])
        };

        let q = out.iter().map(|o| *o * delta_minus_other_delta).collect();
        VerifierSmallVOLEOutputs { delta, q }
    }
    
}

#[cfg(test)]
mod test {
    use itertools::izip;

    use super::*;

    #[test]
    fn test_vole_works() {
        let seed0 = [9u8; 32];
        let seed1 = [2u8; 32];
        let prover_outputs = VOLE::prover_outputs(&seed0, &seed1, 100);
        let verifier_outputs_0 = VOLE::verifier_outputs(&seed0, true, 100);
        let verifier_outputs_1 = VOLE::verifier_outputs(&seed1, false, 100);

        assert!(
            izip!(&prover_outputs.u, &prover_outputs.v, &verifier_outputs_0.q).all(
                |(u, v, q)| u.clone() * verifier_outputs_0.delta + v == q.clone()
            )
        );

        assert!(
            izip!(prover_outputs.u, prover_outputs.v, verifier_outputs_1.q).all(
                |(u, v, q)| u * verifier_outputs_1.delta + v == q
            )
        )
    }
}
