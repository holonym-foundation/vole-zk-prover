//! VOLE with only two options for delta

use blake3::Hasher;
use ff::{PrimeField, Field};
use lazy_static::lazy_static;
use rand::{rngs::{ThreadRng, StdRng}, SeedableRng, RngCore};

use crate::{Fr, vecccom::{expand_seed_to_Fr_vec}, FrVec, utils::{truncate_u8_32_to_254_bit_u64s_be, u64s_overflow_field, fr_from_be_u64_slice, rejection_sample_u8s}, ScalarMul};

lazy_static! {
    // Commented out original logic to generate these delta choices, so it can still be seen / verified
    pub static ref DELTA_CHOICES: [Fr; 2] = {
        let mut first_digest = *blake3::hash("First ∆".as_bytes()).as_bytes();
        let mut second_digest = *blake3::hash("Second ∆".as_bytes()).as_bytes();
        [rejection_sample_u8s(&first_digest), rejection_sample_u8s(&second_digest)]
    };
}
pub struct ProverSmallVOLEOutputs { 
    pub u: FrVec,
    pub v: FrVec,
}
pub struct VerifierSmallVOLEOutputs {
    pub delta: Fr,
    pub q: FrVec,
}

pub struct VOLE;
impl VOLE {
    /// Creates a small VOLE from two seeds and two Deltas
    pub fn prover_outputs(seed1: &[u8; 32], seed2: &[u8; 32], vole_length: usize) -> ProverSmallVOLEOutputs {
        let out1 = expand_seed_to_Fr_vec(seed1.clone(), vole_length);
        let out2 = expand_seed_to_Fr_vec(seed2.clone(), vole_length);
        let zipped = out1.0.iter().zip(out2.0.iter());
        let u = &out1 + &out2;
        let v = zipped.map(|(o1, o2)| Fr::ZERO - (*o1 * DELTA_CHOICES[0] + *o2 * DELTA_CHOICES[1]) ).collect();
        ProverSmallVOLEOutputs { u, v: FrVec(v) }
    }
    /// Verifier shuold call this after (get) to receive their small VOLE output
    pub fn verifier_outputs(seed_i_know: &[u8; 32], idx_i_dont_know: bool, vole_length: usize) -> VerifierSmallVOLEOutputs {
        let out = expand_seed_to_Fr_vec(seed_i_know.clone(), vole_length);
        let (delta, delta_minus_other_delta) = if idx_i_dont_know { 
            (DELTA_CHOICES[1], DELTA_CHOICES[1] - DELTA_CHOICES[0])
        } else { 
            (DELTA_CHOICES[0], DELTA_CHOICES[0] - DELTA_CHOICES[1])
        };

        let q = out.scalar_mul(&delta_minus_other_delta);
        VerifierSmallVOLEOutputs { delta, q }
    }
    
}
/// Construct many small VOLEs and stack into big matrix. This has both prover and verifier output in plaintext
/// TODO: Halve communication cost of sharing the seeds by bringing the seed down to 16 bytes
// #[cfg(test)]
pub struct TestMOLE {
    pub prover_commitment: [u8; 32],
    pub prover_outputs: Vec<ProverSmallVOLEOutputs>,
    pub verifier_outputs: Vec<VerifierSmallVOLEOutputs>
}
// #[cfg(test)]
impl TestMOLE {
    /// For security with the rate 1/2 RAAA code, num_voles should not be less than 1024
    pub fn init(master_seed: [u8; 32], vole_size: usize, num_voles: usize) -> TestMOLE{
        let mut hasher = Hasher::new();
        // let mut seeds = Vec::with_capacity(num_voles);
        let mut prover_outputs = Vec::with_capacity(num_voles);
        let mut verifier_outputs = Vec::with_capacity(num_voles);
        let seed_seed: <StdRng as SeedableRng>::Seed = master_seed;
        let mut r = StdRng::from_seed(seed_seed);
        for i in 0..num_voles {
            let mut seed0 = [0u8; 32];
            let mut seed1 = [0u8; 32];
            hasher.update(blake3::hash(&[seed0, seed1].concat()).as_bytes());
            r.fill_bytes(&mut seed0);
            r.fill_bytes(&mut seed1);
            prover_outputs.push(
                VOLE::prover_outputs(&seed0, &seed1, vole_size)
            );

            verifier_outputs.push(
                VOLE::verifier_outputs(&seed0, true, vole_size)
            )
            
        }
        let prover_commitment = hasher.finalize();//.as_bytes();
        Self { prover_outputs, verifier_outputs, prover_commitment: *prover_commitment.as_bytes() }
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
            izip!(&prover_outputs.u.0, &prover_outputs.v.0, &verifier_outputs_0.q.0).all(
                |(u, v, q)| u.clone() * verifier_outputs_0.delta + v == q.clone()
            )
        );

        assert!(
            izip!(prover_outputs.u.0, prover_outputs.v.0, verifier_outputs_1.q.0).all(
                |(u, v, q)| u * verifier_outputs_1.delta + v == q
            )
        )
    }
}
