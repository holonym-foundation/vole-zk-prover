//! VOLE with only two options for delta

use blake3::Hasher;
use ff::{PrimeField, Field};
use lazy_static::lazy_static;
use rand::{rngs::{ThreadRng}, SeedableRng, RngCore};
use rand_chacha::ChaCha12Rng;

use crate::{vecccom::{expand_seed_to_field_vec}, FVec, PF};

pub struct ProverSmallVOLEOutputs<T: PF> { 
    pub u: FVec<T>,
    pub v: FVec<T>,
}
pub struct VerifierSmallVOLEOutputs<T: PF> {
    pub delta: T,
    pub q: FVec<T>,
}

pub struct VOLE<T: PF> {
    delta_choices: [T; 2]
}
impl<T: PF> VOLE<T> {
    /// Generate uniformly random (not sure that's necessary but it's nice) delta choices in a small subset
    pub fn init() -> Self {
        let mut digest = *blake3::hash("Silk ∆ choices".as_bytes()).as_bytes();
        let mut rng = ChaCha12Rng::from_seed(digest);
        let f1 = T::random(&mut rng); let f2 = T::random(&mut rng);
        debug_assert!(f1 != f2);
        Self { delta_choices: [f1, f2] }
    }
    /// Creates a small VOLE from two seeds and two Deltas
    pub fn prover_outputs(&self, seed1: &[u8; 32], seed2: &[u8; 32], vole_length: usize) -> ProverSmallVOLEOutputs<T> {
        let out1 = expand_seed_to_field_vec(seed1.clone(), vole_length);
        let out2 = expand_seed_to_field_vec(seed2.clone(), vole_length);
        let zipped = out1.0.iter().zip(out2.0.iter());
        let u = &out1 + &out2;
        let v = zipped.map(|(o1, o2)| T::ZERO - (*o1 * self.delta_choices[0] + *o2 * self.delta_choices[1]) ).collect();
        ProverSmallVOLEOutputs { u, v: FVec(v) }
    }
    /// Verifier shuold call this after (get) to receive their small VOLE output
    pub fn verifier_outputs(&self, seed_i_know: &[u8; 32], idx_i_dont_know: bool, vole_length: usize) -> VerifierSmallVOLEOutputs<T> {
        let out = expand_seed_to_field_vec::<T>(seed_i_know.clone(), vole_length);
        let (delta, delta_minus_other_delta) = if idx_i_dont_know { 
            (self.delta_choices[1], self.delta_choices[1] - self.delta_choices[0])
        } else { 
            (self.delta_choices[0], self.delta_choices[0] - self.delta_choices[1])
        };

        let q = out.scalar_mul(delta_minus_other_delta);
        VerifierSmallVOLEOutputs { delta, q }
    }
    
}
/// Construct many small VOLEs and stack into big matrix. This has both prover and verifier output in plaintext
/// TODO: Halve communication cost of sharing the seeds by bringing the seed down to 16 bytes
// #[cfg(test)]
pub struct TestMOLE<T: PF> {
    pub prover_commitment: [u8; 32],
    pub prover_outputs: Vec<ProverSmallVOLEOutputs<T>>,
    pub verifier_outputs: Vec<VerifierSmallVOLEOutputs<T>>
}
// #[cfg(test)]
impl<T: PF> TestMOLE<T> {
    /// For security with the rate 1/2 RAAA code, num_voles should not be less than 1024
    pub fn init(master_seed: [u8; 32], vole_size: usize, num_voles: usize) -> TestMOLE<T> {
        let mut hasher = Hasher::new();
        // let mut seeds = Vec::with_capacity(num_voles);
        let mut prover_outputs = Vec::with_capacity(num_voles);
        let mut verifier_outputs = Vec::with_capacity(num_voles);
        let mut r = ChaCha12Rng::from_seed(master_seed);
        let vole = VOLE::init();
        for i in 0..num_voles {
            let mut seed0 = [0u8; 32];
            let mut seed1 = [0u8; 32];
            hasher.update(blake3::hash(&[seed0, seed1].concat()).as_bytes());
            r.fill_bytes(&mut seed0);
            r.fill_bytes(&mut seed1);
            prover_outputs.push(
                vole.prover_outputs(&seed0, &seed1, vole_size)
            );

            verifier_outputs.push(
                vole.verifier_outputs(&seed0, true, vole_size)
            )
            
        }
        let prover_commitment = hasher.finalize();//.as_bytes();
        Self { prover_outputs, verifier_outputs, prover_commitment: *prover_commitment.as_bytes() }
    }
}

#[cfg(test)]
mod test {
    use itertools::izip;

    use crate::Fr;

    use super::*;

    #[test]
    fn test_vole_works() {
        let seed0 = [9u8; 32];
        let seed1 = [2u8; 32];
        let vole = VOLE::<Fr>::init();
        let prover_outputs = vole.prover_outputs(&seed0, &seed1, 100);
        let verifier_outputs_0 = vole.verifier_outputs(&seed0, true, 100);
        let verifier_outputs_1 = vole.verifier_outputs(&seed1, false, 100);

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
