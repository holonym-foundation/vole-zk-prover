//! Fiat-shamir challenges all in one place
use ff::{PrimeField, Field};
use rand::{RngCore, SeedableRng};
use rand_chacha::ChaCha12Rng;
use crate::{zkp::quicksilver::ZKP, FVec, vecccom::expand_seed_to_field_vec, actors::actors::PublicOpenings, FMatrix, DotProduct, PF};

pub struct Challenges<T: PF> {
    /// Small-field VOLE ∆ indices
    pub delta_choices: Vec<usize>,
    /// VitH ∆'
    pub vith_delta: T,
    // /// TODO, low priority: making this the whole challenge vector that's generated from the quicksilver_challenge scalar
    // pub quicksilver_challenge: Fr,
    /// Consistency check challelnge for the validity of the Q, U and V matrices
    pub subspace_challenge: FVec<T>,
    /// Consistency check challenge for the validity of the S matrix
    pub s_challenge: FVec<T>
}
/// Generates a vector of length `length` from a seed (e.g. from the commitment to the prover's seeds)
/// Be careful not to call this twice the same seed unless that is intended -- it will generate the same randomness
/// Hence, the salt is included to prevent this from easily happening on accident.
pub fn challenge_from_seed<T: PF>(seed: &[u8], salt: &[u8], length: usize) -> FVec<T> {
    let seed = {
        *blake3::hash(&[seed, salt].concat()).as_bytes()
    };
    expand_seed_to_field_vec(seed, length)
}

pub fn calc_quicksilver_challenge<T: PF>(seed_comm: &[u8; 32], witness_comm: &FMatrix<T>) -> T {
    // Universal hash of witness commitment to compress it to one value
    let universal_inner = challenge_from_seed(seed_comm, &"quicksilver_inner".as_bytes(), witness_comm.0.len());
    let universal_outer = challenge_from_seed(seed_comm, &"quicksilver_outer".as_bytes(), witness_comm.0[0].0.len());
    let compressed = universal_outer.dot(
        &(&universal_inner * witness_comm)
    );
    // Hashing may be unecessary but is cheap and removes any potential linear correlation (i have not checekd whether that correlation would be problematic)
    let digest = *blake3::hash(&compressed.to_u8s()).as_bytes();
    T::random(&mut ChaCha12Rng::from_seed(digest))
}

/// Called by Verifier and Prover to calculate the original VOLE ∆s along with the ∆' 
/// seed commitment and ZKP as input
/// Returns (subfield VOLE indices, VitH choice)
/// Important note: if u, v, q, ∆ are known to the prover, the prover can forge another (u, v) pair \
/// that satisfies q = v + u∆
/// therefore, the prover should open the public inputs before learning ∆. In Fiat-Shamir, ∆'s calculation should then include all prover ZKP and public openings
pub fn calc_other_challenges<T: PF>(seed_comm: &[u8; 32], witness_comm: &FMatrix<T>, zkp: &ZKP<T>, vole_length: usize, num_voles: usize, public_openings: &PublicOpenings<T>) -> Challenges<T> {
    // Fiat-Shamir
        // TODO: double check it's fine to skip hashing the witness commitment. I am pretty confident it is:
        // if the prover changes their witness commitment, they will get caught by it either 
        // 1. not being a valid witness
        // 2. not corresponding to a valid VOLE
    let mut frs = vec![zkp.mul_proof.0, zkp.mul_proof.1];
    for i in 0..public_openings.public_inputs.len(){
        frs.push(public_openings.public_inputs[i].0);
        frs.push(public_openings.public_inputs[i].1);

    }
    for i in 0..public_openings.public_outputs.len(){
        frs.push(public_openings.public_outputs[i].0);
        frs.push(public_openings.public_outputs[i].1);

    }
    let mut concatted = &mut seed_comm.to_vec();

    // Concatenate Frs byte representation with seed commitment
    // let mut concatted = Vec::with_capacity(32 * (1 + frs.len()));

    frs.iter_mut().for_each(|f| {
        concatted.append(&mut f.to_u8s())
    });

    let delta_first_try = *blake3::hash(&concatted).as_bytes();
    let vith_delta = T::random(&mut ChaCha12Rng::from_seed(delta_first_try));

    concatted.append(&mut "subspace_vole_challenge".as_bytes().to_vec());
    let subspace_vole_delta_seed = *blake3::hash(&concatted).as_bytes();
    let mut prg = ChaCha12Rng::from_seed(subspace_vole_delta_seed);
    let mut delta_choices: Vec<usize> = Vec::with_capacity(num_voles);
    // This is inefficient but not a bottleneck
    (0..num_voles).for_each(|_|delta_choices.push((prg.next_u32() % 2) as usize));
    
    let subspace_challenge = challenge_from_seed(&concatted, "subspace_vole_consistency".as_bytes(), vole_length);
    assert!(vole_length % 2 == 0, "VOLE length must be a multiple of 2");
    let s_challenge = challenge_from_seed(&concatted, "s_matrix_consistency".as_bytes(), vole_length / 2);
    
    Challenges {
         delta_choices,
         vith_delta, 
        //  quicksilver_challenge: calc_quicksilver_challenge(seed_comm, witness_comm), 
         subspace_challenge, 
         s_challenge
    }

}