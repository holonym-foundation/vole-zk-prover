//! Fiat-shamir challenges all in one place

use ff::PrimeField;
use rand::{rngs::StdRng, SeedableRng, RngCore};

use crate::{zkp::quicksilver::ZKP, FrVec, vecccom::expand_seed_to_Fr_vec, actors::actors::PublicOpenings, Fr, utils::rejection_sample_u8s, FrMatrix, DotProduct};
/// Generates a vector of length `length` from a seed (e.g. from the commitment to the prover's seeds)
/// Be careful not to call this twice the same seed unless that is intended -- it will generate the same randomness
/// Hence, the salt is included to prevent this from easily happening on accident.
pub fn challenge_from_seed(seed: &[u8], salt: &[u8], length: usize) -> FrVec {
    let seed = {
        *blake3::hash(&[seed, salt].concat()).as_bytes()
    };
    expand_seed_to_Fr_vec(seed, length)
}

pub fn calc_quicksilver_challenge(seed_comm: &[u8; 32], witness_comm: &FrMatrix) -> Fr {
    let mut concatted = &mut seed_comm.to_vec();
    concatted.append(&mut "quicksilver_challenge".as_bytes().to_vec());
    let digest = *blake3::hash(concatted).as_bytes();
    // Universal hash of witness commitment to compress it to one value
    let universal_inner = challenge_from_seed(&digest, &"quicksilver_inner".as_bytes(), witness_comm.0.len());
    let universal_outer = challenge_from_seed(&digest, &"quicksilver_outer".as_bytes(), witness_comm.0[0].0.len());
    let compressed = universal_outer.dot(
        &(&universal_inner * witness_comm)
    );

    let mut output_preimg = digest.to_vec();
    output_preimg.append(&mut compressed.to_repr().0.to_vec());
    let output_digest = *blake3::hash(&output_preimg).as_bytes();
    rejection_sample_u8s(&output_digest)
}

/// Called by Verifier and Prover to calculate the original VOLE ∆s along with the ∆' 
/// Takes seed commitment and ZKP as input
/// Returns (subfield VOLE indices, VitH choice)
pub fn calc_deltas(seed_comm: &[u8; 32], zkp: &ZKP, num_voles: usize, public_openings: &PublicOpenings) -> (Vec<usize>, Fr) {
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
        concatted.append(
            &mut f.to_repr().0.to_vec()
        )
    });

    let delta_first_try = *blake3::hash(&concatted).as_bytes();
    let vith_delta = rejection_sample_u8s(&delta_first_try);

    concatted.append(&mut "subspace_vole_challenge".as_bytes().to_vec());
    let subspace_vole_delta_seed = *blake3::hash(&concatted).as_bytes();
    let prg_seed: <StdRng as SeedableRng>::Seed = subspace_vole_delta_seed;
    let mut prg = StdRng::from_seed(prg_seed);
    let mut delta_choices: Vec<usize> = Vec::with_capacity(num_voles);
    // This is inefficient but not a bottleneck
    (0..num_voles).for_each(|_|delta_choices.push((prg.next_u32() % 2) as usize));

    (delta_choices, vith_delta)

}