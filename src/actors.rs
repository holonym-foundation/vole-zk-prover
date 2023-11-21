///! Provides the prover and verifier structs
mod actors {
    use anyhow::Error;
    use ff::PrimeField;
    use rand::{SeedableRng, rngs::StdRng};

    use crate::{subspacevole::RAAACode, FrVec, FrMatrix, Fr, zkp::R1CS, vecccom::{u64s_overflow_field, unchecked_fr_from_be_u64_slice}};


pub struct Prover {
    pub code: RAAACode,
    pub vole_length: usize,
    pub num_voles: usize,
    pub witness: FrVec,
    pub circuit: R1CS,
    /// Starts as None, added when the prover makes the VOLE
    pub secret_artifcats: Option<ProverSecretArtifacts>,
}
pub struct Verifier {
    pub code: RAAACode, 
    /// Starts as None, set during Fiat Shamir
    pub subspace_vole_deltas: Option<FrVec>,
    /// Starts as None, set during Fiat Shamir
    pub vith_delta: Option<Fr>
}

/// Seeds, U, V, etc. Anything the prover made and needs to keep to itself
pub struct ProverSecretArtifacts {
    seeds: Vec<[u8; 32]>,
    u: FrMatrix,
    v: FrMatrix,
}
pub struct ProverCommitment {
    /// Hash of every pair of seed's respective hashes for the seeds used to create the VOLEs. We are just using two seeds per VOLE!
    /// Can/should be used for Fiat-Shamir of subspace VOLE consistency check
    pub seed_comm: [u8; 32],
    /// Witness split into vectors of the same length as the code's dimension
    pub witness_comm: FrMatrix,
    /// subsapce VOLE consistency check of U and V's check values, respectively
    pub consistency_check: (FrVec, FrVec)
}
pub struct ZKP {
    /// Quicksilver multiplication proof of one field element
    pub mul_proof: Fr,
    /// Values for the final gate in the ZKP
    pub last_gate_opening: (Fr, Fr)
}

pub struct Proof {
    pub zkp: ZKP,
    pub prover_commitment: ProverCommitment,
    pub prover_opening: SubspaceVOLEOpening,
    /// The VitH S matrix 
    pub s_matrix: FrMatrix,
}

pub struct SubspaceVOLEOpening {
    /// Openings of one seed per pair
    pub seed_opens: Vec<[u8; 32]>,
    /// Proofs that the openings were done correctly
    pub seed_proofs: Vec<[u8; 32]>,
    // /// S matrix from the final VOLE in the head
    // pub vith_s: FrMatrix,
    // 
    // pub final_gate: (Fr, Fr)
}

impl Prover {
    /// Called first
    /// Mutates self to contain secret artifacts, returning a commitment
    fn mkvole(&mut self) -> ProverCommitment {
        todo!()
    }
    /// Called second
    /// Creates a ZKP that its multiplications are all correct
    fn zkp(&self) -> ZKP {
        todo!()
    }
    /// Calculates Reveals the 
    fn s_matrix(&self) -> FrMatrix {
        todo!()
    }
}

impl Verifier {
    fn check() -> Result<(), Error> {
        todo!()
    }
    /// Mutates self to include choices for delta
    fn set_deltas(&mut self, proof: &Proof) -> FrVec {
        todo!()
    }
    /// Generates a vector of length `length` from a seed (e.g. from the commitment to the prover's seeds)
    /// Be careful not to call this twice the same seed unless that is intended -- it will generate the same randomness
    fn challenge_from_seed(seed: <StdRng as SeedableRng>::Seed, length: usize) -> FrVec {
        todo!()
    }
    fn check_zkp() -> Result<(), Error> {
        todo!()
    }
}

/// Called by Verifier and Prover to calculate the ∆' 
/// WARNING: if Proof struct is modified and this line of code is not, Fiat-Shamir could be rendered insecure !
pub fn calc_vith_delta(proof: &Proof) -> Fr {
    // Fiat-Shamir
        // TODO: double check it's fine to skip hashing the witness commitment. I am pretty confident it is:
        // if the prover changes their witness commitment, they will get caught by it either 
        // 1. not being a valid witness
        // 2. not corresponding to a valid VOLE

    let mut frs = vec![proof.zkp.last_gate_opening.0, proof.zkp.last_gate_opening.1, proof.zkp.mul_proof];
    let mut seed_comm = proof.prover_commitment.seed_comm;
    // Concatenate Frs byte representation with seed commitment
    let mut concatted = Vec::with_capacity(32 * (1 + frs.len()));
    concatted.append(&mut seed_comm.to_vec());
    frs.iter_mut().for_each(|f| {
        concatted.append(
            &mut f.to_repr().0.to_vec()
        )
    });

    let mut delta = blake3::hash(&concatted).as_bytes();
    let mut delta_u64_4 = &u8_32_to_u64_4(*delta);
    // Rejection sample
    while u64s_overflow_field(delta_u64_4) {
        delta = blake3::hash(delta).as_bytes();
        delta_u64_4 = &u8_32_to_u64_4(*delta);
    }
    
    unchecked_fr_from_be_u64_slice(delta_u64_4)

}
// Probably best in different file:
pub fn u8_32_to_u64_4(x: [u8; 32]) -> [u64; 4] {
    [
        u64::from_be_bytes(x[0..8].try_into().unwrap()), 
        u64::from_be_bytes(x[8..16].try_into().unwrap()), 
        u64::from_be_bytes(x[16..24].try_into().unwrap()), 
        u64::from_be_bytes(x[24..32].try_into().unwrap())
    ]
}
}