///! Provides the prover and verifier structs
mod actors {
    use anyhow::{Error, anyhow};
    use ff::{PrimeField, Field};
    use rand::{SeedableRng, rngs::{StdRng, ThreadRng}, RngCore};

    use crate::{subspacevole::{RAAACode, calc_consistency_check}, FrVec, FrMatrix, Fr, zkp::R1CS, vecccom::{expand_seed_to_Fr_vec, commit_seeds, commit_seed_commitments}, utils::{truncate_u8_32_to_254_bit_u64s_be, rejection_sample_u8s}, smallvole::{ProverSmallVOLEOutputs, self}, ScalarMul};


pub struct Prover {
    pub code: RAAACode,
    pub vole_length: usize,
    pub num_voles: usize,
    pub witness: FrMatrix,
    pub circuit: R1CS,
    /// Starts as None, added when the prover makes the subsapce VOLE
    pub subspace_vole_secrets: Option<SubspaceVOLESecrets>,
}
pub struct Verifier {
    pub code: RAAACode, 
    /// Starts as None, set during Fiat Shamir
    pub subspace_vole_deltas: Option<FrVec>,
    /// Starts as None, set during Fiat Shamir
    pub vith_delta: Option<Fr>
}

/// Anything that the prover has learned by the time of the subspace VOLE's completion that it must keep hidden:
pub struct SubspaceVOLESecrets {
    seeds: Vec<([[u8; 32]; 2])>,
    // u: FrMatrix,
    // v: FrMatrix,
    /// First half of u_1s rows
    u1: FrMatrix,
    /// Second half of u_1s rows
    u2: FrMatrix,
    /// First half of v1_s rows
    v1: FrMatrix,
    /// Second half of v1_s rows
    v2: FrMatrix
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
    // THOROUGHLY CHECK AND TEST IT GETS THE DIMENSIONS OF U, V, U1, U2, V1, V2, WITNESS, ETC. CORRECT
    fn mkvole(&mut self) -> Result<ProverCommitment, Error> {
        if self.num_voles < 1024 { eprintln!("Less than 1024 VOLEs could result in <128 bits of soundness with current parameters for linear codes"); }
        let mut rng = ThreadRng::default();
        let mut seeds: Vec<[[u8; 32]; 2]> = vec![[[0u8; 32]; 2]; self.num_voles];
        let mut seed_commitments = Vec::with_capacity(self.num_voles);
        let mut vole_outputs = Vec::with_capacity(self.num_voles);
        for i in 0..self.num_voles {
            rng.fill_bytes(&mut seeds[i][0]);
            rng.fill_bytes(&mut seeds[i][1]);
            seed_commitments.push(commit_seeds(&seeds[i][0], &seeds[i][0]));
            vole_outputs.push(smallvole::VOLE::prover_outputs(&seeds[i][0], &seeds[i][1], self.vole_length));
        }

        let seed_comm = commit_seed_commitments(&seed_commitments);

        let u_prime_cols = FrMatrix(vole_outputs.iter().map(|o|o.u.clone()).collect::<Vec<FrVec>>());
        let v_cols = FrMatrix(vole_outputs.iter().map(|o|o.v.clone()).collect::<Vec<FrVec>>());

        let u_prime_rows = u_prime_cols.transpose();
        let v_rows = v_cols.transpose();

        let (new_u_rows, correction) = self.code.get_prover_correction(&u_prime_rows);
        
        // If I encounter an error when testing, length/dim seems likely culprit here
        let witness_comm = &FrMatrix(new_u_rows.0[0..self.witness.0.len()].to_vec()) 
                    - 
                    &self.witness;


        if !(self.code.q % self.num_voles == 0) {return Err(anyhow!("invalid num_voles param")) };
        let challenge_hash = challenge_from_seed(seed_comm, self.vole_length);
        let consistency_check = calc_consistency_check(&challenge_hash, &new_u_rows.transpose(), &v_cols);
        


        /// Before storing the secrets, split them in half, which will make reteiving the individual halves easier
        
        let u_len = new_u_rows.0.len();
        let v_len = v_rows.0.len();

        if !(u_len % 2 == 0) { return Err(anyhow!("Number of u's rows must be even")) }
        if !(v_len % 2 == 0) { return Err(anyhow!("Number of v's rows must be even")) }

        let half_u_len = u_len / 2;
        let half_v_len = v_len / 2;

        let u1 = FrMatrix(new_u_rows.0[0..half_u_len].to_vec());
        let u2 = FrMatrix(new_u_rows.0[half_u_len..u_len].to_vec());

        let v1 = FrMatrix(v_rows.0[0..half_v_len].to_vec());
        let v2 = FrMatrix(v_rows.0[half_v_len..v_len].to_vec());


        self.subspace_vole_secrets = Some(SubspaceVOLESecrets {
            seeds,
            u1,
            u2,
            v1,
            v2

        });
        Ok(ProverCommitment { 
            seed_comm, 
            witness_comm,
            consistency_check
        })
    }
    /// Called second
    /// Creates a ZKP that its multiplications are all correct
    fn zkp(&self) -> ZKP {
        todo!()
    }
    /// Called last
    /// Calculates the S matrix to reveal to the verifier once it learns ∆'
    /// Returns none if it lacks 
    fn s_matrix(&self, vith_delta: &Fr) -> Result<FrMatrix, Error> {
        let svs = self.subspace_vole_secrets.as_ref().ok_or(anyhow!("VOLE (and ZKP) must be completed before this step"))?;
        Ok(&svs.u1.scalar_mul(vith_delta) + &svs.u2)

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
    fn check_zkp() -> Result<(), Error> {
        todo!()
    }
}

/// Generates a vector of length `length` from a seed (e.g. from the commitment to the prover's seeds)
/// Be careful not to call this twice the same seed unless that is intended -- it will generate the same randomness
fn challenge_from_seed(seed: <StdRng as SeedableRng>::Seed, length: usize) -> FrVec {
    expand_seed_to_Fr_vec(seed, length)
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
    // Concatenate Frs byte representation with seed commitment
    let mut concatted = Vec::with_capacity(32 * (1 + frs.len()));
    concatted.append(&mut proof.prover_commitment.seed_comm.to_vec());

    frs.iter_mut().for_each(|f| {
        concatted.append(
            &mut f.to_repr().0.to_vec()
        )
    });

    let delta_first_try = *blake3::hash(&concatted).as_bytes();
    rejection_sample_u8s(&delta_first_try)
}
}