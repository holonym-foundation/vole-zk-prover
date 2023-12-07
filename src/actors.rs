///! Provides the prover and verifier structs
pub mod actors {
    use anyhow::{Error, anyhow, Ok};
    use ff::{PrimeField, Field};
    use rand::{SeedableRng, rngs::{StdRng, ThreadRng}, RngCore};

    use crate::{subspacevole::{RAAACode, calc_consistency_check, LinearCode}, FrVec, FrMatrix, Fr, zkp::{R1CS, quicksilver::{ZKP, self}, R1CSWithMetadata}, vecccom::{expand_seed_to_Fr_vec, commit_seeds, commit_seed_commitments, proof_for_revealed_seed, reconstruct_commitment}, utils::{truncate_u8_32_to_254_bit_u64s_be, rejection_sample_u8s}, smallvole::{ProverSmallVOLEOutputs, self, DELTA_CHOICES, VOLE}, ScalarMul, challenges::{challenge_from_seed, calc_deltas, calc_quicksilver_challenge}, NUM_VOLES};



pub struct Prover {
    pub code: RAAACode,
    pub vole_length: usize,
    pub num_voles: usize,
    pub witness: FrMatrix,
    /// Commitment to the witness set after the prover makes the subspace VOLE
    pub witness_comm: Option<FrMatrix>,
    pub circuit: R1CSWithMetadata,
    /// Starts as None, added when the prover makes the subsapce VOLE
    pub subspace_vole_secrets: Option<SubspaceVOLESecrets>,
    /// Starts as None, added when the prover makes the subsapce VOLE
    pub seed_commitment: Option<[u8; 32]>
}
pub struct Verifier {
    pub circuit: R1CSWithMetadata,
    pub code: RAAACode, 
    pub num_voles: usize,
    pub vole_length: usize,
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
    /// l x k Witness split into vectors of the same length as the code's dimension k and commmited by substracting them from the first l rows of u1
    pub witness_comm: FrMatrix,
    pub subspace_vole_correction: FrMatrix,
    /// subsapce VOLE consistency check of U and V's check values, respectively
    pub consistency_check: (FrVec, FrVec)
}

pub struct Proof {
    pub zkp: ZKP,
    // pub prover_commitment: ProverCommitment,
    /// Opening of the seeds the verifier needs for subspace VOLE
    pub seed_openings: SubspaceVOLEOpening,
    /// Public input and output (u, v) tuples
    pub public_openings: PublicOpenings,
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
    /// Pads a witness and circuit to dimensions compatible with VitH and the linear code, then creates a prover
    /// Witness of length w is padded to length l where l is a multiple of a linear code's input length. creates a VOLE of length 2l+2
    /// Mutates and destroys its inputs by padding them and taking ownership of them
    pub fn from_witness_and_circuit_unpadded(mut witness: FrVec, mut circuit: R1CSWithMetadata) -> Self {
        let code = RAAACode::rand_default();
        let k = code.k();
        let pp = circuit.r1cs.calc_padding_needed(k);
        witness.zero_pad(pp.pad_len);
        circuit.r1cs.zero_pad(pp.pad_len);

        let mut witness_rows = Vec::with_capacity(pp.num_padded_wtns_rows);
        let mut start_idx = 0;
        for _i in 0..pp.num_padded_wtns_rows {
            witness_rows.push(FrVec(
                witness.0.get(start_idx .. start_idx + k).expect("This panic should not be reached").to_vec()
            ));
            start_idx += k;
        }

        Self {
            num_voles: code.n(),
            /// One extra row for the hiding of the linear combination of the relevant values in the consistency check
            /// 2x extra rows to convert subsapce VOLE into VitH. Overall, we require 2 * `num_padded_witness_rows` + 2 rows
            vole_length: 2 * (pp.num_padded_wtns_rows + 1),
            code,
            circuit,
            witness: FrMatrix(witness_rows),
            seed_commitment: None,
            subspace_vole_secrets: None,
            witness_comm: None,
        }
    }
    /// Called first
    /// Mutates self to contain secret artifacts, returning a commitment
    // THOROUGHLY CHECK AND TEST IT GETS THE DIMENSIONS OF U, V, U1, U2, V1, V2, WITNESS, ETC. CORRECT
    pub fn mkvole(&mut self) -> Result<ProverCommitment, Error> {
        if self.num_voles < 1024 { eprintln!("Less than 1024 VOLEs could result in <128 bits of soundness with current parameters for linear codes"); }
        let mut rng = ThreadRng::default();
        let mut seeds: Vec<[[u8; 32]; 2]> = vec![[[0u8; 32]; 2]; self.num_voles];
        let mut seed_commitments = Vec::with_capacity(self.num_voles);
        let mut vole_outputs = Vec::with_capacity(self.num_voles);
        for i in 0..self.num_voles {
            rng.fill_bytes(&mut seeds[i][0]);
            rng.fill_bytes(&mut seeds[i][1]);
            seed_commitments.push(commit_seeds(&seeds[i][0], &seeds[i][1]));
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

        self.witness_comm = Some(witness_comm.clone());
        if self.num_voles % self.code.q != 0 { return Err(anyhow!("invalid num_voles param")) };
        let challenge_hash = challenge_from_seed(&seed_comm, "vole_consistency_check".as_bytes(), self.vole_length);
        let consistency_check = calc_consistency_check(&challenge_hash, &new_u_rows.transpose(), &v_cols);
        


        // Before storing the secrets, split them in half, which will make reteiving the individual halves easier
        
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

        self.seed_commitment = Some(seed_comm.clone());
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
            consistency_check,
            subspace_vole_correction: correction
        })
    }
    // /// Called as part of proof()
    // /// Creates a ZKP that its multiplications are all correct
    // // /// Also calculates and stores vith_delta when done
    // fn zkp(&self) -> ZKP {
    //     todo!()
    // }

    /// Called as part of proof()
    /// Calculates the S matrix to reveal to the verifier once it learns ∆'
    /// Returns none if it lacks 
    fn s_matrix(&self, vith_delta: &Fr) -> Result<FrMatrix, Error> {
        let svs = self.subspace_vole_secrets.as_ref().ok_or(anyhow!("VOLE must be completed before this step"))?;
        Ok(&svs.u1.scalar_mul(vith_delta) + &svs.u2)
    }

    /// Wrapper for all other prover functions
    pub fn prove(&mut self) -> Result<Proof, Error> {
        let err_uncompleted = ||anyhow!("VOLE must be completed before this step");
        let svs = self.subspace_vole_secrets.as_ref().ok_or(err_uncompleted())?;
        let seed_comm = self.seed_commitment.as_ref().ok_or(err_uncompleted())?;
        let witness_comm = self.witness_comm.as_ref().ok_or(err_uncompleted())?;
        // TODO: without so much cloning
        let prover = quicksilver::Prover::from_vith(svs.u1.clone(), svs.u2.clone(), self.witness.clone(), self.circuit.clone());
        let challenge = calc_quicksilver_challenge(seed_comm, &witness_comm);
        let zkp = prover.prove(&challenge); 

        let public_openings = PublicOpenings {
            public_inputs: prover.open_public(&self.circuit.public_inputs_indices),
            public_outputs: prover.open_public(&self.circuit.public_outputs_indices)
        };
        
        let (subspace_deltas, vith_delta) = calc_deltas(seed_comm, &zkp, self.num_voles, &public_openings);
        let s_matrix = self.s_matrix(&vith_delta)?;

        let mut openings = Vec::with_capacity(self.num_voles);
        let mut opening_proofs = Vec::with_capacity(self.num_voles);
        for i in 0..svs.seeds.len() {
            openings.push(svs.seeds[i][subspace_deltas[i]]);
            opening_proofs.push(
                proof_for_revealed_seed(&svs.seeds[i][1 - subspace_deltas[i]])
            );
        };

        Ok(
            Proof { 
                zkp, 
                s_matrix,
                public_openings,
                seed_openings : SubspaceVOLEOpening { seed_opens: openings, seed_proofs: opening_proofs },
            }
        )
    }

}

impl Verifier {
    /// Calculates the dimensions of the vole and pads the circuit. 
    pub fn from_circuit(mut circuit: R1CSWithMetadata) -> Self {
        let code = RAAACode::rand_default();
        let pp = circuit.r1cs.calc_padding_needed(code.k());
        circuit.r1cs.zero_pad(pp.pad_len);
        Verifier { 
            circuit, 
            num_voles: code.n(), 
            /// One extra row for the hiding of the linear combination of the relevant values in the consistency check
            /// 2x extra rows to convert subsapce VOLE into VitH. Overall, we require 2 * `num_padded_witness_rows` + 2 rows
            vole_length: 2 * (pp.num_padded_wtns_rows + 1),
            code, 
            subspace_vole_deltas: None, 
            vith_delta: None 
        }
        

        
    //     pub circuit: R1CSWithMetadata,
    // pub code: RAAACode, 
    // pub num_voles: usize,
    // pub vole_length: usize,
    // /// Starts as None, set during Fiat Shamir
    // pub subspace_vole_deltas: Option<FrVec>,
    // /// Starts as None, set during Fiat Shamir
    // pub vith_delta: Option<Fr>
    }

    /// TODO: ensure every value in the ProverCommitment and Proof is checked in some way by this function:
    pub fn verify(&self, comm: &ProverCommitment, proof: &Proof) -> Result<PublicOpenings, Error> {
        let (delta_choices, vith_delta) = calc_deltas(&comm.seed_comm, &proof.zkp, self.num_voles, &proof.public_openings);
        let mut deltas = Vec::<Fr>::with_capacity(self.num_voles);
        let mut q = Vec::<FrVec>::with_capacity(self.num_voles);
        // Calculate small VOLE outputs then check they were all commited to in comm.seed_comm 
        let mut hasher = blake3::Hasher::new();
        for i in 0..self.num_voles {
            let rec = reconstruct_commitment(
                &proof.seed_openings.seed_opens[i], 
                delta_choices[i] != 0, // Convert usize that should be 0 or 1 to bool
                &proof.seed_openings.seed_proofs[i]
            );
            hasher.update(&rec);
            let vole_outs = VOLE::verifier_outputs(&proof.seed_openings.seed_opens[i], delta_choices[i] == 0, self.vole_length);
            deltas.push(vole_outs.delta);
            q.push(vole_outs.q);
            
        }

        if !(*hasher.finalize().as_bytes() == comm.seed_comm) { return Err(anyhow!("Seed commitment is not a commitment to the seeds")) }
        
        // Construct the subspace VOLE
        let deltas = &FrVec(deltas);
        let new_q = self.code.correct_verifier_qs(&FrMatrix(q).transpose(), deltas, &comm.subspace_vole_correction);
        // Check that its outputs are in the subspace 
        let challenge_hash = &challenge_from_seed(&comm.seed_comm, "vole_consistency_check".as_bytes(), self.vole_length);
        self.code.verify_consistency_check(challenge_hash, &comm.consistency_check, deltas, &new_q)?;
        // Check S matrix is in the subspace as well 
        // IF TESTS FAIL CHECK THESE ARE ROWS VS. COLUMNS
        self.code.check_parity_batch(&proof.s_matrix.0)?;
        // Verify the ZKP
        let zk_verifier = quicksilver::Verifier::from_vith(new_q.transpose(), vith_delta.clone(), self.circuit.clone());
        let quicksilver_challenge = calc_quicksilver_challenge(&comm.seed_comm, &comm.witness_comm);
        zk_verifier.verify(&quicksilver_challenge, &proof.zkp)?;
        zk_verifier.verify_public(&proof.public_openings)?;
        Ok(proof.public_openings.clone())
    }

    
}

/// Values of the witness that the prover opens
#[derive(Clone, Debug)]
pub struct PublicOpenings {
    pub public_inputs: Vec<(Fr, Fr)>,
    pub public_outputs: Vec<(Fr, Fr)>
}

}

#[cfg(test)]
mod test {
    use ff::PrimeField;
    use crate::{subspacevole::RAAACode, zkp, actors::actors::{Prover, Verifier}, Fr, FrVec};
    use super::*;

    #[test]
    fn prover_verifier_full_integration_tiny_circuit() {
        let mut circuit = zkp::test::TEST_R1CS_WITH_METADA.clone();
        let mut correct_witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());
        let mut prover = Prover::from_witness_and_circuit_unpadded(correct_witness, circuit.clone());
        let vole_comm = prover.mkvole().unwrap();
        let proof = prover.prove().unwrap();


        let verifier = Verifier::from_circuit(circuit);
        let result = verifier.verify(&vole_comm, &proof);
        println!("result: {:?}", result);
        todo!()

    }

    #[test]
    fn prover_verifier_full_integration_circuit_gt_1024_constraints() {
        let circuit = zkp::test::TEST_R1CS_WITH_METADA.clone();
        todo!()

    }
}