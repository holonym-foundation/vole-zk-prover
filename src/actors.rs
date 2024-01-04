///! Provides the prover and verifier structs
pub mod actors {
    // use std::time::Instant;
    use anyhow::{Error, anyhow, Ok};
    use ff::{PrimeField, Field};
    use rand::{rngs::ThreadRng, RngCore};
    use serde::{Serialize, Deserialize};
    use crate::{subspacevole::{RAAACode, calc_consistency_check, LinearCode}, FrVec, FrMatrix, Fr, zkp::{R1CS, quicksilver::{ZKP, self}, R1CSWithMetadata}, vecccom::{expand_seed_to_Fr_vec, commit_seeds, commit_seed_commitments, proof_for_revealed_seed, reconstruct_commitment}, utils::{truncate_u8_32_to_254_bit_u64s_be, rejection_sample_u8s}, smallvole::{ProverSmallVOLEOutputs, self, DELTA_CHOICES, VOLE}, ScalarMul, challenges::{challenge_from_seed, calc_quicksilver_challenge, calc_other_challenges}, NUM_VOLES};



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

#[derive(Clone, Debug, Serialize, Deserialize)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Proof {
    pub zkp: ZKP,
    // pub prover_commitment: ProverCommitment,
    /// Opening of the seeds the verifier needs for subspace VOLE
    pub seed_openings: SubspaceVOLEOpening,
    /// Public input and output (u, v) tuples
    pub public_openings: PublicOpenings,
    /// The VitH S matrix 
    pub s_matrix: FrMatrix,
    /// Proof S was constructed correctly
    pub s_consistency_check: FrVec,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CommitAndProof {
    pub commitment: ProverCommitment,
    pub proof: Proof
}

#[derive(Clone, Debug, Serialize, Deserialize)]
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
        let pp = circuit.calc_padding_needed(k);

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
        
        let witness_comm = &self.witness - &FrMatrix(new_u_rows.0[0..self.witness.0.len()].to_vec());

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
            v2,
        });
        Ok(ProverCommitment { 
            seed_comm, 
            witness_comm,
            consistency_check,
            subspace_vole_correction: correction
        })
    }

    /// Called as part of proof()
    /// Calculates the S matrix to reveal to the verifier once it learns ∆' and challenge
    /// Returns (S, constency check value)
    fn s_matrix_with_consistency_proof(&self, vith_delta: &Fr, challenge: &FrVec) -> Result<(FrMatrix, FrVec), Error> {
        let svs = self.subspace_vole_secrets.as_ref().ok_or(anyhow!("VOLE must be completed before this step"))?;
        let s = &svs.u1.scalar_mul(vith_delta) + &svs.u2;
        let proof = challenge * &(&svs.v1.scalar_mul(vith_delta) + &svs.v2).transpose();
        Ok((s, proof))
    }

    /// Wrapper for all other prover functions
    pub fn prove(&mut self) -> Result<Proof, Error> {
        // let mut start = Instant::now();
        let err_uncompleted = ||anyhow!("VOLE must be completed before this step");
        let svs = self.subspace_vole_secrets.as_ref().ok_or(err_uncompleted())?;
        let seed_comm = self.seed_commitment.as_ref().ok_or(err_uncompleted())?;
        let witness_comm = self.witness_comm.as_ref().ok_or(err_uncompleted())?;

        // println!("Commited {}", start.elapsed().as_micros()); start = Instant::now();
        // TODO: without so much cloning
        let prover = quicksilver::Prover::from_vith(svs.u1.clone(), svs.u2.clone(), self.witness.clone(), self.circuit.clone());
        
        // println!("made prover from VitH {}", start.elapsed().as_micros()); start = Instant::now();

        let challenge = calc_quicksilver_challenge(seed_comm, &witness_comm);
        let zkp = prover.prove(&challenge); 

        // println!("made proof {}", start.elapsed().as_micros()); start = Instant::now();

        let public_openings = PublicOpenings {
            public_inputs: prover.open_public(&self.circuit.public_inputs_indices),
            public_outputs: prover.open_public(&self.circuit.public_outputs_indices)
        };

        // println!("made public openings {}", start.elapsed().as_micros()); start = Instant::now();

        
        let challenges = calc_other_challenges(seed_comm, witness_comm, &zkp, self.vole_length, self.num_voles, &public_openings);
        let (s_matrix, s_consistency_check) = self.s_matrix_with_consistency_proof(&challenges.vith_delta, &challenges.s_challenge)?;

        
        let mut openings = Vec::with_capacity(self.num_voles);
        let mut opening_proofs = Vec::with_capacity(self.num_voles);
        for i in 0..svs.seeds.len() {
            openings.push(svs.seeds[i][challenges.delta_choices[i]]);
            opening_proofs.push(
                proof_for_revealed_seed(&svs.seeds[i][1 - challenges.delta_choices[i]])
            );
        };
        // println!("challenges, consistency check, opening proofs: {}", start.elapsed().as_micros()); start = Instant::now();


        Ok(
            Proof { 
                zkp, 
                s_matrix,
                s_consistency_check,
                public_openings,
                seed_openings : SubspaceVOLEOpening { seed_opens: openings, seed_proofs: opening_proofs },
            }
        )
    }

    pub fn commit_and_prove(&mut self) -> Result<CommitAndProof, Error> {
        let commitment = self.mkvole()?;
        let proof = self.prove()?;
        Ok(CommitAndProof { commitment, proof })
    }

}

impl Verifier {
    /// Calculates the dimensions of the vole and pads the circuit. 
    pub fn from_circuit(mut circuit: R1CSWithMetadata) -> Self {
        let code = RAAACode::rand_default();
        let pp = circuit.calc_padding_needed(code.k());
        circuit.r1cs.zero_pad(pp.pad_len);
        Verifier { 
            circuit, 
            num_voles: code.n(), 
            // One extra row for the hiding of the linear combination of the relevant values in the consistency check
            // 2x extra rows to convert subsapce VOLE into VitH. Overall, we require 2 * `num_padded_witness_rows` + 2 rows
            vole_length: 2 * (pp.num_padded_wtns_rows + 1),
            code, 
            subspace_vole_deltas: None, 
            vith_delta: None 
        }
    }

    /// TODO: ensure every value in the ProverCommitment and Proof is checked in some way by this function:
    pub fn verify(&self, cnp: &CommitAndProof) -> Result<PublicUOpenings, Error> {
        let comm = &cnp.commitment;
        let proof = &cnp.proof;
        let challenges = calc_other_challenges(&comm.seed_comm, &comm.witness_comm, &proof.zkp, self.vole_length, self.num_voles, &proof.public_openings);
        let mut deltas = Vec::<Fr>::with_capacity(self.num_voles);
        let mut q_cols = Vec::<FrVec>::with_capacity(self.num_voles);
        // Calculate small VOLE outputs then check they were all commited to in comm.seed_comm 
        let mut hasher = blake3::Hasher::new();
        for i in 0..self.num_voles {
            let rec = reconstruct_commitment(
                &proof.seed_openings.seed_opens[i], 
                challenges.delta_choices[i] != 0, // Convert usize that should be 0 or 1 to bool
                &proof.seed_openings.seed_proofs[i]
            );
            hasher.update(&rec);
            let vole_outs = VOLE::verifier_outputs(&proof.seed_openings.seed_opens[i], challenges.delta_choices[i] == 0, self.vole_length);
            deltas.push(vole_outs.delta);
            q_cols.push(vole_outs.q);
            
        }

        if !(*hasher.finalize().as_bytes() == comm.seed_comm) { return Err(anyhow!("Seed commitment is not a commitment to the seeds")) }
        
        // Construct the subspace VOLE
        let q_rows = FrMatrix(q_cols).transpose();
        let deltas = FrVec(deltas);
        
        let new_q_rows = self.code.correct_verifier_qs(&q_rows, &deltas, &comm.subspace_vole_correction);
        // Check that its outputs are in the subspace 
        let challenge_hash = &challenge_from_seed(&comm.seed_comm, "vole_consistency_check".as_bytes(), self.vole_length);


        self.code.verify_consistency_check(challenge_hash, &comm.consistency_check, &deltas, &new_q_rows.transpose())?;

        // Perhaps this is better in a separate function since this is long but it is different to uncouple all the components of verification
        // Doing the mutability like the prover may help split large functions:
        // Check S matrix is constructed properly
        debug_assert!((new_q_rows.0.len() == self.vole_length) && (self.vole_length % 2 == 0), "Q must be vole_length and even");
        let half_len = self.vole_length / 2;
        let q1 = FrMatrix(new_q_rows.0[0..half_len].to_vec());
        let q2 = FrMatrix(new_q_rows.0[half_len..self.vole_length].to_vec());
        let sgc_diag_delta = self.code.batch_encode(&proof.s_matrix.0).iter().map( |row| row * &deltas ).collect::<Vec<FrVec>>();
        let lhs = &challenges.s_challenge * &(&q1.scalar_mul(&challenges.vith_delta) + &q2).transpose();
        let rhs = &proof.s_consistency_check + &(&challenges.s_challenge * &FrMatrix(sgc_diag_delta).transpose());
        if lhs != rhs { return Err(anyhow!("failed to verify S matrix")) }

        // Verify the ZKP
        let zk_verifier = quicksilver::Verifier::from_vith(&proof.s_matrix, challenges.vith_delta.clone(), &comm.witness_comm, self.circuit.clone());
        let quicksilver_challenge = calc_quicksilver_challenge(&comm.seed_comm, &comm.witness_comm);
        zk_verifier.verify(&quicksilver_challenge, &proof.zkp)?;
        zk_verifier.verify_public(&proof.public_openings)?;

        // Return the witness (u) values from the public openings (v isn't useful as a public value except for verifiying the proof)
        Ok(proof.public_openings.u_values())
    }

    
}

/// Values of the witness that the prover opens
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PublicOpenings {
    pub public_inputs: Vec<(Fr, Fr)>,
    pub public_outputs: Vec<(Fr, Fr)>
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PublicUOpenings {
    pub public_inputs: Vec<Fr>,
    pub public_outputs: Vec<Fr>
}
impl PublicOpenings {
    pub fn u_values(&self) -> PublicUOpenings {
        PublicUOpenings {
            public_inputs: self.public_inputs.iter().map(|(x, _)|x.clone()).collect(),
            public_outputs: self.public_outputs.iter().map(|(x, _)|x.clone()).collect()
        }
    }
}



}

pub mod test_helpers {
    use anyhow::Error;

    use crate::{FrVec, zkp::R1CSWithMetadata};

    use super::actors::{Prover, Verifier, PublicUOpenings};

    pub fn e2e_test(witness: FrVec, circuit: R1CSWithMetadata) -> Result<PublicUOpenings, Error>{
        let mut prover = Prover::from_witness_and_circuit_unpadded(witness.clone(), circuit.clone());
        // let vole_comm = prover.mkvole().unwrap();
        // let proof = prover.prove().unwrap();
        let cnp = prover.commit_and_prove().unwrap();
        let verifier = Verifier::from_circuit(circuit);
        verifier.verify(&cnp)
    }
}
#[cfg(test)]
mod test {
    use anyhow::anyhow;
    use ff::{PrimeField, Field};
    use crate::{subspacevole::RAAACode, zkp::{self, R1CSWithMetadata}, actors::{actors::{Prover, Verifier, CommitAndProof}, test_helpers::e2e_test}, Fr, FrVec, circom::witness};
    use super::*;

    
    #[test]
    fn prover_verifier_full_integration_tiny_circuit() {
        let circuit = zkp::test::TEST_R1CS_WITH_METADA.clone();
        let correct_witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());
        let len = correct_witness.0.len();

        assert!(e2e_test(correct_witness.clone(), circuit.clone()).is_ok());

        // Test every value in this small witness is accounted for (assuming it is constrained)
        for i in 0..len {
            let mut incorrect_witness = correct_witness.clone();
            incorrect_witness.0[i] += Fr::ONE;
            assert!(e2e_test(incorrect_witness, circuit.clone()).is_err());
        }
    }

    // /// This is already covered in the circom tests
    // #[test]
    // fn prover_verifier_full_integration_circuit_gt_1024_constraints() {
    //     let circuit = zkp::test::TEST_R1CS_WITH_METADA.clone();
    //     todo!()
    // }

    #[test]
    fn public_values() {
        let circuit = zkp::test::TEST_R1CS_WITH_METADA.clone();
        let witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());

        let mut prover = Prover::from_witness_and_circuit_unpadded(witness.clone(), circuit.clone());
        let vole_comm = prover.mkvole().unwrap();
        let correct_proof = prover.prove().unwrap();

        let verifier = Verifier::from_circuit(circuit);
        assert!(verifier.verify(&CommitAndProof { commitment: vole_comm.clone(), proof: correct_proof.clone() }).is_ok());
        
        // Test every value in this small array of public values is accounted for (assuming it is constrained)
        for i in 0..correct_proof.public_openings.public_inputs.len() {
            let mut incorrect_proof = correct_proof.clone();
            
            incorrect_proof.public_openings.public_inputs[i].0 += Fr::ONE;
            assert!(verifier.verify(&CommitAndProof { commitment: vole_comm.clone(), proof: incorrect_proof.clone() }).is_err());
            
            incorrect_proof = correct_proof.clone();

            incorrect_proof.public_openings.public_inputs[i].1 += Fr::ONE;
            assert!(verifier.verify(&CommitAndProof { commitment: vole_comm.clone(), proof: incorrect_proof.clone() }).is_err());
        }
        // Test every value in this small array of public values is accounted for (assuming it is constrained)
        for i in 0..correct_proof.public_openings.public_outputs.len() {
            let mut incorrect_proof = correct_proof.clone();
            
            incorrect_proof.public_openings.public_outputs[i].0 += Fr::ONE;
            assert!(verifier.verify(&CommitAndProof { commitment: vole_comm.clone(), proof: incorrect_proof.clone() }).is_err());

            incorrect_proof = correct_proof.clone();

            incorrect_proof.public_openings.public_outputs[i].0 += Fr::ONE;
            assert!(verifier.verify(&CommitAndProof { commitment: vole_comm.clone(), proof: incorrect_proof.clone() }).is_err());
        }
    }
}