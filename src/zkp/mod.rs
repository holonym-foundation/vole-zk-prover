use crate::{FrMatrix, Fr, FrVec};

#[derive(Clone)]
pub struct R1CS {
    a_rows: FrMatrix,
    b_rows: FrMatrix,
    c_rows: FrMatrix,
}
impl R1CS {
    /// Checks whether it is satisfiable by the witness
    fn witness_check(&self, witness: &FrVec) -> bool {
        (witness * &self.a_rows)
        *
        (witness * &self.b_rows)
        ==
        (witness * &self.c_rows)
    }
}
#[derive(Clone)]
pub struct R1CSWithMetadata {
    r1cs: R1CS,
    public_inputs_indices: Vec<usize>,
    public_outputs_indices: Vec<usize>
}

pub mod quicksilver {
    use crate::{FrVec, Fr, FrMatrix, DotProduct};

    use super::{R1CS, R1CSWithMetadata};

    #[derive(Debug)]
    pub struct ZKP {
        /// Quicksilver multiplication proof of two field elements
        pub mul_proof: (Fr, Fr),
        /// Opening (u, v) of public input wires
        pub public_input_openings: Vec<(Fr, Fr)>,
        /// Opening (u, v) of public output wires
        pub public_output_openings: Vec<(Fr, Fr)>
    }
    pub struct Prover {
        pub u: FrVec,
        pub v: FrVec,
        pub r1cs_with_metadata: R1CSWithMetadata
    }
    impl Prover {
        /// Creates a prover from VitH U1 and R matrices of equal dimension with 2l+2 rows where the witness is split into l chunks of length vole_length
        /// Takes ownership and mutates most of its inputs to something useless
        pub fn from_vith(mut u1_rows: FrMatrix, mut r_rows: FrMatrix, mut witness_rows: FrMatrix, r1cswm: R1CSWithMetadata) -> Prover {
            let r1cs = &r1cswm.r1cs;
            assert!((u1_rows.0.len() == r_rows.0.len()) && (u1_rows.0[0].0.len() == r_rows.0[0].0.len()), "u and v must be same dimension");
            assert!(witness_rows.0.len() == u1_rows.0.len() - 1, "witness must have one fewer row than u1");
            assert!(witness_rows.0[0].0.len() == u1_rows.0[0].0.len() - 1, "witness must have same number of columns as u1");
            assert!((r1cs.a_rows.0.len() == u1_rows.0.len()) && (r1cs.a_rows.0[0].0.len() == u1_rows.0[0].0.len()), "VOLE dimensions must match R1CS dimensions");

            let vith_size = u1_rows.0.len() * u1_rows.0[0].0.len();
            let mut u = Vec::with_capacity(vith_size);
            let mut v = Vec::with_capacity(vith_size);
            witness_rows.0.iter_mut().map(|row| u.append(&mut row.0));
            u.append(&mut u1_rows.0.last().unwrap().0.clone());
            r_rows.0.iter_mut().map(|row| v.append(&mut row.0));

            Self { u: FrVec(u), v: FrVec(v), r1cs_with_metadata: r1cswm }
        }
        /// TODO: explore efficiency gains for polynomial Quicksilver rather than gate-by-gate Quicksilver
        /// 
        /// 1. Calculates the outputs of linear gates, i.e. the dot product of witness with each R1CS row
        /// 2. Uses those outputs as the inputs and outputs of multiplication gates (one multiplication per R1CS row)
        /// 3. Computes and, if it is 0, returns the final gate's decommitment, + a Quicksilver multiplcation proof 
        /// NOTE: According to the Quicksilver paper, `challenge_hash` should be given after the values are determined.
        /// Think about it this way: if the prover knows `challenge` before he commits to u and v (including the witness), 
        /// The prover can find a 'collision'. This is as simple as changing the witnesss
        /// so u is different but still produces the same Quicksilver check value. Note this would not affect the underlying subspace VOLE if used with VitH since a different witness would still
        /// lay in the correct subspace. Therefore, it's important `challenge` depends on the witness.
        pub fn prove(&self, challenge: &Fr) -> ZKP {
            let l = self.u.0.len();
            let r1cs = &self.r1cs_with_metadata.r1cs;
            // Can calculate all linear gates by just dot product of the prover's values with the A, B, and C R1CS rows. These are not multiplication in & out wires
            let u_a = &self.u * &r1cs.a_rows;
            let v_a = &self.v * &r1cs.a_rows;

            let u_b = &self.u * &r1cs.b_rows;
            let v_b = &self.v * &r1cs.b_rows;

            let u_c = &self.u * &r1cs.c_rows;
            let v_c = &self.v * &r1cs.c_rows;

            // Quicksilver protocol to transform VOLE into a new VOLE for linear gates
            let new_u = &(&u_b * &v_a + &u_a * &v_b) - &v_c;
            let new_v = &v_a * &v_b;
            
            let mut challenge_vec = Vec::with_capacity(l);
            challenge_vec[0] = challenge.clone();
            for i in 1..l {
                // TODO: posisble very slight performance gain by cacheing i-1
                challenge_vec[i] = challenge_vec[i-1] * challenge;
            }
            let challenge_vec = FrVec(challenge_vec);

            ZKP {
                mul_proof: (new_u.dot(&challenge_vec), new_v.dot(&challenge_vec)),
                public_input_openings: self.r1cs_with_metadata.public_inputs_indices.iter().map(
                    |i|(self.u.0[*i], self.v.0[*i])
                ).collect(),
                public_output_openings: self.r1cs_with_metadata.public_outputs_indices.iter().map(
                    |i|(self.u.0[*i], self.v.0[*i])
                ).collect(),
            }
            

        }
    }
    pub struct Verifier {
        pub delta: Fr,
        pub q: FrVec
    }
    impl Verifier {
        /// Creates a prover from VitH U1 and R matrices of equal dimension with 2l+2 rows where the witness is split into l chunks of length vole_length
        /// Takes ownership and mutates most of its inputs to something useless
        pub fn from_vith(mut q_rows: FrMatrix, delta: Fr, r1cs: R1CS) -> Verifier {
            assert!((r1cs.a_rows.0.len() == q_rows.0.len()) && (r1cs.a_rows.0[0].0.len() == q_rows.0[0].0.len()), "VOLE dimensions must match R1CS dimensions");
            let vith_size = q_rows.0.len() * q_rows.0[0].0.len();
            let mut q = Vec::with_capacity(vith_size);
            q_rows.0.iter_mut().map(|row| q.append(&mut row.0));
            let q = FrVec(q);
            Self { delta, q }
        }
    }
}

#[cfg(test)]
mod test {
    use ff::{Field, PrimeField};
    use lazy_static::lazy_static;
    use rand::rngs::ThreadRng;
    use crate::{FrVec, Fr, ScalarMul};
    use super::{*, quicksilver::Prover};

    lazy_static! {
        pub static ref TEST_R1CS: R1CS = {
            let a_rows = vec![
                FrVec(vec![1, 1, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FrVec(vec![2, 0, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            ];
            let b_rows = vec![
                FrVec(vec![0, 2, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect())
            ];
            let c_rows = vec![
                FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FrVec(vec![0, 0, 0, 1].iter().map(|x|Fr::from_u128(*x)).collect())
            ];

            R1CS {
                a_rows: crate::FrMatrix(a_rows),
                b_rows: crate::FrMatrix(b_rows),
                c_rows: crate::FrMatrix(c_rows),
            }
        };
        pub static ref TEST_R1CS_WITH_METADA: R1CSWithMetadata = R1CSWithMetadata { 
            r1cs: TEST_R1CS.clone(), 
            public_inputs_indices: vec![0,2], 
            public_outputs_indices: vec![3] 
        };
    }
    #[test]
    fn circuit_satisfiability() {
        let witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());
        assert!(TEST_R1CS.witness_check(&witness));
        assert!(!TEST_R1CS.witness_check(&FrVec(vec![Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ONE])));
    }

    #[test]
    pub fn circuit_satisfiability_proof() {
        let witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());

        // Prove it in ZK this time:
        let delta = Fr::random(&mut ThreadRng::default());
        let v = FrVec::random(witness.0.len());
        let u = witness.clone();
        let q = &u.scalar_mul(&delta) + &v;
        
        let prover = Prover {
            u,
            v: v.clone(),
            r1cs_with_metadata: TEST_R1CS_WITH_METADA.clone()
        };
        let proof = prover.prove(&Fr::from_u128(0987654321));
        println!("proof is {:?}", proof);
        assert_eq!(proof.public_input_openings, vec![
            (witness.0[0].clone(), v.0[0].clone()),
            (witness.0[2].clone(), v.0[2].clone())
        ]);
        assert_eq!(proof.public_output_openings, vec![(witness.0[3].clone(), v.0[3].clone())]);
        todo!()
    }

    #[test]
    pub fn from_vith() {

    }
}