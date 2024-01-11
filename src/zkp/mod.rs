use std::ops::Mul;
use blake3::hash;
use serde::{Serialize, Deserialize};
use crate::{FMatrix, FVec, NUM_VOLES, SparseVec, SparseFMatrix, PF};
#[derive(Clone, Serialize, Deserialize)]
pub struct FullR1CS<T: PF> {
    pub a_rows: FMatrix<T>,
    pub b_rows: FMatrix<T>,
    pub c_rows: FMatrix<T>,
}
#[derive(Clone, Serialize, Deserialize)]
pub struct SparseR1CS<T: PF> {
    pub a_rows: SparseFMatrix<T>,
    pub b_rows: SparseFMatrix<T>,
    pub c_rows: SparseFMatrix<T>
}
#[derive(Clone, Serialize, Deserialize)]
pub enum R1CS<T: PF> {
    Sparse(SparseR1CS<T>),
    Full(FullR1CS<T>)
}

#[derive(Clone, Serialize, Deserialize)]
pub struct R1CSWithMetadata<T: PF> {
    pub r1cs: R1CS<T>,
    pub public_inputs_indices: Vec<usize>,
    pub public_outputs_indices: Vec<usize>,
    pub unpadded_wtns_len: usize
}
#[derive(Debug)]
pub struct PadParams {
    pub orig_wtns_len: usize,
    pub padded_wtns_len: usize,
    /// How long the actual padding is, i.e. `padded_wtns_len` - `orig_wtns_len`
    pub pad_len: usize,
    /// When a matrix is formed via chunking the padded witness, this value represents how many rows it has
    pub num_padded_wtns_rows: usize
}
impl<T: PF> R1CS<T> {
    /// Returns Av, Bv, Cv for a vector v
    fn vec_mul(&self, v: &FVec<T>) -> (FVec<T>, FVec<T>, FVec<T>) {
        match self {
            Self::Sparse(s) => {
                (v * &s.a_rows, v * &s.b_rows, v * &s.c_rows)
            }
            Self::Full(f) => {
                (v * &f.a_rows, v * &f.b_rows, v * &f.c_rows)
            }
        }
    }

    /// Checks whether it is satisfiable by the witness
    fn witness_check(&self, witness: &FVec<T>) -> bool {
        let (wA, wB, wC) = self.vec_mul(witness);
        &wA * &wB == wC
    }

    pub fn zero_pad(&mut self, pad_len: usize) {
        match self {
            Self::Full(f) => {
                for i in 0..f.a_rows.0.len() {
                    f.a_rows.0[i].zero_pad(pad_len);
                    f.b_rows.0[i].zero_pad(pad_len);
                    f.c_rows.0[i].zero_pad(pad_len);
                }
            },
            Self::Sparse(_) => { 
                // Nothing needed since sprase matrix ignores zero values
            }
        }
        
    }
}

impl<T: PF> R1CSWithMetadata<T> {
    /// Given self and number of desired columns i.e. linear code `k`, returns the amount of padding required
    pub fn calc_padding_needed(&self, k: usize) -> PadParams {
        // Pad witness so its length is a product of NUM_VOLES
        // note this pads with a whole new row if it is a product. this is neither intentional nor important.
        let pad_len = k - (self.unpadded_wtns_len % k);
        let padded_len = self.unpadded_wtns_len + pad_len;
        debug_assert_eq!(padded_len % k, 0);

        let num_padded_wtns_rows = padded_len / k; 
        
        PadParams {
            orig_wtns_len: self.unpadded_wtns_len,
            padded_wtns_len: padded_len,
            pad_len,
            num_padded_wtns_rows,
        }
    }
    // pub fn circuit_id(&self) -> Result<[u8; 32], anyhow::Error> {
    //     let serialized = bincode::serialize(&self)?;
    //     let hashed = blake3::hash(&serialized);
    //     Ok(*hashed.as_bytes())
    // }
}
pub mod quicksilver {
    use itertools::Itertools;
    // use std::time::Instant;
    use serde::{Serialize, Deserialize};
    use anyhow::{Error, anyhow, bail, Ok};
    use ff::{Field, PrimeField};

    use crate::{FVec, FMatrix, DotProduct, actors::actors::PublicOpenings, PF};

    use super::{R1CS, R1CSWithMetadata};

    #[derive(Clone, Debug, Serialize, Deserialize)]
    pub struct ZKP<T: PF> {
        /// Quicksilver multiplication proof of two field elements
        pub mul_proof: (T, T),
        // Public inputs and outputs should not be checked in the Quicksilver; they should be opened after conveting VitH to subspace VOLE, before VitH ∆ is chosen 
        // It may be possible to securely reveal public inputs after ∆ is known, but why worry about it if we can reveal public inputs before cheating is as big a concern?
        // /// Opening (u, v) of public input wires
        // pub public_input_openings: Vec<(T, T)>,
        // /// Opening (u, v) of public output wires
        // pub public_output_openings: Vec<(T, T)>
    }
    pub struct Prover<T: PF> {
        pub u: FVec<T>,
        pub v: FVec<T>,
        pub r1cs_with_metadata: R1CSWithMetadata<T>
    }
    impl<T: PF> Prover<T> {
        /// Creates a prover Tom VitH U1 and R matrices of equal dimension with 2l+2 rows where the witness is split into l chunks of length vole_length
        /// Takes ownership and mutates most of its inputs to something useless
        pub fn from_vith(u1_rows: FMatrix<T>, mut r_rows: FMatrix<T>, mut witness_rows: FMatrix<T>, r1cswm: R1CSWithMetadata<T>) -> Self {
            // println!("VOLE dimensions: {:?}", (u1_rows.0.len(), u1_rows.0[0].0.len()));
            // println!("R1CS dimensions: {:?}", (r1cs.a_rows.0.len(), r1cs.a_rows.0[0].0.len()));
            assert!((u1_rows.0.len() == r_rows.0.len()) && (u1_rows.0[0].0.len() == r_rows.0[0].0.len()), "u and v must be same dimension");
            assert!(witness_rows.0.len() == u1_rows.0.len() - 1, "witness must have one fewer column than u1");
            assert!(witness_rows.0[0].0.len() == u1_rows.0[0].0.len(), "witness must same number rows as u1 does");
            // assert!((r1cs.a_rows.0.len() == u1_rows.0.len()) && (r1cs.a_rows.0[0].0.len() == u1_rows.0[0].0.len()), "VOLE dimensions must match R1CS dimensions");

            let vith_size = u1_rows.0.len() * u1_rows.0[0].0.len();
            let mut u = Vec::with_capacity(vith_size);
            let mut v = Vec::with_capacity(vith_size);
            witness_rows.0.iter_mut().for_each(|row| u.append(&mut row.0));
            // Append the final u1_row which wasn't included by iterating through the witness rows:
            u.append(&mut u1_rows.0.last().unwrap().0.clone());
            r_rows.0.iter_mut().for_each(|row| v.append(&mut row.0));
            Self { u: FVec(u), v: FVec(v), r1cs_with_metadata: r1cswm }
        }
        /// TODO: explore efficiency gains for polynomial Quicksilver rather than gate-by-gate Quicksilver
        /// 
        /// 1. Calculates the outputs of linear gates, i.e. the dot product of witness with each R1CS row
        /// 2. Uses those outputs as the inputs and outputs of multiplication gates (one multiplication per R1CS row)
        /// 3. Computes and, if it is 0, returns the final gate's decommitment, + a Quicksilver multiplcation proof 
        /// NOTE: According to the Quicksilver paper, `challenge` should be given after the VOLE outputs are determined.
        /// NOTE: Even after these values are determined, there is still some type of malleability --
        /// the prover can find a 'collision'. This is as simple as changing the witnesss
        /// so u is different but still produces the same Quicksilver check value. Note this would not affect the underlying subspace VOLE if used with VitH since a different witness would still
        /// lay in the correct subspace. Therefore, it's important `challenge` depends on the witness.
        pub fn prove(&self, challenge: &T) -> ZKP<T> {
            let l = self.u.0.len();
            let r1cs = &self.r1cs_with_metadata.r1cs;
            // let mut start = Instant::now();

            // Can calculate all linear gates by just dot product of the prover's values with the A, B, and C R1CS rows. These are not multiplication in & out wires
            let (u_a, u_b, u_c) = r1cs.vec_mul(&self.u);
            let (v_a, v_b, v_c) = r1cs.vec_mul(&self.v);
            // println!("QuickSilver Linear gates {}", start.elapsed().as_micros()); start = Instant::now();
            // Quicksilver protocol to transform VOLE into a new VOLE for linear gates
            let new_u = &(&u_b * &v_a + &u_a * &v_b) - &v_c;
            let new_v = &v_a * &v_b;
            // println!("QuickSilver Transformation {}", start.elapsed().as_micros()); start = Instant::now();
            let challenge_vec = get_challenge_vec::<T>(challenge, l);
            // println!("QuickSilver Challenge {}", start.elapsed().as_micros()); start = Instant::now();
            let mul_proof = (new_u.dot(&challenge_vec), new_v.dot(&challenge_vec));
            // println!("QuickSilver Multiplciation proof {}", start.elapsed().as_micros()); start = Instant::now();

            ZKP { mul_proof }
            

        }
        /// Opens VOLE correlations at public indices
        pub fn open_public(&self, indices: &Vec<usize>) -> Vec<(T, T)>{
            indices.iter().map(|i|{
                ( self.u.0[*i], self.v.0[*i] )
            }).collect()
        }
     }

    /// Creates a vector [challenge, challenge^2, challenge^3, ..., challenge^length]
    fn get_challenge_vec<T: PF>(challenge: &T, length: usize) -> FVec<T> {
        let mut challenge_vec = Vec::with_capacity(length);
            challenge_vec.push(challenge.clone());
            for i in 1..length {
                // TODO: posisble very slight performance gain by cacheing i-1
                challenge_vec.push(challenge_vec[i-1] * challenge);
            }
        FVec::<T>(challenge_vec)
    }
    pub struct Verifier<T: PF> {
        pub delta: T,
        pub q: FVec<T>,
        pub r1cs_with_metadata: R1CSWithMetadata<T>
    }
    impl<T: PF> Verifier<T> {
        /// Creates a verifier Tom VitH S and D matrices where D is the prover's commitment to the witness
        /// Takes ownership and mutates most of its inputs to something useless
        pub fn from_vith(s_rows: &FMatrix<T>, delta: T, witness_comm: &FMatrix<T>, r1cswm: R1CSWithMetadata<T>) -> Self {
            // Adjust S by adding the witness to its first part
            let mut s_adjustment = witness_comm.scalar_mul(delta);
            let row_len = s_adjustment.0[0].0.len();
            // Performance note: this pushes one beyond capacity
            s_adjustment.0.push(FVec::<T>(vec![T::ZERO; row_len]));

            let mut s_adjusted = s_rows + &s_adjustment;
            
            // assert!((r1cs.a_rows.0.len() == s_adjusted.0.len()) /* && (r1cs.a_rows.0[0].0.len() == q_rows.0[0].0.len()) */, "VOLE dimensions must correspond R1CS");
            let vith_size = s_adjusted.0.len() * s_adjusted.0[0].0.len();
            let mut q = Vec::with_capacity(vith_size);
            s_adjusted.0.iter_mut().for_each(|row| q.append(&mut row.0));
            let q = FVec::<T>(q);
            Self { delta, q, r1cs_with_metadata: r1cswm }
        }

        /// Verifies a (degree 2) Quicksilver proof, returning the public inputs and outputs if successfull. Otherwise, returns an error
        /// NOTE: According to the Quicksilver paper, `challenge` should be given after the values are determined.
        pub fn verify(&self, challenge: &T, proof: &ZKP<T>) -> Result<(), Error>{
            let r1cs = &self.r1cs_with_metadata.r1cs;
            let (q_a, q_b, q_c) = r1cs.vec_mul(&self.q);

            // Quicksilver protocol to transform VOLE into a new VOLE that makes multiplcation gates linear relations
            let new_q = &(&q_a * &q_b) - &q_c.scalar_mul(self.delta);
            let challenge_vec = get_challenge_vec(challenge, self.q.0.len());
            let success = proof.mul_proof.1 + proof.mul_proof.0 * self.delta == new_q.dot(&challenge_vec);
            match success {
                true => Ok(()),
                false => Err(anyhow!("Proof was not verified with success"))
            }
        }
        /// Assuming the VOLE was constructed properly, this verifies the opening of witness VOLE correlations
        pub fn verify_public(&self, pos: &PublicOpenings<T>) -> Result<(), Error> {
            if (!pos.public_inputs.len() == self.r1cs_with_metadata.public_inputs_indices.len()) && ((!pos.public_inputs.len() == self.r1cs_with_metadata.public_inputs_indices.len())) {
                bail!("Public values have the wrong input or output length(s)")
            }

            let mut indices = self.r1cs_with_metadata.public_inputs_indices.clone();
            let mut public = pos.public_inputs.clone();

            indices.extend(&self.r1cs_with_metadata.public_outputs_indices);
            public.extend(&pos.public_outputs);

            for (i, (u, v)) in indices.iter().zip(public.iter()) {
                // TODO: consider giving index of which input was invalid.  This could impact performance slightly as it would not be static but dynamic
                if !(*u * &self.delta + v == self.q.0[*i]) { bail!("Invaliding opening of a public input") }
            }
            Ok(())
        }
    } 
}

#[cfg(test)]
pub mod test {
    use ff::{Field, PrimeField};
    use lazy_static::lazy_static;
    use rand::rngs::ThreadRng;
    use crate::{FVec, Fr, zkp::quicksilver::Verifier};
    use super::{*, quicksilver::Prover};

    lazy_static! {
        pub static ref TEST_R1CS: FullR1CS<Fr> = {
            let a_rows = vec![
                FVec(vec![1, 1, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FVec(vec![2, 0, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            ];
            let b_rows = vec![
                FVec(vec![0, 2, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect())
            ];
            let c_rows = vec![
                FVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
                FVec(vec![0, 0, 0, 1].iter().map(|x|Fr::from_u128(*x)).collect())
            ];

            FullR1CS {
                a_rows: FMatrix(a_rows),
                b_rows: FMatrix(b_rows),
                c_rows: FMatrix(c_rows),
            }
        };
        pub static ref TEST_R1CS_WITH_METADA: R1CSWithMetadata<Fr> = R1CSWithMetadata { 
            r1cs: R1CS::Full(TEST_R1CS.clone()), 
            public_inputs_indices: vec![0,2], 
            public_outputs_indices: vec![3],
            unpadded_wtns_len: TEST_R1CS.a_rows.0.len(),
        };
    }
    #[test]
    fn circuit_satisfiability() {
        let witness = FVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect::<Vec<Fr>>());
        assert!(TEST_R1CS_WITH_METADA.r1cs.witness_check(&witness));
        assert!(!TEST_R1CS_WITH_METADA.r1cs.witness_check(&FVec(vec![Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ONE])));
    }

    #[test]
    pub fn circuit_satisfiability_proof() {
        let witness = FVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());

        // Prove it in ZK this time:
        let delta = Fr::random(&mut ThreadRng::default());
        let v = FVec::<Fr>::random(witness.0.len());
        let u = witness.clone();
        let q = &u.scalar_mul(delta) + &v;
        
        let prover = Prover {
            u,
            v: v.clone(),
            r1cs_with_metadata: TEST_R1CS_WITH_METADA.clone()
        };
        let challenge = &Fr::from_u128(123);
        let proof = prover.prove(challenge);

        let verifier = Verifier {
            q,
            delta,
            r1cs_with_metadata: TEST_R1CS_WITH_METADA.clone()
        };
        assert!(verifier.verify(challenge, &proof).is_ok());
        assert!(verifier.verify(&Fr::from_u128(69), &proof).is_err());
        // TODO: assert a bad witness fails (is this necessary tho bc ZK protocol will catch that lol)
    }

    // /// This is covered by practiaclly every single integration tests so commenting it instead of implementing it
    // #[test]
    // pub fn Tom_vith() {
    //     todo!()
    // }
}