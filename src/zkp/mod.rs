use crate::{FrMatrix, Fr, FrVec};

pub struct R1CS {
    a: FrMatrix,
    b: FrMatrix,
    c: FrMatrix,
}
impl R1CS {
    /// Checks whether it is satisfiable by the witness
    fn witness_check(&self, witness: &FrVec) -> bool {
        (witness * &self.a)
        *
        (witness * &self.b)
        ==
        (witness * &self.c)
    }
}

pub mod quicksilver {
    use crate::{FrVec, Fr, FrMatrix};

    use super::R1CS;
    pub struct ZKP {
        /// Quicksilver multiplication proof of one field element
        pub mul_proof: Fr,
        /// Values for the final gate in the ZKP
        pub last_gate_opening: (Fr, Fr)
    }
    pub struct Prover {
        pub u: FrVec,
        pub v: FrVec,
        pub r1cs: R1CS
    }
    impl Prover {
        /// Creates a prover from VitH U1 and R matrices of equal dimension with 2l+2 rows where the witness is split into l chunks of length vole_length
        /// Takes ownership and mutates most of its inputs to something useless
        pub fn from_vith(mut u1_rows: FrMatrix, mut r_rows: FrMatrix, mut witness_rows: FrMatrix, r1cs: R1CS) -> Prover {
            assert!((u1_rows.0.len() == r_rows.0.len()) && (u1_rows.0[0].0.len() == r_rows.0[0].0.len()), "u and v must be same dimension");
            assert!(witness_rows.0.len() == u1_rows.0.len() - 1, "witness must have one fewer row than u1");
            assert!(witness_rows.0[0].0.len() == u1_rows.0[0].0.len() - 1, "witness must have same number of columns as u1");
            assert!((r1cs.a.0.len() == u1_rows.0.len()) && (r1cs.a.0[0].0.len() == u1_rows.0[0].0.len()), "VOLE dimensions must match R1CS dimensions");

            let vith_size = u1_rows.0.len() * u1_rows.0[0].0.len();
            let mut u = Vec::with_capacity(vith_size);
            let mut v = Vec::with_capacity(vith_size);
            witness_rows.0.iter_mut().map(|row| u.append(&mut row.0));
            u.append(&mut u1_rows.0.last().unwrap().0.clone());
            r_rows.0.iter_mut().map(|row| v.append(&mut row.0));

            Self { u: FrVec(u), v: FrVec(v), r1cs }
        }
        pub fn prove() -> ZKP {
            todo!()
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
            assert!((r1cs.a.0.len() == q_rows.0.len()) && (r1cs.a.0[0].0.len() == q_rows.0[0].0.len()), "VOLE dimensions must match R1CS dimensions");
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
    use rand::rngs::ThreadRng;
    use crate::{FrVec, Fr, ScalarMul};
    use super::*;

    #[test]
    fn circuit_satisfiability() {
        let witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());
        let a = vec![
            FrVec(vec![1, 1, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![2, 0, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            
        ];
        let b = vec![
            FrVec(vec![0, 2, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect())
        ];
        let c = vec![
            FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![0, 0, 0, 1].iter().map(|x|Fr::from_u128(*x)).collect())
        ];

        let r = R1CS {
            a: crate::FrMatrix(a),
            b: crate::FrMatrix(b),
            c: crate::FrMatrix(c),
        };

        assert!(r.witness_check(&witness));
        assert!(!r.witness_check(&FrVec(vec![Fr::ONE, Fr::ZERO, Fr::ZERO, Fr::ONE])));
    }

    #[test]
    pub fn circuit_satisfiability_proof() {
        let witness = FrVec(vec![5, 2, 28, 280].iter().map(|x|Fr::from_u128(*x)).collect());
        let a = vec![
            FrVec(vec![1, 1, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![2, 0, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            
        ];
        let b = vec![
            FrVec(vec![0, 2, 0, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect())
        ];
        let c = vec![
            FrVec(vec![0, 0, 1, 0].iter().map(|x|Fr::from_u128(*x)).collect()),
            FrVec(vec![0, 0, 0, 1].iter().map(|x|Fr::from_u128(*x)).collect())
        ];

        let r = R1CS {
            a: crate::FrMatrix(a),
            b: crate::FrMatrix(b),
            c: crate::FrMatrix(c),
        };

        // Prove it in ZK this time:
        let delta = Fr::random(&mut ThreadRng::default());
        let v = FrVec::random(witness.0.len());
        let u = witness;
        let q = u.scalar_mul(&delta) + v;
        todo!()

    }

    #[test]
    pub fn from_vith() {

    }
}