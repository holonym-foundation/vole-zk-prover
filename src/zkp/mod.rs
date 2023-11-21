use crate::{FrMatrix, Fr, FrVec};

pub struct R1CS {
    a: FrMatrix,
    b: FrMatrix,
    c: FrMatrix,
}
impl R1CS {
    /// Checks whether it is satisfiable by the witness
    fn witness_check(&self, witness: &FrVec) -> bool {
        (&self.a * witness)
        *
        (&self.b * witness)
        ==
        (&self.c * witness)
    }
}
#[cfg(test)]
mod test {
    use ff::{Field, PrimeField};

    use crate::{FrVec, Fr};

    use super::R1CS;

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
}

pub mod quicksilver {
    use crate::{FrVec, Fr};

    use super::R1CS;

    struct Prover {
        pub u: FrVec,
        pub v: FrVec,
    }
    struct Verifier {
        pub delta: Fr,
        pub q: FrVec
    }
    impl Verifier {
        pub fn addition_gates(&self, r1cs: &R1CS) {
            
        }
    }
}
