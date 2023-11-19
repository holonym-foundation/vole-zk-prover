pub mod zkp;
pub mod vith;
pub mod subspacevole;
pub mod smallvole;
pub mod vecccom;

use std::{ops::{Add, Mul, AddAssign, Neg, Sub, SubAssign, MulAssign}, process::Output};

use ff::Field;
use nalgebra::{ClosedMul, SMatrix, DMatrix};
use subspacevole::RAAACode;
// use subspacevole::{ElementaryColumnOp, ElementaryColumnOpComposition};
// use num_traits::Zero;
#[macro_use]
extern crate ff;
use crate::ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
// Important this matches the endianness of MODULUS_AS_U128s
#[PrimeFieldReprEndianness = "big"]
pub struct Fr([u64; 4]);

#[derive(Debug, Clone)]
pub struct FrVec(pub Vec<Fr>);

impl Mul for FrVec {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).collect())
    }
}


impl<'a, 'b> Mul<&'b FrVec> for &'a FrVec {
    type Output = FrVec;
    fn mul(self, rhs: &'b FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).collect())
    }
}
impl Add for FrVec {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a + *b).collect())
    }
}
impl<'a, 'b> Add<&'b FrVec> for &'a FrVec {
    type Output = FrVec;
    fn add(self, rhs: &'b FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a + *b).collect())
    }
}
impl<'a, 'b> Add<&'b FrVec> for &'a mut FrVec {
    type Output = FrVec;
    fn add(self, rhs: &'b FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a + *b).collect())
    }
}

impl<'a, 'b> Sub<&'b FrVec> for &'a FrVec {
    type Output =  FrVec;
    fn sub(self, rhs: &'b FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a - *b).collect())
    }
}
impl<'a> SubAssign<FrVec> for &'a mut FrVec{
    fn sub_assign(&mut self, rhs: FrVec) {
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}

impl<'a, 'b> Sub<&'b FrVec> for &'a mut FrVec {
    type Output =  FrVec;
    fn sub(self, rhs: &'b FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a - *b).collect())
    }
}

impl<'a, 'b> Sub<&'b mut FrVec> for &'a mut FrVec {
    type Output =  FrVec;
    fn sub(self, rhs: &'b mut FrVec) -> FrVec {
        FrVec(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a - *b).collect())
    }
}

impl<'a, 'b> SubAssign<&'b FrVec> for FrVec {
    fn sub_assign(&mut self, rhs: &'b FrVec) {
        // *self = FrVec(vec![Fr::ONE]);
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}
impl<'a, 'b> SubAssign<&'b mut FrVec> for FrVec {
    fn sub_assign(&mut self, rhs: &'b mut FrVec) {
        // *self = FrVec(vec![Fr::ONE]);
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}

impl SubAssign for FrVec {
    fn sub_assign(&mut self, rhs: Self) {
        // *self = FrVec(vec![Fr::ONE]);
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}

impl<'a> Neg for &'a FrVec {
    type Output = FrVec;
    fn neg(self) -> FrVec {
        FrVec(self.0.iter().map(|a| -*a).collect())
    }
}

pub trait DotProduct {
    type Inner;
    fn dot(&self, rhs: &Self) -> Self::Inner;
}

// impl DotProduct for Vec<Fr> {
//     type Inner = Fr;
//     fn dot(&self, rhs: &Self) -> Self::Inner {
//         self.iter().zip(rhs.iter()).map(|(a, b)| *a * *b).sum::<Fr>()
//     }
// }
impl PartialEq for FrVec {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a, b)| a == b)
    }
}
impl DotProduct for FrVec {
    type Inner = Fr;
    fn dot(&self, rhs: &Self) -> Self::Inner {
        self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).sum::<Fr>()
    }
}

impl FrVec {
    fn scalar_mul(&self, rhs: &Fr) -> Self {
        Self(self.0.iter().map(|a| *a * *rhs).collect())
    }
}

#[derive(Debug, Clone)]
pub struct FrMatrix(pub Vec<FrVec>);
impl FrMatrix {
    pub fn transpose(&self) -> Self {
        let outer_len = self.0.len();
        let inner_len  = self.0[0].0.len();
        let mut res = Vec::with_capacity(inner_len);
        for i in 0..inner_len {
            let mut new = Vec::with_capacity(outer_len);
            for j in 0..outer_len {
                new.push(self.0[j].0[i]);
            }
            res.push(FrVec(new));
        }
        Self(res)
    }
    pub fn dim(&self) -> (usize, usize) {
        (self.0[0].0.len(), self.0.len())
    }
}

// impl MulAssign<ElementaryColumnOp> for FrMatrix {
//     fn mul_assign(&mut self, rhs: ElementaryColumnOp) {
//         match rhs {
//             ElementaryColumnOp::Swap(i, j) => {
//                 self.0.swap(i, j);
//             },
//             ElementaryColumnOp::Scale(s, i) => {
//                 self.0[i] = self.0[i].scalar_mul(&s);
//             },
//             ElementaryColumnOp::AddMultiple(s, i, j) => {
//                 self.0[j] = &self.0[j] + &(self.0[i].scalar_mul(&s));
//             },
//         }
//     }
// }

// impl MulAssign<ElementaryColumnOpComposition> for FrMatrix {
//     fn mul_assign(&mut self, rhs: ElementaryColumnOpComposition) {
//         for op in rhs.0 {
//             self.mul_assign(op);
//         }
//     }
// }

impl<'a, 'b> Mul<&'b FrVec> for &'a FrMatrix {
    type Output = FrVec;
    fn mul(self, rhs: &'b FrVec) -> FrVec {
        FrVec(
            self.0.iter().map(|vec| vec.dot(rhs)).collect()
        )
    }
}

impl PartialEq for FrMatrix {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a, b)| a == b)
    }
}

impl polynomen::One for Fr {
    fn one() -> Self {
        Fr::ONE
    }
    fn is_one(&self) -> bool {
       self.eq(&Fr::ONE)
    }
}


impl polynomen::Zero for Fr {
    fn zero() -> Self {
        Fr::ZERO
    }
    fn is_zero(&self) -> bool {
       self.eq(&Fr::ZERO)
    }
}

impl num_traits::One for Fr {
    fn one() -> Self {
        Fr::ONE
    }
    fn is_one(&self) -> bool {
       self.eq(&Fr::ONE)
    }
}
impl num_traits::Zero for Fr {
    fn zero() -> Self {
        Fr::ZERO
    }
    fn is_zero(&self) -> bool {
       self.eq(&Fr::ZERO)
    }
}

// impl ClosedMul for Fr {}
// impl std::ops::Mul for Fr {
//     type Output = Self;
//     fn mul(self, rhs: Self) -> Self {
//         let mut out = self;
//         out.mul_assign(rhs);
//         out
//     }
// }
// impl MulAssign for Fr {
//     fn mul_assign(&mut self, rhs: Self) {
//         *self = *self * rhs;
//     }
// }
// // Crates such as nalgebra do not allow 
// fn matmul_fr() {

// }

pub struct ProverCommitment {
    /// Hash of every pair of seed's respective hashes for the seeds used to create the VOLEs. We are just using two seeds per VOLE!
    /// Can/should be used for Fiat-Shamir of subspace VOLE consistency check
    pub seed_comm: [u8; 32],
    /// Witness split into vectors of the same length as the code's dimension
    pub witness_comm: FrMatrix,
    /// subsapce VOLE consistency check of U and V's check values, respectively
    pub consistency_check: (FrVec, FrVec)
}
pub struct ProverOpening {
    /// Openings of one seed per pair
    pub seed_opens: Vec<[u8; 32]>,
    /// Proofs that the openings were done correctly
    pub seed_proofs: Vec<[u8; 32]>,
    /// S matrix from the final VOLE in the head
    pub vith_s: FrMatrix,
    /// R1 and U1 values for the final gate in the ZKP
    pub final_gate: (Fr, Fr)
}

pub struct Prover {
    pub code: RAAACode,
    pub vole_length: usize,
    pub num_voles: usize,
    pub witness: FrVec
}
impl Prover {
    // fn committ(&self) -> ProverCommitment {

    // }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_transpose() {
        let x = FrMatrix(vec![
            FrVec(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]),
            FrVec(vec![Fr::from(4u64), Fr::from(5u64), Fr::from(6u64)]),
            FrVec(vec![Fr::from(7u64), Fr::from(8u64), Fr::from(9u64)]),
        ]);
        let x_t = FrMatrix(vec![
            FrVec(vec![Fr::from(1u64), Fr::from(4u64), Fr::from(7u64)]),
            FrVec(vec![Fr::from(2u64), Fr::from(5u64), Fr::from(8u64)]),
            FrVec(vec![Fr::from(3u64), Fr::from(6u64), Fr::from(9u64)]),
        ]);
        assert_eq!(x.transpose(), x_t);
    }
}