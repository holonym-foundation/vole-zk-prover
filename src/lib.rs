pub mod circom;
pub mod zkp;
pub mod utils;
pub mod vith;
pub mod subspacevole;
pub mod challenges;
pub mod smallvole;
pub mod vecccom;
pub mod actors;
pub mod format;
pub mod codeparams;

use std::{ops::{Add, Mul, AddAssign, Neg, Sub, SubAssign, MulAssign}, process::Output, fmt::Display};
use wasm_bindgen::prelude::*;

use ff::Field;
use nalgebra::{ClosedMul, SMatrix, DMatrix};
use rand::rngs::ThreadRng;
use subspacevole::RAAACode;
use serde::{Serialize, Deserialize};
// use subspacevole::{ElementaryColumnOp, ElementaryColumnOpComposition};
// use num_traits::Zero;
#[macro_use]
extern crate ff;
use crate::ff::PrimeField;

/// Important that it is the block size of the linear code
const NUM_VOLES: usize = 1024;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
// Important this matches the endianness of MODULUS_AS_U128s
#[PrimeFieldReprEndianness = "big"]
pub struct Fr([u64; 4]);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrVec(pub Vec<Fr>);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseVec<T: Mul + Add>(pub Vec<(usize, T)>);

/// Pretty display
impl Display for FrVec {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[ {} ]", self.0.iter().map(|fr|{
                let repr = fr.to_repr().0;
                format!("0x{:02x}..{:02x}", repr[0], repr[31])
            }).collect::<Vec<String>>().join(", ")
        )
    }
}
/// Pretty display
impl Display for FrMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.0.iter().map(|fv|fv.to_string()).collect::<Vec<String>>();
        write!(f, "Matrix in row major order:\n[\n\t{}\n]", s.join("\n\t"))
    }
}
impl Mul for FrVec {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).collect())
    }
}

// TODO: clean up this ridiculous math trait derivation :p

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
    fn sparse_dot(&self, rhs: &SparseVec<Fr>) -> Self::Inner;
}
impl DotProduct for FrVec {
    type Inner = Fr;
    fn dot(&self, rhs: &Self) -> Self::Inner {
        self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).sum::<Fr>()
    }
    // TODO: see whether this can be optimized
    fn sparse_dot(&self, rhs: &SparseVec<Fr>) -> Self::Inner {
        rhs.0.iter().fold(Fr::ZERO, |acc, (idx, val)|{
            acc + &(self.0[*idx] * val)
        })
    }
}

pub trait ScalarMul<T> {
    fn scalar_mul(&self, rhs: T) -> Self;
}

impl<'b> ScalarMul<&'b Fr> for FrMatrix {
    fn scalar_mul(&self, rhs: &'b Fr) -> Self {
        Self(
            self.0.iter().map(|x|x.scalar_mul(rhs)).collect()
        )
    }
}

impl<'b> ScalarMul<&'b Fr> for FrVec {
    fn scalar_mul(&self, rhs: &'b Fr) -> Self {
        Self(self.0.iter().map(|a| *a * *rhs).collect())
    }
}

impl PartialEq for FrVec {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a, b)| a == b)
    }
}


impl FrVec {
    /// Appends `len` zeroes
    pub fn zero_pad(&mut self, len: usize) {
        self.0.append(&mut vec![Fr::ZERO; len]);
    }
    pub fn random(len: usize) -> Self {
        let mut r = &mut ThreadRng::default();
        Self(
            (0..len).map(|_|Fr::random(&mut r)).collect()
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseFrMatrix(pub Vec<SparseVec<Fr>>);


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

impl<'a, 'b> Add<&'b FrMatrix> for &'a FrMatrix {
    type Output = FrMatrix;
    fn add(self, rhs: &'b FrMatrix) -> FrMatrix {
        FrMatrix(
            self.0.iter().zip(rhs.0.iter()).map(|(a, b)| a + b).collect()
        )
    }
}

impl<'a, 'b> Sub<&'b FrMatrix> for &'a FrMatrix {
    type Output = FrMatrix;
    fn sub(self, rhs: &'b FrMatrix) -> FrMatrix {
        FrMatrix(
            self.0.iter().zip(rhs.0.iter()).map(|(a, b)| a - b).collect()
        )
    }
}

impl<'a, 'b> Mul<&'b FrMatrix> for &'a FrVec {
    type Output = FrVec;
    fn mul(self, rhs: &'b FrMatrix) -> FrVec {
        FrVec(
            rhs.0.iter().map(|row_or_col| self.dot(row_or_col)).collect()
        )
    }
}

impl<'a, 'b> Mul<&'b SparseFrMatrix> for &'a FrVec {
    type Output = FrVec;
    fn mul(self, rhs: &'b SparseFrMatrix) -> FrVec {
        FrVec(
            rhs.0.iter().map(|row_or_col| self.sparse_dot(row_or_col)).collect()
        )
    }
}

// // --
// /// This is weird since the rhs should be the matrix. The reason it is done this way too is so 
// /// that the R1CS or matrices in general can be an enum of either Sparse or Full, both of which 
// /// implement Mul<FrVec>
// impl<'a, 'b> Mul<&'b FrVec> for &'a FrMatrix {
//     type Output = FrVec;
//     fn mul(self, lhs: &'b FrVec) -> FrVec {
//         FrVec(
//             self.0.iter().map(|row_or_col| lhs.dot(row_or_col)).collect()
//         )
//     }
// }

// /// This is weird since the rhs should be the matrix. The reason it is done this way too is so 
// /// that the R1CS or matrices in general can be an enum of either Sparse or Full, both of which 
// /// implement Mul<FrVec>
// impl<'a, 'b> Mul<&'b FrVec> for &'a Vec<SparseVec<Fr>> {
//     type Output = FrVec;
//     fn mul(self, lhs: &'b FrVec) -> FrVec {
//         FrVec(
//             self.iter().map(|row_or_col| lhs.sparse_dot(row_or_col)).collect()
//         )
//     }
// }


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
    
    // Could cover more edge cases
    #[test]
    fn test_sparse_vec() {
        let a = FrVec(vec![Fr::ZERO, Fr::ZERO, Fr::ZERO, Fr::from_u128(69)]);
        let b = SparseVec(vec![(3, Fr::from_u128(100)), (2, Fr::from_u128(5))]);
        assert!(a.sparse_dot(&b) == Fr::from_u128(6900));
    }
}