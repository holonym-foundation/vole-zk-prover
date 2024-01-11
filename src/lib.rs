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
const NUM_VOLES: u32 = 1024;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
// Important this matches the endianness of MODULUS_AS_U128s
#[PrimeFieldReprEndianness = "big"]
pub struct Fr([u64; 4]);

/// Alias for types suitable for the prime field element
pub trait PF: PrimeField + Add + Sub + Mul + FromU8s + ToU8s {}
impl<T: PrimeField + Add + Sub + Mul + FromU8s + ToU8s> PF for T {}

/// A vector of field elements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FVec<T: PF>(pub Vec<T>);
// pub struct FVec<T><T: PrimeField>(pub Vec<T>);

// pub trait Field {
//     fn fast_secure_rand(&mut rng: ) {}
// }
// pub struct F<T: PrimeField>

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseVec<T: Mul + Add>(pub Vec<(usize, T)>);

pub trait FromU8s {
    fn from_u8s(u: &Vec<u8>) -> Self;
}
pub trait ToU8s {
    fn to_u8s(&self) -> Vec<u8>;
}
impl FromU8s for Fr {
    fn from_u8s(u: &Vec<u8>) -> Self {
        if u.len() != 32 { panic!("field element bust must be 32-byte") }
        Fr::from_repr(FrRepr(u[0..32].try_into().unwrap())).unwrap()
    }
}
impl ToU8s for Fr {
    fn to_u8s(&self) -> Vec<u8> {
        self.to_repr().0.try_into().unwrap()
    }
}
// impl FromBytes for Fr {
//     type Bytes = [u8; 32];
//     fn from_be_bytes(bytes: &Self::Bytes) -> Self {
//         Fr::from_repr(FrRepr(bytes.clone())).unwrap()
//     }
//     fn from_le_bytes(bytes: &Self::Bytes) -> Self {
//         let mut b = bytes.clone();
//         b.reverse();
//         Fr::from_repr(FrRepr(b)).unwrap()
//     }
// }
// impl ToBytes for Fr {
//     type Bytes = [u8; 32];
//     fn to_be_bytes(&self) -> Self::Bytes {
//         self.to_repr().0
//     }
//     fn to_le_bytes(&self) -> Self::Bytes {
//         let mut b = self.to_repr().0;
//         b.reverse();
//         b
//     }
// }


/// Pretty display
impl Display for FVec<Fr> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[ {} ]", self.0.iter().map(|fr|{
                let repr = fr.to_repr().0;
                format!("0x{:02x}..{:02x}", repr[0], repr[31])
            }).collect::<Vec<String>>().join(", ")
        )
    }
}
/// Pretty display
impl Display for FMatrix<Fr> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = self.0.iter().map(|fv|fv.to_string()).collect::<Vec<String>>();
        write!(f, "Matrix in row major order:\n[\n\t{}\n]", s.join("\n\t"))
    }
}

// TODO: clean up this ridiculous math trait derivation :p

impl<'a, 'b, T: PF> Mul<&'b FVec<T>> for &'a FVec<T> {
    type Output = FVec<T>;
    fn mul(self, rhs: &'b FVec<T>) -> FVec<T> {
        FVec::<T>(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).collect())
    }
}
impl<T: PF> Add for FVec<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a + *b).collect())
    }
}
impl<'a, 'b, T: PF> Add<&'b FVec<T>> for &'a FVec<T> {
    type Output = FVec<T>;
    fn add(self, rhs: &'b FVec<T>) -> FVec<T> {
        FVec::<T>(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a + *b).collect())
    }
}

impl<'a, 'b, T: PF> Sub<&'b FVec<T>> for &'a FVec<T> {
    type Output =  FVec<T>;
    fn sub(self, rhs: &'b FVec<T>) -> FVec<T> {
        FVec::<T>(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a - *b).collect())
    }
}
impl<'a, T: PF> SubAssign<FVec<T>> for &'a mut FVec<T>{
    fn sub_assign(&mut self, rhs: FVec<T>) {
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}

impl<'a, 'b, T: PF> Sub<&'b FVec<T>> for &'a mut FVec<T> {
    type Output =  FVec<T>;
    fn sub(self, rhs: &'b FVec<T>) -> FVec<T> {
        FVec::<T>(self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a - *b).collect())
    }
}

impl<'a, 'b, T: PF> SubAssign<&'b mut FVec<T>> for FVec<T> {
    fn sub_assign(&mut self, rhs: &'b mut FVec<T>) {
        // *self = FVec<T>(vec![Fr::ONE]);
        self.0.iter_mut().zip(rhs.0.iter()).for_each(|(a, b)| *a -= *b);
    }
}

impl<'a, T: PF> Neg for &'a FVec<T> {
    type Output = FVec<T>;
    fn neg(self) -> FVec<T> {
        FVec::<T>(self.0.iter().map(|a| -*a).collect())
    }
}

pub trait DotProduct<T: PF> {
    type Inner;
    fn dot(&self, rhs: &Self) -> Self::Inner;
    fn sparse_dot(&self, rhs: &SparseVec<T>) -> Self::Inner;
}
impl<T: PF> DotProduct<T> for FVec<T> {
    type Inner = T;
    fn dot(&self, rhs: &Self) -> Self::Inner {
        self.0.iter().zip(rhs.0.iter()).map(|(a, b)| *a * *b).sum::<T>()
    }
    // TODO: see whether this can be optimized
    fn sparse_dot(&self, rhs: &SparseVec<T>) -> Self::Inner {
        rhs.0.iter().fold(T::ZERO, |acc, (idx, val)|{
            acc + &(self.0[*idx] * val)
        })
    }
}

impl<T: PF> PartialEq for FVec<T> {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a, b)| a == b)
    }
}


impl<T: PF> FVec<T> {
    fn scalar_mul(&self, rhs: T) -> Self {
        Self(self.0.iter().map(|a| *a * rhs).collect())
    }
    /// Appends `len` zeroes
    pub fn zero_pad(&mut self, len: usize) {
        self.0.append(&mut vec![T::ZERO; len]);
    }
    pub fn random(len: usize) -> Self {
        let mut r = &mut ThreadRng::default();
        Self(
            (0..len).map(|_|T::random(&mut r)).collect()
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FMatrix<T: PF>(pub Vec<FVec<T>>);
impl<T: PF> FMatrix<T> {
    pub fn transpose(&self) -> Self {
        let outer_len = self.0.len();
        let inner_len  = self.0[0].0.len();
        let mut res = Vec::with_capacity(inner_len);
        for i in 0..inner_len {
            let mut new = Vec::with_capacity(outer_len);
            for j in 0..outer_len {
                new.push(self.0[j].0[i]);
            }
            res.push(FVec::<T>(new));
        }
        Self(res)
    }

    fn scalar_mul(&self, rhs: T) -> Self {
        Self(
            self.0.iter().map(|x|x.scalar_mul(rhs)).collect()
        )
    }

    pub fn dim(&self) -> (usize, usize) {
        (self.0[0].0.len(), self.0.len())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseFMatrix<T: PF>(pub Vec<SparseVec<T>>);


// impl MulAssign<ElementaryColumnOp> for FMatrix {
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

// impl MulAssign<ElementaryColumnOpComposition> for FMatrix {
//     fn mul_assign(&mut self, rhs: ElementaryColumnOpComposition) {
//         for op in rhs.0 {
//             self.mul_assign(op);
//         }
//     }
// }

impl<'a, 'b, T: PF> Add<&'b FMatrix<T>> for &'a FMatrix<T> {
    type Output = FMatrix<T>;
    fn add(self, rhs: &'b FMatrix<T>) -> FMatrix<T> {
        FMatrix(
            self.0.iter().zip(rhs.0.iter()).map(|(a, b)| a + b).collect()
        )
    }
}

impl<'a, 'b, T: PF> Sub<&'b FMatrix<T>> for &'a FMatrix<T> {
    type Output = FMatrix<T>;
    fn sub(self, rhs: &'b FMatrix<T>) -> FMatrix<T> {
        FMatrix::<T>(
            self.0.iter().zip(rhs.0.iter()).map(|(a, b)| a - b).collect()
        )
    }
}

impl<'a, 'b, T: PF> Mul<&'b FMatrix<T>> for &'a FVec<T> {
    type Output = FVec<T>;
    fn mul(self, rhs: &'b FMatrix<T>) -> FVec<T> {
        FVec::<T>(
            rhs.0.iter().map(|row_or_col| self.dot(row_or_col)).collect()
        )
    }
}

impl<'a, 'b, T: PF> Mul<&'b SparseFMatrix<T>> for &'a FVec<T> {
    type Output = FVec<T>;
    fn mul(self, rhs: &'b SparseFMatrix<T>) -> FVec<T> {
        FVec::<T>(
            rhs.0.iter().map(|row_or_col| self.sparse_dot(row_or_col)).collect()
        )
    }
}

impl<T: PF> PartialEq for FMatrix<T> {
    fn eq(&self, rhs: &Self) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a, b)| a == b)
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_transpose() {
        let x = FMatrix(vec![
            FVec(vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)]),
            FVec(vec![Fr::from(4u64), Fr::from(5u64), Fr::from(6u64)]),
            FVec(vec![Fr::from(7u64), Fr::from(8u64), Fr::from(9u64)]),
        ]);
        let x_t = FMatrix(vec![
            FVec(vec![Fr::from(1u64), Fr::from(4u64), Fr::from(7u64)]),
            FVec(vec![Fr::from(2u64), Fr::from(5u64), Fr::from(8u64)]),
            FVec(vec![Fr::from(3u64), Fr::from(6u64), Fr::from(9u64)]),
        ]);
        assert_eq!(x.transpose(), x_t);
    }
    
    // Could cover more edge cases
    #[test]
    fn test_sparse_vec() {
        let a = FVec(vec![Fr::ZERO, Fr::ZERO, Fr::ZERO, Fr::from_u128(69)]);
        let b = SparseVec(vec![(3, Fr::from_u128(100)), (2, Fr::from_u128(5))]);
        assert!(a.sparse_dot(&b) == Fr::from_u128(6900));
    }
}