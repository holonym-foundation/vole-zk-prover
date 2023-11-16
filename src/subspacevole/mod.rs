use std::ops::{Add, Mul, Neg, MulAssign};
use std::time::Instant;
use std::usize;

use ff::PrimeField;
use nalgebra::{self, SMatrix, };
use polynomen::Poly;
use crate::{Fr, FrRepr, DotProduct, FrVec, FrMatrix};
use crate::ff::Field;


pub fn sender_correction<const N: usize, const K: usize, const L: usize, const NminusK: usize>(U: &SMatrix<Fr, L, N>) -> SMatrix<Fr, N, NminusK> {
    assert!(NminusK == N - K, "NminusK must be N - K");
    let start = Instant::now();
    let tc_inverse = ReedSolomonCode::construct_tc_inverse::<N>();
    println!("inversion time: {:?}", start.elapsed());

    let start = Instant::now();
    let U_tcinv = U * tc_inverse;
    println!("multiplication time: {:?}", start.elapsed());
    U_tcinv.fixed_slice::<N, NminusK>(K, 0).into()
}
// pub trait SubspaceVOLESender {
//     /// Returns the U matrix
//     fn u(&self);
//     /// Returns the V matrix
//     fn v(&self);
//     /// Returns the correction vector so that the Receiver can correct their ∆, Q
//     fn correction(&self);
// }

// // Very inefficient (128+32 VOLES to get 32/2 elements of the witness), just for experiment with while i'm writing the library
// pub struct ExampleSender {
//     U: SMatrix<Fr, 128, 32>,
//     V: SMatrix<Fr, 128, 32>,
// }
// impl SubspaceVOLESender for ExampleSender {
//     fn u(&self) {
//         self.U
//     }
//     fn v(&self) {
//         self.V
//     }
//     fn correction() {
//         unimplemented!()
//     }
// }

// pub trait LinearCode<T, const N: usize, const K: usize>  {
//     fn generator() -> SMatrix<T, K, N>;
//     /// Tc matrix from the paper, extra rows appended to the generator matrix 
//     fn extended_generator() -> SMatrix<T, N, N>;
// }

// struct ReedSolomonCode;
// /// Creates a Vandermode Matrix for given dimensions
// impl<T, constN, constK> LinearCode<T, constN, constK> for ReedSolomonCode {
//     fn generator() -> SMatrix<T, N, K> {
//         for i in 0..K {
//             for j in 0..N {
//                 let mut tmp = Fr::one();
//                 for k in 0..K {
//                     if k != i {
//                         tmp *= Fr::from(j) - Fr::from(k);
//                     }
//                 }
//                 tmp = tmp.invert().unwrap();
//                 generator[i][j] = tmp;
//             }
//         }
//     }
// }

/// Represents an an elementary operation on the columns of a matrix
pub enum ElementaryColumnOp {
    /// Scales the column at the given index by the given scalar
    Scale(Fr, usize),
    /// Scales the column at the first index and adds it to the column at the second index
    AddMultiple(Fr, usize, usize),
    /// Swaps two columns at the given indices
    Swap(usize, usize),
}
impl ElementaryColumnOp {
    fn invert(&self) -> Option<Self> {
        match self {
            ElementaryColumnOp::Scale(scalar, idx) => {
                let inv = scalar.invert();
                if inv.is_none().unwrap_u8() == 1 { return None }
                Some(
                    ElementaryColumnOp::Scale(scalar.invert().unwrap(), *idx)
                )
            },
            ElementaryColumnOp::AddMultiple(scalar, idx1, idx2) => {
                Some(
                    ElementaryColumnOp::AddMultiple(scalar.neg(), *idx1, *idx2)
                )
            },
            ElementaryColumnOp::Swap(idx1, idx2) => {
                Some(
                    ElementaryColumnOp::Swap(*idx1, *idx2)
                )
            },
        }
    }
}
pub struct ElementaryColumnOpComposition(pub Vec<ElementaryColumnOp>);
impl ElementaryColumnOpComposition {
    pub fn inverse(&self) -> Option<Self> {
        let res: Option<Vec<_>> = self.0.iter().rev().map(|op| { 
            op.invert()
        }).collect();

        match res {
            None => None,
            Some(res) => Some(Self(res))
        }
    }
}

// Converts a matrix to systematic form 
pub fn convert_to_systematic() {
    todo!()
    // Great starting point (but need to multiply the row before subtracting, then also need to do the same process from the bottom up)
    // Also need to return the elementary column operations performed
    //     assert!(N >= K, "N must >= K");
    //     let mut g = ReedSolomonCode::construct_generator_quickly::<N, K>();
    //     // Convert to CEF
    //     // for (idx, mut col) in g.0[1..K].iter_mut().enumerate() {
    //     for i in 1..K {
    //         // let i = idx + 1;
    //         // This clone is likely slightly expensive but I don't see an easy way around it in Rust because it won't allow mutable borrowing twice
    //         for j in 0..i {
    //             let prev_row = &g.0[j].clone();
    //             g.0[i] -= prev_row;
    //             g.0[i].0[i].invert().unwrap(); // Unwrap is fine because it will never be zero for a Vandermode matrix, and thus invert() should always be Some
    //         }
    //     }
}
/// Gets idx'th Lagrange basis for a polynomial of degree deg
/// idx is 1-indexed
/// TODO: make sure Poly keeps 0 coefficients rather than return a shorter coefficient vector. This would cause a panic
/// although it does not seem a crash would leak any info since these matrices are public...
/// TODO 0-idx instead of 1-idx because i don't see a need to 1-idx
pub fn lagrange_basis_coeffs(idx: u64, deg: u64) -> Vec<Fr> {
        let roots = (1..deg+1)
            .filter(|x| *x != idx)
            .map(|x| Fr::from(x));

        let poly = Poly::new_from_roots_iter(roots.clone());
        // A lazy but expensive way to calculate the scaling factor lol (TODO: implement a better way, if it affects performance)
        let scaling_factor = poly.eval(&Fr::from(idx)).invert().unwrap();
        // Here is the 'right' way but for some reason it is not working properly:
        // let scaling_factor = roots.reduce(|prev: Fr, next| prev * (fr_idx - next)).unwrap().invert().unwrap();

        let coeffs: Vec<Fr> = poly.coeffs().iter().map(|x| *x * scaling_factor).collect();
        debug_assert!((1..deg+1).all(|x|{
            if x == idx {
                Poly::new_from_coeffs(&coeffs).eval(&Fr::from(x)) == Fr::ONE
            } else {
                Poly::new_from_coeffs(&coeffs).eval(&Fr::from(x)) == Fr::ZERO
            }
        }));
    coeffs   
}

/// 0 indexed
pub fn lagrange_basis(idx: u64, deg: u64) -> Poly<Fr> {
    let roots = (0..deg)
        .filter(|x| *x != idx)
        .map(|x| Fr::from(x));

    let poly = Poly::new_from_roots_iter(roots.clone());
    // A lazy but expensive way to calculate the scaling factor lol (TODO: implement a better way, if it affects performance)
    let scaling_factor = poly.eval(&Fr::from(idx)).invert().unwrap();
    // Here is the 'right' way but for some reason it is not working properly:
    // let scaling_factor = roots.reduce(|prev: Fr, next| prev * (fr_idx - next)).unwrap().invert().unwrap();

    let coeffs: Vec<Fr> = poly.coeffs().iter().map(|x| *x * scaling_factor).collect();
    debug_assert!((0..deg).all(|x|{
        if x == idx {
            Poly::new_from_coeffs(&coeffs).eval(&Fr::from(x)) == Fr::ONE
        } else {
            Poly::new_from_coeffs(&coeffs).eval(&Fr::from(x)) == Fr::ZERO
        }
    })); 
    Poly::new_from_coeffs(&coeffs)
}
pub struct ReedSolomonCode;
impl ReedSolomonCode {
    // Cosntructs a Vandermode Matrix for N and K
    pub fn construct_generator<const N: usize, const K: usize>() -> SMatrix<Fr, K, N> {
        SMatrix::<Fr, K, N>::from_fn(|row, col| {
            let base = Fr::from(row as u64 + 1);
            base.pow([col as u64])
        })
    }

    pub fn construct_generator_quickly<const N: usize, const K: usize>() -> FrMatrix {
        assert!(N >= K, "N must be greater than or equal to K");
        let mut out = Vec::with_capacity(K);
        let zeroth_col = FrVec(vec![Fr::ONE; N]);
        let first_col = FrVec((1..N+1).map(|x| Fr::from(x as u64)).collect());
        out.push(zeroth_col);
        out.push(first_col.clone());
            for i in 2..K {
               out.push(&out[i-1] * &first_col);
            }
        // FrMatrix(out)
        FrMatrix(out)
    }

    pub fn construct_systematic_generator<const N: usize, const K: usize>() -> FrMatrix {
        let mut res: Vec<FrVec> = Vec::with_capacity(K);
        for i in 0..K {
            let p = lagrange_basis(i as u64, K as u64);
            let evals = (0..N).map(|x| p.eval(&Fr::from(x as u64))).collect::<Vec<_>>();
            res.push(FrVec(evals));
        }
        FrMatrix(res)
        // println!("eval at [0,1,2,3,N-1]: {:?}", (0..N).map(|x| p.eval(&Fr::from(x as u64))).collect::<Vec<_>>());
    }

    /// Creates an invertible NxN matrix starting with the NxK generator matrix
    pub fn create_tc_for_systematic_genetor<const N: usize, const K: usize>() -> FrMatrix {
        let mut g = ReedSolomonCode::construct_systematic_generator::<N, K>();
        println!("Dimensions: {:?}", g.dim());
        for i in 0..(N-K) {
            let mut new_col = g.0[i].clone();
            new_col.0[0] += Fr::ONE;
        }
        todo!()
    }
    // // /// Creates the generator and converts it to reduced column echelon form. I think converting CCEF will be quicker constructing a bunch of Lagrange polynomials
    // pub fn construct_systematic_generator_quickly<const N: usize, const K: usize>() -> FrMatrix {
    //     assert!(N >= K, "N must >= K");
    //     let mut g = ReedSolomonCode::construct_generator_quickly::<N, K>();
    //     // Convert to CEF
    //     // for (idx, mut col) in g.0[1..K].iter_mut().enumerate() {
    //     for i in 1..K {
    //         // let i = idx + 1;
    //         // This clone is likely slightly expensive but I don't see an easy way around it in Rust because it won't allow mutable borrowing twice
    //         for j in 0..i {
    //             let prev_row = &g.0[j].clone();
    //             g.0[i] -= prev_row;
    //             g.0[i].0[i].invert().unwrap(); // Unwrap is fine because it will never be zero for a Vandermode matrix, and thus invert() should always be Some
    //         }
    //     }
    //     println!("CEF {:?}", g);
    //     // Convert to RCEF
    //     todo!()
    //     // FrMatrix(Vec<FrVec::)
    // }
    // Exploits the fact that an NxN Vandermonde matrix also satisfies the requirements for Tc matrix
    pub fn construct_tc<const N: usize>() -> SMatrix<Fr, N, N> {
        ReedSolomonCode::construct_generator::<N, N>()
    }
    // Exploits the fact that the columns of the inverse Vandermodnde matrix are the coefficients of the Lagrange bases
    pub fn construct_tc_inverse<const N: usize>() -> SMatrix<Fr, N, N> {
        let mut out = Vec::with_capacity(N * N);
        for i in 1..N+1 {
            let coeffs = lagrange_basis_coeffs(i as u64, N as u64);
            out.append(&mut coeffs.clone());
        }
        SMatrix::<Fr, N, N>::from_column_slice(&out) // Switch to columns for transpose
    }
}

#[cfg(test)]
mod test {
    use std::{ops::Mul, time::Instant};

    use ff::{Field, PrimeField};
    use nalgebra::{Matrix2x4, Matrix4x2};
    use rand::rngs::ThreadRng;

    use crate::FrRepr;

    use super::*;
    #[test]
    fn rs_is_vandermonde() {
        let g = ReedSolomonCode::construct_generator::<5,3>();
        assert_eq!(
            g,
            SMatrix::<Fr, 3, 5>::from_row_slice(&[
                Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
                Fr::from(1u64), Fr::from(2u64), Fr::from(4u64), Fr::from(8u64), Fr::from(16u64),
                Fr::from(1u64), Fr::from(3u64), Fr::from(9u64), Fr::from(27u64), Fr::from(81u64), 
            ])
        );
    }
    #[test]
    fn tc_times_inv_eq_identity() {
        let tc = ReedSolomonCode::construct_tc::<5>();
        let tcinv = ReedSolomonCode::construct_tc_inverse::<5>();
        let should_be_identity = tc * tcinv;
        println!("should_be_identity: {:?}", should_be_identity);
        let identity = SMatrix::<Fr, 5, 5>::identity();
        assert_eq!(
            should_be_identity,
            identity
        );
    }

    #[test]
    fn test_construct_generator_quickly() {
        let slow = ReedSolomonCode::construct_generator::<4, 3>();
        let fast = ReedSolomonCode::construct_generator_quickly::<4, 3>();
        println!("slow {:?}", slow);
        println!("fast {:?}", fast);
        // assert_eq!(slow., fast);
    }
    #[test]
    fn test_elementary_op_inverses() {
        todo!()
    }
    #[test]
    fn test_elementary_op_composite_inverses() {
        todo!()
    }
    #[test]
    fn test_systematic_code_creation() {
        todo!("check converting a known matrix to its known systematic form");
        todo!("alternatively, check systematic form spans the same subspace")
    }
    
    #[test]
    fn asdf() {
        // let start = Instant::now();
        // let tc = ReedSolomonCode::construct_tc::<64>();
        // println!("elapsed time: {:?}", start.elapsed());
        // let start = Instant::now();
        // let tcinv = ReedSolomonCode::construct_tc_inverse::<64>();
        // println!("elapsed time: {:?}", start.elapsed());
        // let U = SMatrix::<Fr, 45, 30>::from_fn(|row, col| {
        //     Fr::random(&mut ThreadRng::default())
        // });
        // sender_correction::<30, 20, 45, 10>(&U);

        ReedSolomonCode::create_tc_for_systematic_genetor::<8,4>();
    }
}