use std::ops::{Add, Mul};

use nalgebra::{self, SMatrix, };
use polynomen::Poly;
use crate::Fr;
use crate::ff::Field;

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

/// Gets idx'th Lagrange basis for a polynomial of degree deg
/// idx is 1-indexed
/// TODO: make sure Poly keeps 0 coefficients rather than return a shorter coefficient vector. This would cause a panic
/// although it does not seem a crash would leak any info since these matrices are public...
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

struct ReedSolomonCode;
impl ReedSolomonCode {
    // Cosntructs a Vandermode Matrix for N and K
    pub fn construct_generator<const N: usize, const K: usize>() -> SMatrix<Fr, K, N> {
        SMatrix::<Fr, K, N>::from_fn(|row, col| {
            let base = Fr::from(row as u64 + 1);
            base.pow([col as u64])
        })
    }
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
        SMatrix::<Fr, N, N>::from_row_slice(&out) // Switch to columns for transpose
    }
}

#[cfg(test)]
mod test {
    use std::ops::Mul;

    use ff::Field;
    use nalgebra::{Matrix2x4, Matrix4x2};

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

        todo!("implement")
    }
}