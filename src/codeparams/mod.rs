//! Calculates a lower bound for security for a Repeat Multiple Accumulate linear code
//! Note this underestimates the security. It simply provides
//! A lower bound on the probaility of getting a code with minimum distance d.
//! It is based off of 
//! "The Serial Concatenation of Rate-1 Codes Through Uniform Random Interleavers" by Henry D. Pfister and Paul H. Siegel [year]
//! "Coding Theorems for Turbo-Like Codes" by Divsalar [year]
//! And a proof the security should be no lower when used prime-field inputs rather than binary-field inputs which may be found in the Holonym V2 whitepaper

use itertools::Itertools;
use nalgebra::{SMatrix, Const, Vector, SVector};
use num_bigint::{BigInt, BigUint};
use num_integer::{binomial, Integer};
use num_traits::{FromPrimitive, ToPrimitive};

/// Batch n-choose-k, exploting patterns in the binomial coefficient for efficient batch calcualtion in a particular range
/// This only does it for square nxn matrices but coudlbe easily modified to suport abritrary matrix dimensions
fn n_choose_k_square_matrix(n: usize) -> Vec<Vec<BigUint>>{
    let mut init_row = vec![BigUint::from_u8(0).unwrap(); n+1];
    init_row[0] = BigUint::from_u8(1).unwrap();
    let mut output = vec![init_row; n+1];

    for i in 1..n+1 {
        for j in 1..n+1 {
            output[i][j] = &output[i-1][j-1] + &output[i-1][j];
        }
    }

    output
}
/// Calculates the Input Output Weight Enumeration for an accumulate code with input hamming weight w, output hamming weight h, and block size n
pub fn calc_iowe_entry(input_hamming: usize, output_hamming: usize, block_size: usize) -> u128 {
    if (output_hamming != 0) && (input_hamming == 0) { return 0 }
    let w = BigInt::from_usize(input_hamming).unwrap();
    let h = BigInt::from_usize(output_hamming).unwrap();
    let n = BigInt::from_usize(block_size).unwrap();

    let two = BigInt::from_usize(2).unwrap();

    let lhs = binomial(n - &h, w.div_floor(&two));
    let rhs = binomial(h - 1, w.div_ceil(&two) - 1);

    let out: BigInt = lhs * rhs;
    out.to_u128().unwrap()
}

/// ESTIMATES MAY BE INACCURATE UNLESS f64 IS CHANGED TO A PRECISSE DECIMAL REPRESENTATION 
pub fn calc_transition_prob(input_hamming: usize, output_hamming: usize, block_size: usize) -> f64 {
   let iowe = calc_iowe_entry(input_hamming, output_hamming, block_size) as f64;
   let w = BigInt::from_usize(input_hamming).unwrap();
   let n = BigInt::from_usize(block_size).unwrap();

   let denominator = binomial(n, w).to_u128().unwrap() as f64;
   iowe / denominator
}

/// Calculates column of IOWE matrix for the accumulate ode
pub fn calc_iowe_column(output_hamming: usize, block_size: usize) -> Vec<u128> {
    (0..block_size+1).map(|ih|{
        calc_iowe_entry(ih, output_hamming, block_size)
    }).collect_vec()
}

pub fn calc_transition_prob_column(output_hamming: usize, block_size: usize) -> Vec<f64> {
    (0..block_size+1).map(|ih|{
        calc_transition_prob(ih, output_hamming, block_size)
    }).collect_vec()
}

/// Calculates IOWE matrix in column-major order  for the accumulate ode
pub fn calc_iowe_matrix(block_size: usize) -> Vec<Vec<u128>> {
    (0..block_size+1).map(|ih|{
        calc_iowe_column(ih, block_size)
    }).collect_vec()
}

pub fn calc_transition_prob_matrix(block_size: usize) -> Vec<Vec<f64>> {
    (0..block_size+1).map(|ih|{
        calc_transition_prob_column(ih, block_size)
    }).collect_vec()
}

// pub fn multiple_acc_tprob_matrix<const BlockSize: usize>(n_acc: usize) -> SMatrix<f64, BlockSize, BlockSize> {
//     let columns = calc_transition_prob_matrix(BlockSize).iter().map(|col|{
//         SVector::from_vec(*col)
//     }).collect_vec();
//     let m = SMatrix::from_column_slice_generic(Const::<BlockSize>, Const::<BlockSize>, columns.as_slice());
//     m
// }
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn iowe_matrix() {
        let m = calc_iowe_matrix(3);
        // println!("C is {:?}", m);
        assert_eq!(m, vec![vec![1, 0, 0, 0], vec![0, 1, 2, 0], vec![0, 1, 1, 1], vec![0, 1, 0, 0]]);

        let c = calc_iowe_matrix(6);
        println!("c is {:?}", c);

        todo!("test against correct answer for c")
    }

    #[test]
    fn tprob_matrix() {
        let m = calc_transition_prob_matrix(3);
        assert_eq!(m, vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 0.3333333333333333, 0.6666666666666666, 0.0], vec![0.0, 0.3333333333333333, 0.3333333333333333, 1.0], vec![0.0, 0.3333333333333333, 0.0, 0.0]]);
    }

    #[test]
    fn test_batch_n_choose_k() {
        let m = n_choose_k_square_matrix(12);
        for (i, row) in m.iter().enumerate() {
            for (j, item) in row.iter().enumerate() {
                assert_eq!(
                    item.to_usize().unwrap(), 
                    binomial(i, j)
                );
            }
        }
    }
}