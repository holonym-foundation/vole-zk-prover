//! Calculates a lower bound for security for a Repeat Multiple Accumulate linear code
//! Note this underestimates the security. It simply provides
//! A lower bound on the probaility of getting a code with minimum distance d.
//! It is based off of 
//! "The Serial Concatenation of Rate-1 Codes Through Uniform Random Interleavers" by Henry D. Pfister and Paul H. Siegel [year]
//! "Coding Theorems for Turbo-Like Codes" by Divsalar [year]
//! And a proof the security should be no lower when used prime-field inputs rather than binary-field inputs which may be found in the Holonym V2 whitepaper

use itertools::Itertools;
use num_bigint::BigInt;
use num_integer::{binomial, Integer};
use num_traits::{FromPrimitive, ToPrimitive};

/// Calculates the Input Output Weight Enumeration for an accumulate code with input hamming weight w, output hamming weight h, and block size n
fn calc_iowe_entry(input_hamming: usize, output_hamming: usize, block_size: usize) -> u128 {
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

/// Calculates column of IOWE matrix for the accumulate ode
fn calc_iowe_column(output_hamming: usize, block_size: usize) -> Vec<u128> {
    (0..block_size+1).map(|ih|{
        calc_iowe_entry(ih, output_hamming, block_size)
    }).collect_vec()
}

/// Calculates IOWE matrix in column-major order  for the accumulate ode
fn calc_iowe_matrix(block_size: usize) -> Vec<Vec<u128>> {
    (0..block_size+1).map(|ih|{
        calc_iowe_column(ih, block_size)
    }).collect_vec()
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn col() {
        let m = calc_iowe_matrix(3);
        // println!("C is {:?}", m);
        assert_eq!(m, vec![vec![1, 0, 0, 0], vec![0, 1, 2, 0], vec![0, 1, 1, 1], vec![0, 1, 0, 0]]);

        let c = calc_iowe_matrix(6);
        println!("c is {:?}", c);

        todo!("test against correct answer for c")
    }
}