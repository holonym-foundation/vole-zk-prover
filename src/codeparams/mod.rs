//! Calculates a lower bound for security for a Repeat Multiple Accumulate linear code
//! Note this underestimates the security. It simply provides
//! A lower bound on the probaility of getting a code with minimum distance d.
//! It is based off of 
//! "The Serial Concatenation of Rate-1 Codes Through Uniform Random Interleavers" by Henry D. Pfister and Paul H. Siegel [year]
//! "Coding Theorems for Turbo-Like Codes" by Divsalar [year]
//! And a proof the security should be no lower when used prime-field inputs rather than binary-field inputs which may be found in the Holonym V2 whitepaper
use itertools::Itertools;
use lazy_static::lazy_static;
use nalgebra::{SMatrix, Const, Vector, SVector};
use num_bigint::{BigInt, BigUint};
use num_traits::{FromPrimitive, ToPrimitive};

/// Batch n-choose-k, exploting patterns in the binomial coefficient for efficient batch calcualtion in a particular range
/// This only does it for square nxn matrices but coudlbe easily modified to suport abritrary matrix dimensions
pub fn n_choose_k_square_matrix(n: usize) -> Vec<Vec<BigUint>>{
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
pub fn calc_iowe_entry(input_hamming: usize, output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> u128 {
    if (input_hamming == 0) { 
        return if (output_hamming == 0) { 1 } else { 0 }
    } else if (output_hamming == 0) { 
        return 0
    }
    let w = input_hamming; //BigInt::from_usize(input_hamming).unwrap();
    let h = output_hamming; //BigInt::from_usize(output_hamming).unwrap();
    let n = block_size; //BigInt::from_usize(block_size).unwrap();
    // let two = BigInt::from_usize(2).unwrap();

    let lhs = &binomial_coeffs[n - &h][w.div_floor(2)]; // binomial(n - &h, w.div_floor(&two));
    let rhs = &binomial_coeffs[h - 1][w.div_ceil(2) - 1]; // binomial(h - 1, w.div_ceil(&two) - 1);

    let out: BigUint = lhs * rhs;
    out.to_u128().unwrap()
}

/// ESTIMATES MAY BE INACCURATE UNLESS f64 IS CHANGED TO A PRECISSE DECIMAL REPRESENTATION 
pub fn calc_transition_prob(input_hamming: usize, output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> f64 {
   let iowe = calc_iowe_entry(input_hamming, output_hamming, block_size, binomial_coeffs) as f64;
   let w = input_hamming; // BigInt::from_usize(input_hamming).unwrap();
   let n = block_size; // BigInt::from_usize(block_size).unwrap();

   let denominator = binomial_coeffs[n][w].to_f64().unwrap();// binomial(n, w).to_u128().unwrap() as f64;
   iowe / denominator
}

/// Calculates column of IOWE matrix for the accumulate ode
pub fn calc_iowe_column(output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> Vec<u128> {
    (0..block_size+1).map(|ih|{
        calc_iowe_entry(ih, output_hamming, block_size, binomial_coeffs)
    }).collect_vec()
}

pub fn calc_transition_prob_column(output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> DecimalVec {
    let v = (0..block_size+1).map(|ih|{
        calc_transition_prob(ih, output_hamming, block_size, binomial_coeffs)
    }).collect_vec();
    DecimalVec(v)
}

/// Calculates IOWE matrix in column-major order  for the accumulate ode
pub fn calc_iowe_matrix(block_size: usize) -> Vec<Vec<u128>> {
    let bcm = &n_choose_k_square_matrix(block_size);
    (0..block_size+1).map(|ih|{
        calc_iowe_column(ih, block_size, bcm)
    }).collect_vec()
}

pub fn calc_transition_prob_matrix(block_size: usize) -> DecimalMatrix {
    let bcm = &n_choose_k_square_matrix(block_size);
    let m = (0..block_size+1).map(|ih|{
        calc_transition_prob_column(ih, block_size, bcm)
    }).collect_vec();
    DecimalMatrix(m)
}

/// Calcualtes the transition probability
pub fn calc_multi_transition_prob_matrix(block_size: usize, num_accumulators: usize) -> DecimalMatrix {
    
}

#[derive(PartialEq, Debug)]
pub struct DecimalVec(pub Vec<f64>);
impl DecimalVec {
    pub fn dot(&self, rhs: &Self) -> f64 {
        self.0.iter().zip(rhs.0.iter()).map(|(a,b)| a*b).sum()
    }
}
#[derive(PartialEq, Debug)]
pub struct DecimalMatrix(pub Vec<DecimalVec>);
impl DecimalMatrix {
    pub fn mul(&self, rhs: &Self) -> Self {
        let transposed = rhs.transpose();
        let mut output_rows = Vec::with_capacity(self.0.len());
        
        for row in self.0.iter() {
            output_rows.push(
                DecimalVec(transposed.0.iter().map(|col|col.dot(row)).collect_vec())
            );
        }
        Self(output_rows)
    }
    pub fn transpose(&self) -> Self {
        let outer_len = self.0.len();
        let inner_len  = self.0[0].0.len();
        let mut res = Vec::with_capacity(inner_len);
        
        for i in 0..inner_len {
            let mut new = Vec::with_capacity(outer_len);
            for j in 0..outer_len {
                new.push(self.0[j].0[i]);
            }
            res.push(DecimalVec(new));
        }

        Self(res)
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use num_integer::binomial;

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
        let ans = DecimalMatrix(
            vec![
                DecimalVec(vec![1.0, 0.0, 0.0, 0.0]), 
                DecimalVec(vec![0.0, 0.3333333333333333, 0.6666666666666666, 0.0]), 
                DecimalVec(vec![0.0, 0.3333333333333333, 0.3333333333333333, 1.0]), 
                DecimalVec(vec![0.0, 0.3333333333333333, 0.0, 0.0])
            ]
        );
        assert_eq!(m, ans);
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

    #[test]
    fn matmul() {
        let a = DecimalMatrix(vec![
            DecimalVec(vec![ 1.0, 2.0, 3.0, 4.0, 5.0, 0.0 ]),
            DecimalVec(vec![ 0.0, 2.0, 3.0, 0.0, 5.0, 1.0 ]),
        ]);
        let b = DecimalMatrix(vec![
            DecimalVec(vec![ 1.0, 2.0, 3.0, ]),
            DecimalVec(vec![ 0.0, 5.0, 0.0 ]),
            DecimalVec(vec![ 0.0, 5.0, 0.0 ]),
            DecimalVec(vec![ 0.0, 2.0, 0.0 ]),
            DecimalVec(vec![ 0.0, 0.0, 1.0 ]),
            DecimalVec(vec![ 2.0, 5.0, 6.96969 ]),
        ]);

        let c = DecimalMatrix(vec![
            DecimalVec(vec![ 1.0, 35.0, 8.0 ]),
            DecimalVec(vec![ 2.0, 30.0, 11.96969]),

        ]);
        assert_eq!(a.mul(&b),c);
        todo!("Test edge cases")
    }
}