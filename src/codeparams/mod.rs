//! Calculates a lower bound for security for a Repeat Multiple Accumulate linear code
//! Note this underestimates the security. It simply provides
//! A lower bound on the probaility of getting a code with minimum distance d.
//! It is based off of 
//! "The Serial Concatenation of Rate-1 Codes Through Uniform Random Interleavers" by Henry D. Pfister and Paul H. Siegel [year]
//! "Coding Theorems for Turbo-Like Codes" by Divsalar [year]
//! And a proof the security should be no lower when used prime-field inputs rather than binary-field inputs which may be found in the Holonym V2 whitepaper

use std::str::FromStr;
use itertools::Itertools;
use bigdecimal::BigDecimal;
use nalgebra::{SMatrix, Const, Vector, SVector};
use num_bigint::{BigInt, BigUint};
use num_traits::{FromPrimitive, ToPrimitive};

/// This is easy: the IOWE of the repetition code. The rest of this file is for the accumulate code
/// rate is 1/q
pub fn repeat_iowe(block_size: usize, q: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> DecimalMatrix{
    let l = block_size+1;
    assert_eq!(block_size % q, 0, "block_size must be divisible by q");
    let mut rows = Vec::with_capacity(l);
    
    for w in 0..l {
        let mut row = Vec::with_capacity(l);

        for h in 0..l {
            if q * w == h {
                row.push(BigDecimal::from_str(&binomial_coeffs[block_size][w].to_string()).unwrap());
            } else {
                row.push(BigDecimal::from(0));
            }
        }
        rows.push(DecimalVec(row));

    }
    DecimalMatrix(rows)
}

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
pub fn calc_iowe_entry(input_hamming: usize, output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> BigUint {
    if (input_hamming == 0) { 
        return if (output_hamming == 0) { BigUint::from_u8(1).unwrap() } else { BigUint::from_u8(0).unwrap() }
    } else if (output_hamming == 0) { 
        return BigUint::from_u8(0).unwrap()
    }
    let w = input_hamming; //BigInt::from_usize(input_hamming).unwrap();
    let h = output_hamming; //BigInt::from_usize(output_hamming).unwrap();
    let n = block_size; //BigInt::from_usize(block_size).unwrap();
    // let two = BigInt::from_usize(2).unwrap();

    // floor and ceiling divided by 2 
    let floor = w / 2; let ceil = (w + 1) / 2;
    let lhs = &binomial_coeffs[n - &h][floor]; // binomial(n - &h, w.div_floor(&two));
    let rhs = &binomial_coeffs[h - 1][ceil - 1]; // binomial(h - 1, w.div_ceil(&two) - 1);

    lhs * rhs
}

/// ESTIMATES MAY BE INACCURATE UNLESS f64 IS CHANGED TO A PRECISSE DECIMAL REPRESENTATION 
pub fn calc_transition_prob(input_hamming: usize, output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> BigDecimal {
   let iowe = calc_iowe_entry(input_hamming, output_hamming, block_size, binomial_coeffs);
   let w = input_hamming; // BigInt::from_usize(input_hamming).unwrap();
   let n = block_size; // BigInt::from_usize(block_size).unwrap();

   let denominator = &binomial_coeffs[n][w];// binomial(n, w).to_u128().unwrap() as f64;
   BigDecimal::from_str(&iowe.to_string()).unwrap() / BigDecimal::from_str(&denominator.to_string()).unwrap()
}

/// Calculates column of IOWE matrix for the accumulate ode
pub fn calc_iowe_column(output_hamming: usize, block_size: usize, binomial_coeffs: &Vec<Vec<BigUint>>) -> Vec<BigUint> {
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
pub fn calc_iowe_matrix_cols(block_size: usize) -> Vec<Vec<BigUint>> {
    let bcm = &n_choose_k_square_matrix(block_size);
    (0..block_size+1).map(|ih|{
        calc_iowe_column(ih, block_size, bcm)
    }).collect_vec()
}

pub fn calc_transition_prob_matrix_cols(block_size: usize) -> DecimalMatrix {
    let bcm = &n_choose_k_square_matrix(block_size);
    let m = (0..block_size+1).map(|ih|{
        calc_transition_prob_column(ih, block_size, bcm)
    }).collect_vec();
    DecimalMatrix(m)
}
pub fn calc_transition_prob_matrix(block_size: usize) -> DecimalMatrix {
    calc_transition_prob_matrix_cols(block_size).transpose()
}

/// Calcualtes the transition probability
pub fn calc_multi_transition_prob_matrix(block_size: usize, num_accumulators: usize) -> DecimalMatrix {
    if num_accumulators < 1 { panic!("num_accumulators must be >= 1") }
    let pm = calc_transition_prob_matrix(block_size);
    let mut res = pm.clone();
    for _i in (1..num_accumulators) {
        res = res.mul(
            &pm.clone()
        );
    }
    res
}

/// Calculates the expected value of number of outputs with hamming weight h of an RMA code with rate 1/`q`, block size `block_size`, and `num_accumulators` rate-1 accumulators preceded by rate-1 interleavers
pub fn expected_num_outputs_with_weight(q: usize, block_size: usize, num_accumulators: usize, h: usize) -> BigDecimal {
    assert!(h > 0, "h must be > 0");
    assert!(block_size % q == 0, "block size must be divisible by q");
    let k = block_size / q;

    let bcm = &n_choose_k_square_matrix(block_size);
    let iowe_rep = repeat_iowe(block_size, q, &bcm);
    let outer_cols = iowe_rep.transpose();
    let pm = calc_multi_transition_prob_matrix(block_size, num_accumulators);
    let inner_cols = pm.transpose();

    let mut res: BigDecimal = BigDecimal::from(0);
    for i in 1..k+1 {
        // The expected number of outputs of Hamming weight h given input of Hamming weight i
        let for_input_weight_i = iowe_rep.0[i].dot(&inner_cols.0[h]);
        res += for_input_weight_i;
    }
    res
}

/// Calculates an upper *bound* on the probabilty that a RMA code with rate 1/`q`, block size `block_size` and `num_accumulators` rate-1 accumulators preceded by rate-1 interleavers has minimum distance less than `d`
/// Also returns all upper bounds on this probability for weights [1..d)
pub fn max_prob_distance_lt(q: usize, block_size: usize, num_accumulators: usize, d: usize) -> (BigDecimal, Vec<BigDecimal>) {
    let mut upper_bounds = Vec::<BigDecimal>::with_capacity(d-1);
    let mut upper_bound_d = BigDecimal::from(0);
    for i in 1..d {
        let a_h = expected_num_outputs_with_weight(q, block_size, num_accumulators, i);
        upper_bound_d += a_h;
        upper_bounds.push(upper_bound_d.clone());
    }
    (upper_bound_d, upper_bounds)
}

/// Entry point
pub fn main() {
    let d = 80;
    let (for_d, for_all_til_d) = max_prob_distance_lt(2, 64, 3, d);
    println!("prob of minimum distance under {} is at most {}", d, for_d);
    println!("prob of minimum distance under all numbers precenting {} are at most {:?}", d, for_all_til_d); 
}
#[derive(PartialEq, Debug, Clone)]
pub struct DecimalVec(pub Vec<BigDecimal>);
impl DecimalVec {
    pub fn dot(&self, rhs: &Self) -> BigDecimal {
        self.0.iter().zip(rhs.0.iter()).map(|(a,b)| a*b).sum()
    }
    pub fn from_f64_vec(v: Vec<f64>) -> Self {
        Self(v.iter().map(|x|BigDecimal::from_f64(*x).unwrap()).collect_vec())
    }
    pub fn is_close_to(&self, rhs: &Self, tol: f64) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a,b)| (a-b).abs() < BigDecimal::from_f64(tol).unwrap())
    }
}
#[derive(PartialEq, Debug, Clone)]
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
                new.push(self.0[j].0[i].clone());
            }
            res.push(DecimalVec(new));
        }

        Self(res)
    }
    pub fn is_close_to(&self, rhs: &Self, tol: f64) -> bool {
        self.0.iter().zip(rhs.0.iter()).all(|(a,b)| a.is_close_to(b, tol))
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use num_integer::binomial;

    #[test]
    fn iowe_matrix() {
        let m = calc_iowe_matrix_cols(3);
        // println!("C is {:?}", m);
        assert_eq!(m, vec![
            vec![1, 0, 0, 0].iter().map(|x|BigUint::from_u8(*x).unwrap()).collect_vec(),
            vec![0, 1, 2, 0].iter().map(|x|BigUint::from_u8(*x).unwrap()).collect_vec(),
            vec![0, 1, 1, 1].iter().map(|x|BigUint::from_u8(*x).unwrap()).collect_vec(),
            vec![0, 1, 0, 0].iter().map(|x|BigUint::from_u8(*x).unwrap()).collect_vec()
        ]
        );

        let c = calc_iowe_matrix_cols(6);
        todo!("test against correct answer for c")
    }

    #[test]
    fn tprob_matrix() {
        let m = calc_transition_prob_matrix_cols(3);
        let ans = DecimalMatrix(
            vec![
                DecimalVec::from_f64_vec(vec![1.0, 0.0, 0.0, 0.0]), 
                DecimalVec::from_f64_vec(vec![0.0, 0.3333333333333333, 0.6666666666666666, 0.0]), 
                DecimalVec::from_f64_vec(vec![0.0, 0.3333333333333333, 0.3333333333333333, 1.0]), 
                DecimalVec::from_f64_vec(vec![0.0, 0.3333333333333333, 0.0, 0.0])
            ]
        );
        assert!(m.is_close_to(&ans, 1e-10));
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
            DecimalVec::from_f64_vec(vec![ 1.0, 2.0, 3.0, 4.0, 5.0, 0.0 ]),
            DecimalVec::from_f64_vec(vec![ 0.0, 2.0, 3.0, 0.0, 5.0, 1.0 ]),
        ]);
        let b = DecimalMatrix(vec![
            DecimalVec::from_f64_vec(vec![ 1.0, 2.0, 3.0, ]),
            DecimalVec::from_f64_vec(vec![ 0.0, 5.0, 0.0 ]),
            DecimalVec::from_f64_vec(vec![ 0.0, 5.0, 0.0 ]),
            DecimalVec::from_f64_vec(vec![ 0.0, 2.0, 0.0 ]),
            DecimalVec::from_f64_vec(vec![ 0.0, 0.0, 1.0 ]),
            DecimalVec::from_f64_vec(vec![ 2.0, 5.0, 6.96969 ]),
        ]);

        let c = DecimalMatrix(vec![
            DecimalVec::from_f64_vec(vec![ 1.0, 35.0, 8.0 ]),
            DecimalVec::from_f64_vec(vec![ 2.0, 30.0, 11.96969]),

        ]);
        assert!(a.mul(&b).is_close_to(&c, 1e-10));
        todo!("Test edge cases")
    }
    #[test]
    fn transpose() {
        let x = DecimalMatrix(vec![
            DecimalVec::from_f64_vec(vec![ 1.0, 3.0, 5.0 ]),
            DecimalVec::from_f64_vec(vec![ 2.0, 4.0, 6.7]),

        ]);
        let y = DecimalMatrix(vec![
            DecimalVec::from_f64_vec(vec![ 1.0, 2.0 ]),
            DecimalVec::from_f64_vec(vec![ 3.0, 4.0 ]),
            DecimalVec::from_f64_vec(vec![ 5.0, 6.7 ]),

        ]);
        assert_eq!(x, y.transpose());
        assert_eq!(y, x.transpose());
    }
    #[test]
    fn multiple_acc() {
        let once = calc_multi_transition_prob_matrix(3, 1);
        assert!(
            once.is_close_to(
                &DecimalMatrix(
                    vec![
                        DecimalVec::from_f64_vec(vec![1.0, 0.0, 0.0, 0.0]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.6666666666666666, 0.3333333333333333, 0.0]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.0, 1.0, 0.0])
                    ]
                ),
                1e-10
            )
        );

        let twice = calc_multi_transition_prob_matrix(3, 2);
        assert!(
            twice.is_close_to( 
                &DecimalMatrix(
                    vec![
                        DecimalVec::from_f64_vec(vec![1.0, 0.0, 0.0, 0.0]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.3333333333333333, 0.5555555555555556, 0.1111111111111111]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.4444444444444444, 0.3333333333333333, 0.2222222222222222]), 
                        DecimalVec::from_f64_vec(vec![0.0, 0.6666666666666666, 0.3333333333333333, 0.])
                    ]
                ),
                1e-10
            )
        );
        // println!("Large probabilty mat {:?}", calc_multi_transition_prob_matrix(1024, 3));
    }

    #[test]
    fn repetition_iowe() {
        todo!("test against correct answer")
    }
}