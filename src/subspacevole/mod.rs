use std::{usize, result};
use anyhow::{Error, anyhow};
use num_traits::ToBytes;
use rand::{Rng, SeedableRng};
use rand::rngs::{ThreadRng, StdRng};
use crate::{Fr, FrVec, FrMatrix, NUM_VOLES};
use crate::ff::Field;

// lazy_static! {
//     // pub static ref RAAA_CODE: RAAACode = RAAACode::deserialize(bytes)
// }

pub trait LinearCode {
    fn k(&self) -> usize;
    fn n(&self) -> usize;
    fn encode(&self, vec: &FrVec) -> FrVec;
    fn encode_extended(&self, vec: &FrVec) -> FrVec;
    fn check_parity(&self, putative_codeword: &FrVec) -> bool;
    fn check_parity_batch(&self, putative_codewords: &Vec<FrVec>) -> Result<(), Error> {
        match putative_codewords.iter().all(|pc|self.check_parity(pc)) {
            true => Ok(()),
            false => Err(anyhow!("Parity check failure"))
        }
    }

    fn mul_vec_by_extended_inverse(&self, u: &FrVec) -> FrVec;
        fn batch_encode_extended(&self, matrix: &Vec<FrVec>) -> Vec<FrVec> {
        matrix.iter().map(|x|self.encode_extended(x)).collect()
    }
    /// Calculates the prover's correction value for the whole U matrix
    fn mul_matrix_by_extended_inverse(&self, old_us: &FrMatrix) -> Vec<FrVec> {
        old_us.0.iter().map(|u|self.mul_vec_by_extended_inverse(u)).collect()
    }
    /// Returns (U, C) where U is the prover's correct U and C the correction value to send to verifier
    /// k is the dimension of the code
    fn get_prover_correction(&self, old_us: &FrMatrix) -> (FrMatrix, FrMatrix) {
        let l = old_us.0[0].0.len();
        let start_idx = self.k();
        let full_size = self.mul_matrix_by_extended_inverse(old_us);
        
        (
            FrMatrix(full_size.iter().map(|u|FrVec(u.0[0..start_idx].to_vec())).collect()),
            FrMatrix(full_size.iter().map(|u|FrVec(u.0[start_idx..].to_vec())).collect())
        )
    }

    /// Corrects the verifier's Q matrix give the prover's correction
    fn correct_verifier_qs(&self, old_qs: &FrMatrix, deltas: &FrVec, correction: &FrMatrix) -> FrMatrix {
        // Concatenate zero matrix with C as in the subsapace VOLE protocol:
        let l = old_qs.0[0].0.len();
        let correction_len = correction.0[0].0.len();
        let zero_len = l - correction_len;
        let zeroes_cons_c = (0..old_qs.0.len()).map(|i|{
            let mut out = Vec::with_capacity(l);
            out.append(&mut vec![Fr::ZERO; zero_len]);
            out.append(&mut correction.0[i].0.clone());
            FrVec(out)
        }).collect::<Vec<FrVec>>();
        
        let times_extended_generator = self.batch_encode_extended(&zeroes_cons_c);
        let times_deltas = times_extended_generator.iter().map(|x| {
            x * deltas
            //x.0.iter().zip(deltas.0.iter()).map(|(a, b)| a * b)
        }).collect::<Vec<FrVec>>();
        FrMatrix(
            old_qs.0.iter().zip(&times_deltas).map(|(q, t)|q - t).collect()
        )
    }
    /// `challenge_hash`` is the universal hash 
    /// `consistency_check` is the value returned frmo `calc_consistency_check`
    /// `deltas` and `q` are the verifier's deltas and q
    /// encoder
    /// WARNING If Using a smaller field, it may be important to use a challenge matrix instead of vector for sufficient security! 
    /// TODO: generics instead of RAAACode. And ofc generics for field
    fn verify_consistency_check(&self, challenge_hash: &FrVec, consistency_check: &(FrVec, FrVec), deltas: &FrVec, q_cols: &FrMatrix) -> Result<(), Error> {
        let u_hash = &consistency_check.0;
        let v_hash = &consistency_check.1;
        let q_hash = challenge_hash * q_cols;
        let u_hash_x_generator_x_diag_delta = &self.encode(u_hash) * deltas;
        if (*v_hash != &q_hash - &u_hash_x_generator_x_diag_delta) {
            Err(anyhow!("Consistency check fail!"))
        } else {
            Ok(())
        }
    }

}


#[derive(Debug, PartialEq)]
pub struct RAAACode {
    /// Forward and reverse permutations required for interleave and inverting interleave each time
    /// In order of when the interleaves are applied (e.g. 0th is after repetition and 2nd is before final accumulation)
    pub permutations: [(Vec<usize>, Vec<usize>); 3],
    /// Codeword length over dimension (rate's inverse). Default 2
    /// Exercise caution when changing q as this will affect the minimum distance and therefore security. Default q was selected for roughtly 128 bits of security at block length Fr,
    /// But THIS SECURITY CALCULATION WAS NOT DONE EXTREMELY RIGOROUSLY, rather by glancing at charts on "Coding Theorems for Repeat Multiple
    /// Accumulate Codes" by Kliewer et al
    /// A punctured code will likely perform better for the same security; the standard, unpuctured 1/2 rate RAAA code is used for its simplicity before choosing better codes.
    /// Furthermore, I have not sufficiently analyzed the security of using these binary RAAA codes on prime fields but 
    /// I would imagine it is fine as prime fields do not seem to make outputting a low-hamming-weight vector (and thus reducing distance of the code) any easier than doing so would be in GF2.
    pub q: usize
}
impl RAAACode {
    pub fn repeat(input: &FrVec, num_repeats: usize) -> FrVec {
        let mut out = Vec::with_capacity(num_repeats * input.0.len());
        for _ in 0..num_repeats {
            out.append(&mut input.0.clone());
        }
        FrVec(out)
    }
    // pub fn repeat_inverse(input: &FrVec, num_repeats: usize) -> FrVec {
    //     assert!(input.0.len() % num_repeats == 0, "input length must be divisible by num_repeats");
    //     let out_len = input.0.len() / num_repeats;
    //     FrVec(input.0[0..out_len].to_vec())
    // }
    
    /// Represents an invertible repetition matrix with extra rows so it's nxn rather than kxn
    /// The rows are choosen to be linearly independent and sparse so it is invertible and easy to invert
    /// Since repetition matrix is just 
    /// [ 1 0 0 0 1 0 0 0 ]
    /// [ 0 1 0 0 0 1 0 0 ]
    /// [ 0 0 1 0 0 0 1 0 ]
    /// [ 0 0 0 1 0 0 0 1 ]
    /// Extended reptition can be simply
    /// [ 1 0 0 0 1 0 0 0 ]
    /// [ 0 1 0 0 0 1 0 0 ]
    /// [ 0 0 1 0 0 0 1 0 ]
    /// [ 0 0 0 1 0 0 0 1 ]
    /// [ 0 0 0 0 1 0 0 0 ]
    /// [ 0 0 0 0 0 1 0 0 ]
    /// [ 0 0 0 0 0 0 1 0 ]
    /// [ 0 0 0 0 0 0 0 1 ]
    /// with its inverse sparsely described as 
    /// [ 1 0 0 0 -1 0 0 0 ]
    /// [ 0 1 0 0 0 -1 0 0 ]
    /// [ 0 0 1 0 0 0 -1 0 ]
    /// [ 0 0 0 1 0 0 0 -1 ]
    /// [ 0 0 0 0 1 0 0 0 ]
    /// [ 0 0 0 0 0 1 0 0 ]
    /// [ 0 0 0 0 0 0 1 0 ]
    /// [ 0 0 0 0 0 0 0 1 ]
    /// The sparse desciptions of these are respectively "clone the input vector and add its latter half with its latter half + first half"
    /// and "clone the input vector and replace its latter half with its latter half - first half"
    ///
    /// Exmaple for q = 3:
    /// [ 1 0 0 1 0 0 1 0 0 ]
    /// [ 0 1 0 0 1 0 0 1 0 ]
    /// [ 0 0 1 0 0 1 0 0 1 ]
    /// 
    /// extended = 
    /// [ 1 0 0 1 0 0 1 0 0 ]
    /// [ 0 1 0 0 1 0 0 1 0 ]
    /// [ 0 0 1 0 0 1 0 0 1 ]
    /// [ 0 0 0 1 0 0 0 0 0 ]
    /// [ 0 0 0 0 1 0 0 0 0 ]
    /// [ 0 0 0 0 0 1 0 0 0 ]
    /// [ 0 0 0 0 0 0 1 0 0 ]
    /// [ 0 0 0 0 0 0 0 1 0 ]
    /// [ 0 0 0 0 0 0 0 0 1 ]
    /// 
    ///  extended and inverted = 
    /// [ 1 0 0 -1 0 0 -1 0 0 ]
    /// [ 0 1 0 0 -1 0 0 -1 0 ]
    /// [ 0 0 1 0 0 -1 0 0 -1 ]
    /// [ 0 0 0 1 0 0 0 0 0 ]
    /// [ 0 0 0 0 1 0 0 0 0 ]
    /// [ 0 0 0 0 0 1 0 0 0 ]
    /// [ 0 0 0 0 0 0 1 0 0 ]
    /// [ 0 0 0 0 0 0 0 1 0 ]
    /// [ 0 0 0 0 0 0 0 0 1 ]
    /// The sparse description of these matrices are "clone the input vector. Let its 2nd and 3rd thirds += its first third"
    /// And "clone the input vector. Let its 2nd and 3rd thirds -= its first thirds."
    pub fn repeat_extended(input: &FrVec, q: usize) -> FrVec {
        let len = input.0.len();
        assert!(len % q == 0, "length must be divisible by q");
        let section_len = len / q;
        let zeroth_section = FrVec(input.0[0..section_len].to_vec());
        let mut out = Vec::with_capacity(len);
        out.append(&mut zeroth_section.0.clone());
        for i in 1..q {
            let start_idx = section_len * i;
            let new = &mut (&zeroth_section + 
                        &FrVec(input.0[start_idx..start_idx+section_len].to_vec()));
            out.append(&mut new.0);
        }
        FrVec(out)

    }

    pub fn repeat_extended_inverse(input: &FrVec, q: usize) -> FrVec {
        let len = input.0.len();
        assert!(len % q == 0, "length must be divisible by q");
        let section_len = len / q;
        let zeroth_section = FrVec(input.0[0..section_len].to_vec());
        let mut out = Vec::with_capacity(len);
        out.append(&mut zeroth_section.0.clone());
        for i in 1..q {
            let start_idx = section_len * i;
            let new = &mut (
                &FrVec(input.0[start_idx..start_idx+section_len].to_vec()) - 
                &zeroth_section
            );
            out.append(&mut new.0);
        }
        FrVec(out)
    }

    /// Permutation is not checked to be uniform. It simply contains a vec of new indices
    /// Interleave inverse is just interleave with the inverse of `permutation`
    pub fn interleave(input: &FrVec, permutation: &Vec<usize>) -> FrVec {
        let len = input.0.len();
        assert!(len == permutation.len(), "input length {} must match number of swaps {}", len, permutation.len());
        let mut out = vec![Fr::ZERO; len];
        for i in 0..len {
            out[permutation[i]] = input.0[i];
        }
        FrVec(out)
    }

    pub fn accumulate(input: &FrVec) -> FrVec {
        let l = input.0.len();
        let mut out = Vec::with_capacity(l);
        out.push(input.0[0]);
        for i in 1..l {
            out.push(input.0[i] + out[i-1]);
        }
        // let out = input.0.iter().reduce(|a, b| &(*a + b)).unwrap(); // Shouldn't panic because its simply addition...
        FrVec(out)
    }
    pub fn accumulate_inverse(input: &FrVec) -> FrVec {
        let l = input.0.len();
        let mut out = Vec::with_capacity(l);
        out.push(input.0[0]);
        for i in 1..l {
            out.push(input.0[i] - input.0[i-1]);
        }
        // let out = input.0.iter().reduce(|a, b| &(*a + b)).unwrap(); // Shouldn't panic because its simply addition...
        FrVec(out)
    }
    /// Returns a uniform permutation and its inverse
    /// It will be deterministic if and only if a seed is provided
    pub fn random_interleave_permutations(len: usize, seed: Option<[u8; 32]>) -> (Vec<usize>, Vec<usize>) {
        let range = 0..len;
        let mut rng = match seed { 
            Some(s) => StdRng::from_seed(s),
            None => StdRng::from_rng(&mut rand::thread_rng()).unwrap()
        };
        // let mut rng = ThreadRng::default();
        let mut forward = Vec::with_capacity(len);
        let mut backward = vec![0; len];
        let mut avail_indices = range.clone().collect::<Vec<usize>>();
        for i in range {
            let remove_idx = rng.gen_range(0..avail_indices.len());
            let removed_value = avail_indices.swap_remove(remove_idx);
            forward.push(removed_value);
            backward[removed_value] = i;
        }
        (forward, backward)
    }

    /// Creates an RAAA code of the default parameters 
    pub fn rand_default() -> RAAACode {
        
        let interleave_seeds = (0..3).map(|i|{
            *blake3::hash(
            format!("VOLE in the head RAAA code interleave {}", i).as_bytes()
            ).as_bytes()
        }).collect::<Vec<[u8; 32]>>();
        let permutations = [
            RAAACode::random_interleave_permutations(NUM_VOLES, Some(interleave_seeds[0])),
            RAAACode::random_interleave_permutations(NUM_VOLES, Some(interleave_seeds[1])),
            RAAACode::random_interleave_permutations(NUM_VOLES, Some(interleave_seeds[2])),
        ];
        RAAACode { permutations, q: 2 }
    }
    /// For testing. Note that block size under roughly 1024 for current code may not give 128 bits of security
    pub fn rand_with_parameters(block_size: usize, q: usize) -> Self {
        let permutations = [
            RAAACode::random_interleave_permutations(block_size, None),
            RAAACode::random_interleave_permutations(block_size, None),
            RAAACode::random_interleave_permutations(block_size, None),
        ];
        RAAACode { permutations, q }
    }
    /// Returns an array of u8s. Every four u8s represents a little-endian value. While these values are usizes for indexing, they should be small.
    /// If a usize go beyond the max u32 value, this returns an error.
    /// Codes should not be so large that they overflow a u32 so it is unlikely this will return an error.
    /// The four-byte chunks are as follows
    /// 0th: Number of repetitions for the repetition code, i.e. the code's `q` parameter
    /// 1st: Number of interleave*accumulates. For the foreseeable future 3 seems optimal and this is fixed at 3.
    /// 2nd: Length of codewords, i.e. the code's `n`
    /// [3rd , `n`+3rd): The first interleave permutation
    /// [`n`+3rd, `2n`+3rd): The first interleave permutation's inverse
    /// [2`n`+3rd , 3`n`+3rd): The second interleave permutation
    /// [3`n`+3rd , 4`n`+3rd): The second interleave permutation's inverse
    /// [4`n`+3rd , 5`n`+3rd): The third interleave permutation
    /// [5`n`+3rd , 6`n`+3rd): The third interleave permutation's inverse
    pub fn serialize(&self) -> Result<Vec<u8>, Error> {
        let mut usizes: Vec<usize> = Vec::with_capacity(
            3 + self.permutations[0].0.len() * 2 * self.permutations.len()
        );

        usizes.push(self.q.clone());
        usizes.push(self.permutations.len());
        usizes.push(self.permutations[0].0.len());
        self.permutations.iter().for_each(|x|{
            usizes.append(&mut x.0.clone());
            usizes.append(&mut x.1.clone());
        });

        if !usizes.iter().all(|u| *u <= u32::MAX as usize) {
            return Err(anyhow!("overflow"))
        }

        let mut u8s: Vec<u8> = Vec::with_capacity(usizes.len()*4);

        usizes.iter().for_each(|u|{
            u8s.append(&mut u.to_le_bytes()[0..4].to_vec())
        });
        
        Ok(u8s)
    }
    pub fn deserialize<T: AsRef<[u8]>>(bytes: T) -> Result<Self, Error> {
        let bytes = bytes.as_ref();
        if !(bytes.len() % 4 == 0) { return Err(anyhow!("input length must be divisible by 4")) }
        let l = bytes.len() / 4;

        let mut usizes = Vec::with_capacity(l);
        let mut idx_start = 0;
        for i_ in 0..l {
            usizes.push(u32::from_le_bytes(bytes[idx_start..idx_start+4].try_into().unwrap()) as usize);
            idx_start +=4;
        }  
        
        if usizes[1] != 3 { return  Err(anyhow!("only 3 interleaved accumulators are supported now")) }
        let nperms = usizes[1];
        let codeword_len = usizes[2];

        let perms: Vec<(Vec<usize>, Vec<usize>)> = (0..nperms).map(|i| {
            let start0 = 3 + i*codeword_len*2;
            let start1 = start0 + codeword_len;
            let end = start1 + codeword_len;
            (   
                // TODO: error instead of panic
                usizes.get(start0..start1).expect("Permutation is too short").to_vec(), 
                usizes.get(start1..end).expect("Permutation is too short").to_vec()
            )
        }).collect();
        Ok(Self {
            q: usizes[0],
            permutations: perms.try_into().unwrap() // Shouldn't panic since legnth is gauranteed 3
        })
    }
}

impl LinearCode for RAAACode {
    fn k(&self) -> usize {
        assert!(self.n() % self.q == 0, "n must be a multiple of q");
        return self.n() / self.q;
    }
    fn n(&self) -> usize {
        self.permutations[0].0.len()
    }
    /// Converts a vector to its codeword
    fn encode(&self, vec: &FrVec) -> FrVec {
        let repeated = Self::repeat(vec, self.q);
        let in0 = Self::interleave(&repeated, &self.permutations[0].0);
        let acc0 = Self::accumulate(&in0);
        let in1 = Self::interleave(&acc0, &self.permutations[1].0);
        let acc1 = Self::accumulate(&in1);
        let in2 = Self::interleave(&acc1, &self.permutations[2].0);
        let acc2 = Self::accumulate(&in2);

        acc2
    }

    /// Multiplies a single vector by the Tc matrix, the extended codeword generator to be invertible
    fn encode_extended(&self, vec: &FrVec) -> FrVec {
        let repeated = Self::repeat_extended(vec, self.q);
        let in0 = Self::interleave(&repeated, &self.permutations[0].0);
        let acc0 = Self::accumulate(&in0);
        let in1 = Self::interleave(&acc0, &self.permutations[1].0);
        let acc1 = Self::accumulate(&in1);
        let in2 = Self::interleave(&acc1, &self.permutations[2].0);
        let acc2 = Self::accumulate(&in2);

        acc2
    }

    /// Returns a single u vector multiplied by the Tc^-1 matrix (the extended generator matrix that is invertible). 
    fn mul_vec_by_extended_inverse(&self, u: &FrVec) -> FrVec { 
        let acc2_inv = Self::accumulate_inverse(&u);
        let in2_inv = Self::interleave(&acc2_inv, &self.permutations[2].1);
        let acc1_inv = Self::accumulate_inverse(&in2_inv);
        let in1_inv = Self::interleave(&acc1_inv, &self.permutations[1].1);
        let acc0_inv = Self::accumulate_inverse(&in1_inv);
        let in0_inv = Self::interleave(&acc0_inv, &self.permutations[0].1);         
        let re_inv = Self::repeat_extended_inverse(&in0_inv, self.q);

        re_inv
    }

    /// SECURITY TODO: (for audit?) check this is sufficient for determining whether something is a RAAA codeword
    /// For partity check, you can invert the accumulations and permutations and then check the result is in the subspace of the repetition code
    fn check_parity(&self, putative_codeword: &FrVec) -> bool {
        // Invet all the operations until the initial repetition code
        let acc2_inv = Self::accumulate_inverse(&putative_codeword);
        let in2_inv = Self::interleave(&acc2_inv, &self.permutations[2].1);
        let acc1_inv = Self::accumulate_inverse(&in2_inv);
        let in1_inv = Self::interleave(&acc1_inv, &self.permutations[1].1);
        let acc0_inv = Self::accumulate_inverse(&in1_inv);
        let should_be_repeated = Self::interleave(&acc0_inv, &self.permutations[0].1);     
       // Check that the reuslt is a codeword for the repetition code
        let len = should_be_repeated.0.len();
        assert!(len % self.q == 0, "length must be divisible by q");
        let section_len = len / self.q;
        assert!(self.q > 1, "can't check parity without repetition");
        let zeroth_section = should_be_repeated.0[0..section_len].to_vec();
        for i in 1..self.q {
            let idx_start = section_len * i;
            if should_be_repeated.0[idx_start ..idx_start + section_len].to_vec() != zeroth_section {
                return false;
            }
        }

        true

    }
}


/// `challenge_hash`` is the universal hash 
/// `u` and `v` are the prover's u and v values
/// WARNING If Using a smaller field, it may be important to use a challenge matrix instead of vector for sufficient security! 
/// Returns (challenge_hash*u, challenge_hash*v)
/// 
pub fn calc_consistency_check(challenge_hash: &FrVec, u_cols: &FrMatrix, v_cols: &FrMatrix) -> (FrVec, FrVec) {
    (challenge_hash * u_cols, challenge_hash * v_cols)
}

#[cfg(test)]
mod test {
    use std::{ops::Mul, time::Instant, io::repeat};

    use ff::{Field, PrimeField};
    use itertools::izip;
    use nalgebra::{Matrix2x4, Matrix4x2};
    use rand::rngs::ThreadRng;

    use crate::{FrRepr, smallvole::{self, VOLE, TestMOLE}};

    use super::*;

    #[test]
    fn test_serialize_deserialize() {
        let code = RAAACode {
            permutations: [RAAACode::random_interleave_permutations(6, None), RAAACode::random_interleave_permutations(6, None), RAAACode::random_interleave_permutations(6, None)],
            q: 2
        };
        // let code = RAAACode::rand_default();
        let s = code.serialize().unwrap();
        let d = RAAACode::deserialize(&s).unwrap();
        assert!(d == code);
    }


    #[test]
    fn test_permutation_and_inverse() {
        let (forward, backward) = RAAACode::random_interleave_permutations(5, None);
        let input = (0..5).map(|_|Fr::random(&mut ThreadRng::default())).collect();
        let input = FrVec(input);
        let permuted = RAAACode::interleave(&input, &forward);
        let inverse_permuted = RAAACode::interleave(&permuted, &backward);
        println!("input:\n{:?}\n\nforward:\n{:?}\n\nbackward:\n{:?}\n\npermute:\n{:?}\n\ninverse_permute:\n{:?}\n\n", 
            &input,
            &forward,
            &backward,
            &permuted,
            &inverse_permuted
        );
        assert_eq!(input, inverse_permuted);
    }
    #[test]
    fn test_accumulate_and_inverse() {
        let test0 = FrVec(vec![Fr::ZERO; 5]);
        let test1 = FrVec(vec![Fr::ONE; 5]);
        let test2 = FrVec(vec![Fr::ZERO, Fr::ONE, Fr::from_u128(2), Fr::from_u128(3)]);
        let test3 = FrVec(vec![Fr::ZERO, Fr::from_u128(2), Fr::ONE, Fr::from_u128(3)]);

        assert_eq!(RAAACode::accumulate(&test0).0, vec![Fr::ZERO; 5]);
        assert_eq!(RAAACode::accumulate(&test1).0, vec![Fr::ONE, Fr::from_u128(2), Fr::from_u128(3), Fr::from_u128(4), Fr::from_u128(5)]);
        assert_eq!(RAAACode::accumulate(&test2).0, vec![Fr::ZERO, Fr::ONE, Fr::from_u128(3), Fr::from_u128(6), ]);
        assert_eq!(RAAACode::accumulate(&test3).0, vec![Fr::ZERO, Fr::from_u128(2), Fr::from_u128(3), Fr::from_u128(6)]);
        vec![test0, test1, test2, test3].iter().for_each(|test|{
            let should_be_test = RAAACode::accumulate_inverse(&RAAACode::accumulate(test));
            assert_eq!(test.0, should_be_test.0);
        })
    }
    #[test]
    fn test_repeat() {
        let test0 = FrVec(vec![Fr::ZERO]);
        let test1 = FrVec(vec![Fr::ONE]);
        let test2 = FrVec(vec![Fr::from_u128(10), Fr::from_u128(11), Fr::from_u128(123456)]);
        assert_eq!(RAAACode::repeat(&test0, 2).0, vec![Fr::ZERO, Fr::ZERO]);
        assert_eq!(RAAACode::repeat(&test0, 3).0, vec![Fr::ZERO, Fr::ZERO, Fr::ZERO]);
        assert_eq!(RAAACode::repeat(&test1, 2).0, vec![Fr::ONE, Fr::ONE]);
        assert_eq!(RAAACode::repeat(&test1, 3).0, vec![Fr::ONE, Fr::ONE, Fr::ONE]);
        assert_eq!(RAAACode::repeat(&test2, 2).0, vec![Fr::from_u128(10), Fr::from_u128(11), Fr::from_u128(123456), Fr::from_u128(10), Fr::from_u128(11), Fr::from_u128(123456)]);
    }

    #[test]
    /// Tests the extended repetition matrix and its inverse
    /// TODO: more test cases
    fn test_repeat_extend_and_inverse() {
        let in0 = vec![Fr::from_u128(0), Fr::from_u128(1), Fr::from_u128(2), Fr::from_u128(3), Fr::from_u128(4), Fr::from_u128(5)];
        let in1 = vec![Fr::from_u128(5), Fr::from_u128(10), Fr::from_u128(15), Fr::from_u128(20), Fr::from_u128(25), Fr::from_u128(30)];

        let out0 = RAAACode::repeat_extended(&FrVec(in0.clone()), 2);
        let out1 = RAAACode::repeat_extended(&FrVec(in1.clone()),3);
        
        let inverted0 = RAAACode::repeat_extended_inverse(&out0, 2);
        let inverted1 = RAAACode::repeat_extended_inverse(&out1,3);

        assert_eq!(in0, inverted0.0);
        assert_eq!(in1, inverted1.0);
    }

    // TODO comprehesnvie test cases
    #[test]
    fn test_extended_encode() {
        let input = vec![Fr::from_u128(1), Fr::from_u128(5), Fr::from_u128(10), Fr::from_u128(0)];
        let code = RAAACode::rand_with_parameters(4, 2);
        let codeword = code.encode_extended(&FrVec(input.clone()));
        let inverse = code.mul_vec_by_extended_inverse(&codeword);
        assert_eq!(input, inverse.0);
    }
    
    #[test]
    fn test_prover_correction() {
        let test_mole = TestMOLE::init([123u8; 32], 16, 1024);
        // Check (at least one of the) VOLEs (and therefore likely all of them) was successful
        assert!(
            izip!(&test_mole.prover_outputs[7].u.0, &test_mole.prover_outputs[7].v.0, &test_mole.verifier_outputs[7].q.0).all(
                |(u, v, q)| u.clone() * test_mole.verifier_outputs[7].delta + v == q.clone()
            )
        );

        let u_cols = FrMatrix(test_mole.prover_outputs.iter().map(|o|o.u.clone()).collect::<Vec<FrVec>>());
        let v_cols = FrMatrix(test_mole.prover_outputs.iter().map(|o|o.v.clone()).collect::<Vec<FrVec>>());
        let q_cols = FrMatrix(test_mole.verifier_outputs.iter().map(|o|o.q.clone()).collect::<Vec<FrVec>>());
        let deltas = FrVec(test_mole.verifier_outputs.iter().map(|o|o.delta.clone()).collect());
        
        let u_rows = u_cols.transpose();
        let v_rows = v_cols.transpose();
        let q_rows = q_cols.transpose();

        let code = RAAACode::rand_default();
        // println!("single u length per VOLE {:?}\nNumber of deltas {:?}\n Single q length per VOLE{:?}. ", &test_mole.prover_outputs[0].u.len(), &test_mole.verifier_outputs.len(), &test_mole.verifier_outputs[0].q.len());

        let (new_us, correction) = code.get_prover_correction(
            &u_rows
        );

        let new_qs = code.correct_verifier_qs(
            &q_rows,
            &deltas,
            &correction
        );

        // check that (at least one of the) subspace VOLEs (and therefore likely all of them) is a successfull subspace VOLE:
        assert!(code.encode(&new_us.0[15]) * deltas + v_rows.0[15].clone() == new_qs.0[15]);
    }

    // TODO: more edge cases
    #[test]
    fn check_parity() {
        let code = RAAACode {
            permutations: [RAAACode::random_interleave_permutations(6, None), RAAACode::random_interleave_permutations(6, None), RAAACode::random_interleave_permutations(6, None)],
            q: 2
        };
        let input = FrVec::random(3);
        // let code = RAAACode::rand_default();
        // let input = FrVec::random(512);
        let codeword = code.encode(&input);
        let mut invalid_codeword = codeword.clone();
        invalid_codeword.0[2] = Fr::random(&mut rand::thread_rng());
        let mut invalid_length = codeword.clone();
        invalid_length.0.push(Fr::random(&mut rand::thread_rng()));
        assert!(code.check_parity(&codeword));
        assert!(!code.check_parity(&invalid_codeword));
        // assert!(!code.check_parity(&invalid_length));
    }
    #[test]
    fn check_parity_batch() {
        let code = RAAACode::rand_default();
        let input: Vec<FrVec> = (0..10).map(|_|FrVec::random(512)).collect();
        let mut codewords: Vec<FrVec> = input.iter().map(|x|code.encode(x)).collect();
        assert!(code.check_parity_batch(&codewords).is_ok());
        codewords[2].0[7] = Fr::random(&mut rand::thread_rng());
        assert!(code.check_parity_batch(&codewords).is_err())

    }
    #[test]
    fn consistency_check() {
        todo!()
    }
    
}