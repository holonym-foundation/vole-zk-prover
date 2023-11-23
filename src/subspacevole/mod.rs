use std::usize;
use rand::Rng;
use rand::rngs::ThreadRng;
use crate::{Fr, FrVec, FrMatrix};
use crate::ff::Field;

// lazy_static! {
//     // pub static ref RAAA_CODE: RAAACode = RAAACode::deserialize(bytes)
// }

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
        assert!(len == permutation.len(), "input length must match number of swaps");
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
    pub fn random_interleave_permutations(len: usize) -> (Vec<usize>, Vec<usize>) {
        let range = 0..len;
        let mut rng = ThreadRng::default();
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
        let permutations = [
            RAAACode::random_interleave_permutations(1024),
            RAAACode::random_interleave_permutations(1024),
            RAAACode::random_interleave_permutations(1024),
        ];
        RAAACode { permutations, q: 2 }
    }
    /// Block size under roughly 1024 for current code is unsafe 
    pub fn rand_with_parameters(block_size: usize, q: usize) -> Self {
        let permutations = [
            RAAACode::random_interleave_permutations(block_size),
            RAAACode::random_interleave_permutations(block_size),
            RAAACode::random_interleave_permutations(block_size),
        ];
        RAAACode { permutations, q }
    }
    pub fn serialize<T: AsRef<[u8]>>(&self) -> T {
        todo!()
    }
    pub fn deserialize<T>(bytes: T) -> Self {
        todo!()
    }

    /// Converts a vector to its codeword
    pub fn encode(&self, vec: &FrVec) -> FrVec {
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
    pub fn encode_extended(&self, vec: &FrVec) -> FrVec {
        let repeated = Self::repeat_extended(vec, self.q);
        let in0 = Self::interleave(&repeated, &self.permutations[0].0);
        let acc0 = Self::accumulate(&in0);
        let in1 = Self::interleave(&acc0, &self.permutations[1].0);
        let acc1 = Self::accumulate(&in1);
        let in2 = Self::interleave(&acc1, &self.permutations[2].0);
        let acc2 = Self::accumulate(&in2);

        acc2
    }

    pub fn batch_encode_extended(&self, matrix: &Vec<FrVec>) -> Vec<FrVec> {
        matrix.iter().map(|x|self.encode_extended(x)).collect()
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

    /// Calculates the prover's correction value for the whole U matrix
    fn mul_matrix_by_extended_inverse(&self, old_us: &FrMatrix) -> Vec<FrVec> {
        old_us.0.iter().map(|u|self.mul_vec_by_extended_inverse(u)).collect()
    }

    /// Returns (U, C) where U is the prover's correct U and C the correction value to send to verifier
    /// k is the dimension of the code
    pub fn get_prover_correction(&self, old_us: &FrMatrix) -> (FrMatrix, FrMatrix) {
        let l = old_us.0[0].0.len();
        assert!(l % self.q == 0, "block length is not a product of q");
        let start_idx = l / self.q;
        let full_size = self.mul_matrix_by_extended_inverse(old_us);
        
        (
            FrMatrix(full_size.iter().map(|u|FrVec(u.0[0..start_idx].to_vec())).collect()),
            FrMatrix(full_size.iter().map(|u|FrVec(u.0[start_idx..].to_vec())).collect())
        )
    }

    /// Corrects the verifier's Q matrix give the prover's correction
    pub fn correct_verifier_qs(&self, old_qs: &FrMatrix, deltas: &FrVec, correction: &FrMatrix) -> FrMatrix {
        // Concatenate zero matrix with C as in the subsapace VOLE protocol:
        let l = old_qs.0[0].0.len();
        assert!(l % self.q == 0, "block length is not a product of q");
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

    /// Challenge hash is the universal hash 
    /// WARNING If Using a smaller field, it may be important to use a challenge matrix instead of vector for sufficient security! 
    /// Returns (challenge_hash*u, challenge_hash*v)
    pub fn calc_prover_consistency_check(challenge_hash: &FrVec, u: &FrMatrix, v: &FrMatrix) -> (FrVec, FrVec) {
        (challenge_hash * u, challenge_hash * v)
    }
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
    fn test_permutation_and_inverse() {
        let (forward, backward) = RAAACode::random_interleave_permutations(5);
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

    #[test]
    fn consistency_check() {
        todo!()
    }
    
}