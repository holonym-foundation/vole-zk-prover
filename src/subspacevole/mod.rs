use nalgebra::{self, MatrixMN, SMatrix, RowSVector};
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

struct ReedSolomonCode;
impl ReedSolomonCode {
    // Cosntructs a Vandermode Matrix for N and K
    fn construct_generator<const N: usize, const K: usize>() -> SMatrix<Fr, K, N> {
        SMatrix::<Fr, K, N>::from_fn(|row, col| {
            let base = Fr::from(row as u64 + 1);
            base.pow([col as u64])
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn rs_is_vandermonde() {
        let g = ReedSolomonCode::construct_generator::<5,3>();
        println!("Gen is {:?}", g);
        assert_eq!(
            g,
            SMatrix::<Fr, 3, 5>::from_row_slice(&[
                Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64),
                Fr::from(1u64), Fr::from(2u64), Fr::from(4u64), Fr::from(8u64), Fr::from(16u64),
                Fr::from(1u64), Fr::from(3u64), Fr::from(9u64), Fr::from(27u64), Fr::from(81u64), 
            ])
        );
    }
}