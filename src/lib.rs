pub mod subspacevole;
pub mod smallvole;
pub mod vecccom;

use ff::Field;
use nalgebra::ClosedMul;
// use num_traits::Zero;
#[macro_use]
extern crate ff;
use crate::ff::PrimeField;

#[derive(PrimeField)]
#[PrimeFieldModulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[PrimeFieldGenerator = "7"]
// Important this matches the endianness of MODULUS_AS_U128s
#[PrimeFieldReprEndianness = "big"]
pub struct Fr([u64; 4]);

impl polynomen::One for Fr {
    fn one() -> Self {
        Fr::ONE
    }
    fn is_one(&self) -> bool {
       self.eq(&Fr::ONE)
    }
}


impl polynomen::Zero for Fr {
    fn zero() -> Self {
        Fr::ZERO
    }
    fn is_zero(&self) -> bool {
       self.eq(&Fr::ZERO)
    }
}

impl num_traits::One for Fr {
    fn one() -> Self {
        Fr::ONE
    }
    fn is_one(&self) -> bool {
       self.eq(&Fr::ONE)
    }
}
impl num_traits::Zero for Fr {
    fn zero() -> Self {
        Fr::ZERO
    }
    fn is_zero(&self) -> bool {
       self.eq(&Fr::ZERO)
    }
}

// impl ClosedMul for Fr {}
// impl std::ops::Mul for Fr {
//     type Output = Self;
//     fn mul(self, rhs: Self) -> Self {
//         let mut out = self;
//         out.mul_assign(rhs);
//         out
//     }
// }
// impl MulAssign for Fr {
//     fn mul_assign(&mut self, rhs: Self) {
//         *self = *self * rhs;
//     }
// }
// // Crates such as nalgebra do not allow 
// fn matmul_fr() {

// }