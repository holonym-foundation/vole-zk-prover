//! Modules for reading and writing witness and R1CS from circom format

use std::io::Read;

use anyhow::Error;
use byteorder::{LittleEndian, ReadBytesExt};
use ff::PrimeField;

use crate::{Fr, FrRepr};
pub mod witness;
pub mod r1cs;

/// Reads an Fr from a circom file
fn read_fr<R: Read>(mut reader: R) -> Fr {
    let mut buf = [0u8; 32];
    reader.read_exact(&mut buf).unwrap();
    buf.reverse(); // Convert endianness to big
    Fr::from_repr(FrRepr(buf)).unwrap()
}

/// Reads l Frs from a circom file
/// I believe this should be more performant because it seems the compiler will be able to vectorize easily than doing multiple individual function calls
fn read_fr_vec<R: Read>(mut reader: R, l: usize) -> Vec<Fr> {
    let mut bufs = vec![[0u8; 32]; l];
    bufs.iter_mut().map(|buf|{
        reader.read_exact(buf).unwrap();
        buf.reverse();
        Fr::from_repr(FrRepr(*buf)).unwrap()
    }).collect()
}

/// Reads l u32 wire labels and corresponding Frs from a R1CS file
fn read_constraint_vec<R: Read>(mut reader: R, l: usize) -> (Vec<usize>, Vec<Fr>) {
    // let mut bufs = vec![[0u8; 32]; l];
    let mut wire_labels = Vec::with_capacity(l);
    let mut frs = Vec::with_capacity(l);
    for i in 0..l {
        wire_labels.push(
            reader.read_u32::<LittleEndian>().unwrap() as usize
        );

        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf).unwrap();
        buf.reverse();
        
        frs.push(
            Fr::from_repr(FrRepr(buf)).unwrap(),
        );
    };
    
    (wire_labels, frs)
}