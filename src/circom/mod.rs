//! Modules for reading and writing witness and R1CS from circom format

use std::io::Read;

use anyhow::Error;
use byteorder::{LittleEndian, ReadBytesExt};
use ff::PrimeField;

use crate::{Fr, FrRepr, SparseVec};
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
fn read_constraint_vec<R: Read>(mut reader: R) -> SparseVec<Fr> {
    let l = reader.read_u32::<LittleEndian>().unwrap() as usize;
    let mut constraints = Vec::with_capacity(l);
    for _ in 0..l {
        constraints.push(
            (
                reader.read_u32::<LittleEndian>().unwrap() as usize,
                {
                    let mut buf = [0u8; 32];
                    reader.read_exact(&mut buf).unwrap(); 
                    buf.reverse();
                    Fr::from_repr(FrRepr(buf)).unwrap()
                }
            )
        )
    };
    SparseVec(constraints)
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};

    use crate::{circom::{witness::wtns_from_reader}, actors::test_helpers::e2e_test};

    use super::{*, r1cs::R1CSFile};
    #[test]
    fn e2e_r1cs_wtns_files() {
        let wtns_file = File::open("src/circom/examples/witness.wtns").unwrap();
        let mut wtns_reader = BufReader::new(wtns_file);
        let witness = wtns_from_reader(wtns_reader).unwrap();

        let r1cs_file = File::open("src/circom/examples/test.r1cs").unwrap();
        let mut r1cs_reader = BufReader::new(r1cs_file);
        let r1cs = R1CSFile::from_reader(r1cs_reader).unwrap().to_crate_format();

        assert!(e2e_test(witness, r1cs).is_ok());
        
    }
}