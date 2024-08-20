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
    match reader.read_exact(&mut buf){
        Ok(num) => num,
        Err(e) => panic!("Problem opening the buffer: {e:?}"),
    };
    buf.reverse(); // Convert endianness to big
    Fr::from_repr(FrRepr(buf)).unwrap()
}

/// Reads l Frs from a circom file
/// I believe this should be more performant because it seems the compiler will be able to vectorize easily than doing multiple individual function calls
fn read_fr_vec<R: Read>(mut reader: R, l: usize) -> Vec<Fr> {
    let mut bufs = vec![[0u8; 32]; l];
    bufs.iter_mut().map(|buf|{
        match reader.read_exact(buf){
            Ok(num) => num,
            Err(e) => panic!("Problem opening the buffer: {e:?}"),           
        };
        buf.reverse();
        Fr::from_repr(FrRepr(*buf)).unwrap()
    }).collect()
}

/// Reads l u32 wire labels and corresponding Frs from a R1CS file
fn read_constraint_vec<R: Read>(mut reader: R) -> SparseVec<Fr> {
    let l = match reader.read_u32::<LittleEndian>(){
        Ok(num) => num,
        Err(e) => panic!("Problem reading u32: {e:?}"),
    } as usize;
    let mut constraints = Vec::with_capacity(l);
    for _ in 0..l {
        constraints.push(
            (
                match reader.read_u32::<LittleEndian>(){
                    Ok(num) => num,
                    Err(e) => panic!("Problem reading u32: {e:?}"),
                } as usize,
                {
                    let mut buf = [0u8; 32];
                    match reader.read_exact(&mut buf){
                        Ok(num) => num,
                        Err(e) => panic!("Problem opening the buffer: {e:?}"),
                    }; 
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
        let wtns_file_res = File::open("src/circom/examples/witness.wtns");
        let wtns_file = match wtns_file_res {
            Ok(file) => file,
            Err(e) => panic!("Problem opening the file: {e:?}"),
        };
        let wtns_reader = BufReader::new(wtns_file);
        let witness = match wtns_from_reader(wtns_reader){
            Ok(wtns) => wtns,
            Err(e) => panic!("Problem opening the reader: {e:?}"),  
        };

        let r1cs_file_res = File::open("src/circom/examples/test.r1cs");
        let r1cs_file = match r1cs_file_res {
            Ok(file) => file,
            Err(e) => panic!("Problem opening the file: {e:?}"),
        };
        let r1cs_reader = BufReader::new(r1cs_file);
        let r1cs = match R1CSFile::from_reader(r1cs_reader){
            Ok(r1cs) => r1cs,
            Err(e) => panic!("Problem opening the reader: {e:?}"),
        }.to_crate_format();

        let _x = match e2e_test(witness, r1cs){
            Ok(verif) => verif,
            Err(e) => panic!("e2e test not passing: {e:?}"),
        };
    }
}