use anyhow::{bail, Error};
use ff::PrimeField;
use std::io::Read;
use byteorder::{LittleEndian, ReadBytesExt};

use crate::{FVec, Fr};

use super::{read_fr, read_fr_vec};

/// Parses bytes in a circom .wtns binary format
/// Borrowed extensively from Nova Scotia https://github.com/nalinbhardwaj/Nova-Scotia/blob/main/src/circom/reader.rs
pub fn wtns_from_reader<R: Read>(mut reader: R) -> Result<FVec<Fr>, Error> {
    let mut wtns_header = [0u8; 4];
    reader.read_exact(&mut wtns_header)?;
    if wtns_header != "wtns".as_bytes() {
        bail!("invalid file header");
    }
    let version = reader.read_u32::<LittleEndian>()?;
    // println!("wtns version {}", version);
    if version > 2 {
        bail!("unsupported file version");
    }
    let num_sections = reader.read_u32::<LittleEndian>()?;
    if num_sections != 2 {
        bail!("invalid num sections");
    }
    // read the first section
    let sec_type = reader.read_u32::<LittleEndian>()?;
    if sec_type != 1 {
        bail!("invalid section type");
    }
    let sec_size = reader.read_u64::<LittleEndian>()?;
    if sec_size != 4 + 32 + 4 {
        bail!("invalid section len")
    }
    let field_size = reader.read_u32::<LittleEndian>()?;
    if field_size != 32 {
        bail!("invalid field byte size");
    }
    let mut prime = vec![0u8; field_size as usize];
    reader.read_exact(&mut prime)?;
    if prime != hex::decode("010000f093f5e1439170b97948e833285d588181b64550b829a031e1724e6430").unwrap() {
        bail!("invalid curve prime {:?}", prime);
    }
    let witness_len = reader.read_u32::<LittleEndian>()?;
    // println!("witness len {}", witness_len);
    let sec_type = reader.read_u32::<LittleEndian>()?;
    if sec_type != 2 {
        bail!("invalid section type");
    }
    let sec_size = reader.read_u64::<LittleEndian>()?;
    if sec_size != (witness_len * field_size) as u64 {
        bail!("invalid witness section size {}", sec_size);
    }

    Ok(
        FVec::<Fr>(read_fr_vec(reader, witness_len as usize))
    )
}

#[cfg(test)]
mod test {
    use std::{fs::{self, File}, io::BufReader};

    use super::*;
    #[test]
    fn read_wtns_file() {
        let file = File::open("src/circom/examples/witness.wtns").unwrap();
        let mut buf_reader = BufReader::new(file);
        let witness = wtns_from_reader(buf_reader).unwrap();
        println!("Witness\n{:?}", witness.0);
        println!("Witness\n{}", witness);
    }
}