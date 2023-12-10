//! Borrowed extensively from Nova Scotia https://github.com/nalinbhardwaj/Nova-Scotia/

use anyhow::{bail, Error};
use ff::PrimeField;
use std::{io::{Read, SeekFrom}, collections::HashMap};
use byteorder::{LittleEndian, ReadBytesExt};

use crate::{Fr, FrVec, FrRepr};

type Witness = Vec<Fr>;

// R1CSFile's header
#[derive(Debug, Default)]
pub struct Header {
    pub field_size: u32,
    pub prime_size: Vec<u8>,
    pub n_wires: u32,
    pub n_pub_out: u32,
    pub n_pub_in: u32,
    pub n_prv_in: u32,
    pub n_labels: u64,
    pub n_constraints: u32,
}
/// Parses bytes in a circom .r1cs binary format
fn r1cs_from_reader<R: Read + Seek>(mut reader: R) -> Result<FrVec, Error> {
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if magic != "r1cs".as_bytes() {
        bail!("Invalid magic number");
    }

    let version = reader.read_u32::<LittleEndian>()?;
    if version != 1 {
        bail!("Unsupported version")
    }

    let num_sections = reader.read_u32::<LittleEndian>()?;

    // section type -> file offset
    let mut section_offsets = HashMap::<u32, u64>::new();
    let mut section_sizes = HashMap::<u32, u64>::new();

    // get file offset of each section
    for _ in 0..num_sections {
        let section_type = reader.read_u32::<LittleEndian>()?;
        let section_size = reader.read_u64::<LittleEndian>()?;
        let offset = reader.seek(SeekFrom::Current(0))?;
        section_offsets.insert(section_type, offset);
        section_sizes.insert(section_type, section_size);
        reader.seek(SeekFrom::Current(section_size as i64))?;
    }

    let header_type = 1;
    let constraint_type = 2;
    let wire2label_type = 3;

    reader.seek(SeekFrom::Start(*section_offsets.get(&header_type).unwrap()))?;
    let header = read_header(&mut reader, *section_sizes.get(&header_type).unwrap())?;
    if header.field_size != 32 {
        bail!("This parser only supports 32-byet fields");
    }

    if header.prime_size != hex::decode("010000f093f5e1439170b97948e833285d588181b64550b829a031e1724e6430").unwrap() {
        bail!("This parser only supports bn254");
    }

    reader.seek(SeekFrom::Start(
        *section_offsets.get(&constraint_type).unwrap(),
    ))?;
    let constraints = read_constraints::<&mut R>(
        &mut reader,
        *section_sizes.get(&constraint_type).unwrap(),
        &header,
    )?;

    reader.seek(SeekFrom::Start(
        *section_offsets.get(&wire2label_type).unwrap(),
    ))?;
    let wire_mapping = read_map(
        &mut reader,
        *section_sizes.get(&wire2label_type).unwrap(),
        &header,
    )?;

    Ok(R1CSFile {
        version,
        header,
        constraints,
        wire_mapping,
    })
}

fn read_header<R: Read>(mut reader: R, size: u64) -> Result<Header, Error> {
    let field_size = reader.read_u32::<LittleEndian>()?;
    let mut prime_size = vec![0u8; field_size as usize];
    reader.read_exact(&mut prime_size)?;
    if size != 32 + field_size as u64 {
        bail!("Invalid header section size");
    }

    Ok(Header {
        field_size,
        prime_size,
        n_wires: reader.read_u32::<LittleEndian>()?,
        n_pub_out: reader.read_u32::<LittleEndian>()?,
        n_pub_in: reader.read_u32::<LittleEndian>()?,
        n_prv_in: reader.read_u32::<LittleEndian>()?,
        n_labels: reader.read_u64::<LittleEndian>()?,
        n_constraints: reader.read_u32::<LittleEndian>()?,
    })
}