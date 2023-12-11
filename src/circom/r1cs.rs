//! Borrowed extensively from Nova Scotia https://github.com/nalinbhardwaj/Nova-Scotia/

use anyhow::{bail, Error};
use ff::PrimeField;
use itertools::Itertools;
use std::{io::{Read, Seek, SeekFrom}, collections::HashMap};
use byteorder::{LittleEndian, ReadBytesExt};

use crate::{Fr, FrVec, FrMatrix, zkp::{R1CS, R1CSWithMetadata}};

use super::read_constraint_vec;

type Witness = Vec<Fr>;

// R1CSFile's header
#[derive(Debug)]
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

#[derive(Debug)]
pub struct Constraints {
    a_rows: FrMatrix,
    b_rows: FrMatrix,
    c_rows: FrMatrix,
    a_wires: Vec<Vec<usize>>,
    b_wires: Vec<Vec<usize>>,
    c_wires: Vec<Vec<usize>>
}

#[derive(Debug)]
pub struct R1CSFile {
    pub version: u32,
    pub header: Header,
    pub constraints: Constraints,
    pub wire_mapping: Vec<u64>,
}

impl R1CSFile {
    /// Converts this to the R1CS format used by the rest of this crate
    pub fn to_crate_format(self) -> R1CSWithMetadata {
        println!("WIRE MAP IS {:?}", self.wire_mapping);
        // TODO: wire map
        let r1cs = R1CS {
            a_rows: self.constraints.a_rows,
            b_rows: self.constraints.b_rows,
            c_rows: self.constraints.c_rows
        };
        let pub_in_start = 1 + self.header.n_pub_out as usize;
        let public_outputs_indices = (1 .. pub_in_start).collect_vec();
        let public_inputs_indices = (pub_in_start .. pub_in_start + self.header.n_pub_in as usize).collect_vec();
        R1CSWithMetadata { r1cs, public_inputs_indices, public_outputs_indices }
    }
    /// Parses bytes in a circom .r1cs binary format
    pub fn from_reader<R: Read + Seek>(mut reader: R) -> Result<Self, Error> {
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
            bail!("This parser only supports 32-byte fields");
        }

        if header.prime_size != hex::decode("010000f093f5e1439170b97948e833285d588181b64550b829a031e1724e6430").unwrap() {
            bail!("This parser only supports bn254");
        }

        reader.seek(SeekFrom::Start(
            *section_offsets.get(&constraint_type).unwrap(),
        ))?;

        let constraints = read_constraints(
            &mut reader,
            *section_sizes.get(&constraint_type).unwrap(),
            &header,
        );

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

fn read_constraints<R: Read>(
    mut reader: R,
    size: u64,
    header: &Header,
) -> Constraints {
    
    let mut a_rows = Vec::with_capacity(header.n_constraints as usize);
    let mut b_rows = Vec::with_capacity(header.n_constraints as usize);
    let mut c_rows = Vec::with_capacity(header.n_constraints as usize);

    let mut a_wires = Vec::with_capacity(header.n_constraints as usize);
    let mut b_wires = Vec::with_capacity(header.n_constraints as usize);
    let mut c_wires = Vec::with_capacity(header.n_constraints as usize);


    for i in 0..header.n_constraints {
        let (a_wires_, a_row_) = read_constraint_vec(&mut reader);
        let (b_wires_, b_row_) = read_constraint_vec(&mut reader);
        let (c_wires_, c_row_) = read_constraint_vec(&mut reader);

        a_rows.push(FrVec(a_row_));
        b_rows.push(FrVec(b_row_));
        c_rows.push(FrVec(c_row_));

        a_wires.push(a_wires_);
        b_wires.push(b_wires_);
        c_wires.push(c_wires_);
    }
    let a_rows = FrMatrix(a_rows);
    let b_rows = FrMatrix(b_rows);
    let c_rows = FrMatrix(c_rows);

    Constraints { a_rows, b_rows, c_rows, a_wires, b_wires, c_wires }
}

fn read_map<R: Read>(mut reader: R, size: u64, header: &Header) -> Result<Vec<u64>, Error> {
    if size != header.n_wires as u64 * 8 {
        bail!("Invalid map section size");
    }
    let mut vec = Vec::with_capacity(header.n_wires as usize);
    for _ in 0..header.n_wires {
        vec.push(reader.read_u64::<LittleEndian>()?);
    }
    if vec[0] != 0 {
        bail!("Wire 0 should always be mapped to 0");
    }
    Ok(vec)
}

#[cfg(test)]
mod test {
    use std::{fs::{self, File}, io::BufReader};

    use super::*;
    #[test]
    fn read_r1cs_file() {
        let file = File::open("src/circom/examples/test.r1cs").unwrap();
        let mut buf_reader = BufReader::new(file);
        let r1cs = R1CSFile::from_reader(buf_reader).unwrap();
        println!("R1CS\n{:?}", r1cs);
    }
}