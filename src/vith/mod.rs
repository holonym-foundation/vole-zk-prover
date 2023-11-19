// use crate::{FrMatrix, Fr};

// struct VoleInTheHeadTest {
//     /// A commitment to the prover's witness
//     pub witness_comm: FrMatrix,
//     /// S from the paper
//     pub s: FrMatrix,
//     /// U1 from the paper
//     pub u1: Fr,
//     /// R (U2) from the paper
//     pub r: Fr,
//     /// ∆' from paper
//     pub delta: Fr,
// }

// struct VITHValues {
//     /// A commitment to the prover's witness
//     pub witness_comm: FrMatrix,
//     /// U1 from the paper
//     pub u1: Fr,
//     /// R (U2) from the paper
//     pub r: Fr,

// }
// pub fn from_subspace_vole_to_vith(u_rows: &FrMatrix, v_rows: &FrMatrix, witness: &FrMatrix) -> VITHValues {
//     let num_u_rows = u_rows.0.len();
//     let num_v_rows = v_rows.0.len();
//     assert!(num_u_rows % 2 == 0, "U must have an even number of rows");
//     assert!(num_v_rows % 2 == 0, "V must have an even number of rows");
//     let u_halfway = num_u_rows / 2;
//     let v_halfway = num_v_rows / 2;
//     let u1 = u_rows.0[0..u_halfway].to_vec();
//     let r = u_rows.0[u_halfway..].to_vec();
    
// }