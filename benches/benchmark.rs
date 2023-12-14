use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::{Field, PrimeField};
use lazy_static::lazy_static;
use nalgebra::SMatrix;
use rand::{rngs::ThreadRng, Rng};
use volonym::{vecccom::expand_seed_to_Fr_vec, smallvole::{VOLE, TestMOLE}, Fr, FrRepr, FrVec, DotProduct, subspacevole::RAAACode, FrMatrix, circom::{witness::wtns_from_reader, r1cs::R1CSFile}, actors::actors::Prover, zkp::R1CSWithMetadata, SparseVec};
use std::{fs::File, io::BufReader};

lazy_static! {
    pub static ref WITNESS: FrVec = {
        let wtns_file = File::open("src/circom/examples/witness.wtns").unwrap();
        let mut wtns_reader = BufReader::new(wtns_file);
        wtns_from_reader(wtns_reader).unwrap()
    };
    pub static ref CIRCUIT: R1CSWithMetadata = {
        let r1cs_file = File::open("src/circom/examples/test.r1cs").unwrap();
        let mut r1cs_reader = BufReader::new(r1cs_file);
        R1CSFile::from_reader(r1cs_reader).unwrap().to_crate_format()
    };
}
fn load_and_prove() {
    let mut prover = Prover::from_witness_and_circuit_unpadded(WITNESS.clone(), CIRCUIT.clone());
    let vole_comm = prover.mkvole().unwrap();
    let proof = prover.prove().unwrap();
}

/// Subroutine in matrix multiplication
fn matrix_row_col_dot<const N: usize>(a: &SMatrix<Fr, N, 1>, b: &SMatrix<Fr, 1, N>) -> Fr {
    a.dot(&b.transpose())
}
fn naive_dot(a: &Vec<Fr>, b: &Vec<Fr>) -> Fr {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum::<Fr>()
}
fn u64_dot(a: &[u64], b: &[u64]) -> u64 {
    a.iter().zip(b.iter()).map(|(a, b)| *a * *b).sum::<u64>()
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("slow");
    // let seed = || vec![0u8; 32];
    let seed = [0u8; 32];
    let seed0 = [9u8; 32];
    let seed1 = [2u8; 32];
    let x = SMatrix::<Fr, 1024, 1>::from_fn(|row, col| {
        Fr::random(&mut ThreadRng::default())
    });
    let y = SMatrix::<Fr, 1, 1024>::from_fn(|row, col| {
        Fr::random(&mut ThreadRng::default())
    });
    let nx = (0..1024).map(|_|Fr::random(&mut ThreadRng::default())).collect::<Vec<_>>();
    let ny = (0..1024).map(|_|Fr::random(&mut ThreadRng::default())).collect::<Vec<_>>();

    let ux: Vec<u64> = (0..1024).map(|_|ThreadRng::default().gen()).collect();
    let uy: Vec<u64> = (0..1024).map(|_|ThreadRng::default().gen()).collect();

    let fvx = FrVec::random(640);
    let fvy = FrVec::random(640);

    let a_ = Fr::random(&mut ThreadRng::default());
    let b_ = Fr::random(&mut ThreadRng::default());

    let fr_matrix = FrMatrix(
        (0..1024).map(|_|
            FrVec((0..1024).map(|_|Fr::random(&mut ThreadRng::default())).collect::<Vec<_>>())
        ).collect()
    );

    let non_sparse_vec_1 = FrVec::random(640);
    let non_sparse_vec_2 = FrVec::random(640);
    let sparse_vec = SparseVec(vec![
        (5, Fr::random(&mut ThreadRng::default())), 
        (21, Fr::random(&mut ThreadRng::default())), 
        (615, Fr::random(&mut ThreadRng::default())), 
    ]);
    let mut repr = [0u8; 32];
    ThreadRng::default().fill(&mut repr);
    group.sample_size(10);

    // group.bench_function("512x1024 MOLE", |b|b.iter(|| TestMOLE::init(
    //     black_box([123u8; 32]), 
    //     black_box(512), 
    //     black_box(1024)
    // )));
    
    // group.bench_function("Tc-1 times 1024-bit vector", |b|b.iter(|| RAAACode::repeat_extended_inverse(black_box(&FrVec(nx.clone())), 2)));
    // group.bench_function("add", |b|b.iter(|| black_box(a_) + black_box(b_)));
    // group.bench_function("mul", |b|b.iter(|| black_box(a_)* black_box(b_)));
    // group.bench_function("Creating a Fr from a repr", |b|b.iter(||Fr::from_repr(black_box(FrRepr(repr)))));
    // group.bench_function("Multiplying two 512x512 matrices", |b|b.iter(||matmul::<512>()));
    // group.bench_function("1024 length nalgebra dot", |b|b.iter(||matrix_row_col_dot(black_box(&x), black_box(&y))));
    // group.bench_function("1024 length simple dot", |b|b.iter(||naive_dot(&nx, &ny)));
    // group.bench_function("1024 length u64 dot", |b|b.iter(||u64_dot(&ux, &uy)));
    // group.bench_function("640 x 640 dot", |b|b.iter(||fvx.dot(&fvy)));
    // group.bench_function("Constructing 288 x 256 Reed Solomon Generator matrix", |b|b.iter(move||ReedSolomonCode::construct_generator::<64, 16>()));
    // group.bench_function("Constructing 1152 x 1024 Reed Solomon Generator matrix quickly", |b|b.iter(move||ReedSolomonCode::construct_generator_quickly::<1152, 1024>()));
    // group.bench_function("Constructing 64 x 64 systematic Reed Solomon Generator matrix", |b|b.iter(move||ReedSolomonCode::construct_systematic_generator::<64, 64>())); // 256 x 256 takes 1.4s
    // group.bench_function("Constructing 288 x 256 Reed Solomon Generator matrix", |b|b.iter(move||ReedSolomonCode::construct_tc_inverse::<288>()));

    // group.bench_function("expand to 4 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed), 4)));
    // group.bench_function("expand to 2^20 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed), 1048576)));
    // group.bench_function("1024 x 1024 transpose", |b|b.iter(||black_box(fr_matrix.clone()).transpose()));
    // group.bench_function("smallvole prover 1024 elements", |b|b.iter(move || VOLE::prover_outputs(black_box(&seed0), black_box(&seed1), 1024)));
    // // group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed()), 10)));
    // // c.bench_function("1048576 random Frs", |b| b.iter(|| rand_fr_vec(black_box(20))));
    // group.bench_function("Sparse Vector Dot", |b|b.iter(|| black_box(&non_sparse_vec_1).sparse_dot(black_box(&sparse_vec))));
    // group.bench_function("Vector Dot", |b|b.iter(|| black_box(&non_sparse_vec_1).dot(black_box(&non_sparse_vec_2))));
    group.bench_function("Load R1CS, Witness, and Create the VOLE in the Head Quicksilver proof", |b|b.iter(load_and_prove));



}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);