use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::{Field, PrimeField};
use nalgebra::SMatrix;
use rand::{rngs::ThreadRng, Rng};
use volonym::{vecccom::expand_seed_to_Fr_vec, smallvole::VOLE, subspacevole::ReedSolomonCode, Fr, FrRepr};
// use volonym::rand_fr_vec;

// fn matmul<const N: usize>() {
//     let x1 = SMatrix::<Fr, N, N>::from_fn(|row, col| {
//         Fr::random(&mut ThreadRng::default())
//     });

//     let x1 = SMatrix::<Fr, N, N>::from_fn(|row, col| {
//         Fr::random(&mut ThreadRng::default())
//     });

//     let _ = x1 * x1;
// }

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

    let mut repr = [0u8; 32];
    ThreadRng::default().fill(&mut repr);
    group.sample_size(100);

    // group.bench_function("Creating a Fr from a repr", |b|b.iter(||Fr::from_repr(black_box(FrRepr(repr)))));
    // group.bench_function("Multiplying two 512x512 matrices", |b|b.iter(||matmul::<512>()));
    // group.bench_function("1024 length nalgebra dot", |b|b.iter(||matrix_row_col_dot(black_box(&x), black_box(&y))));
    // group.bench_function("1024 length simple dot", |b|b.iter(||naive_dot(&nx, &ny)));
    // group.bench_function("1024 length u64 dot", |b|b.iter(||u64_dot(&ux, &uy)));
    // group.bench_function("Constructing 288 x 256 Reed Solomon Generator matrix", |b|b.iter(move||ReedSolomonCode::construct_generator::<64, 16>()));
    group.bench_function("Constructing 1152 x 1024 Reed Solomon Generator matrix quickly", |b|b.iter(move||ReedSolomonCode::construct_generator_quickly::<1152, 1024>()));
    // group.bench_function("Constructing 288 x 256 Reed Solomon Generator matrix", |b|b.iter(move||ReedSolomonCode::construct_tc_inverse::<288>()));

    // group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed), 1048576)));
    // // group.bench_function("expand to 2^10 Frs; optimized", |b|b.iter(move ||expand_seed_to_Fr_vec_faster(black_box(seed), 1024)));

    // group.bench_function("smallvole prover 1024 elements", |b|b.iter(move || VOLE::prover_outputs(black_box(&seed0), black_box(&seed1), 1024)));
    // // group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed()), 10)));
    // // c.bench_function("1048576 random Frs", |b| b.iter(|| rand_fr_vec(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);