use criterion::{black_box, criterion_group, criterion_main, Criterion};
use volonym::{vecccom::{expand_seed_to_Fr_vec, expand_seed_to_Fr_vec_faster}, smallvole::VOLE};
// use volonym::rand_fr_vec;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("slow");
    // let seed = || vec![0u8; 32];
    let seed = [0u8; 32];
    let seed0 = [9u8; 32];
    let seed1 = [2u8; 32];
    group.sample_size(10);
    group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(&seed), 1024)));
    group.bench_function("expand to 2^10 Frs; optimized", |b|b.iter(move ||expand_seed_to_Fr_vec_faster(black_box(&seed), 1024)));

    // group.bench_function("smallvole prover 1024 elements", |b|b.iter(move || VOLE::prover_outputs(black_box(&seed0), black_box(&seed1), 1024)));
    // group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed()), 10)));
    // c.bench_function("1048576 random Frs", |b| b.iter(|| rand_fr_vec(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);