use criterion::{black_box, criterion_group, criterion_main, Criterion};
use volonym::expand_seed_to_Fr_vec;
// use volonym::rand_fr_vec;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("slow");
    let seed = || vec![0u8; 32];
    group.sample_size(10);
    group.bench_function("expand to 2^10 Frs", |b|b.iter(move ||expand_seed_to_Fr_vec(black_box(seed()), 10)));
    // c.bench_function("1048576 random Frs", |b| b.iter(|| rand_fr_vec(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);