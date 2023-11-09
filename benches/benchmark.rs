use criterion::{black_box, criterion_group, criterion_main, Criterion};
use volonym::rand_fr_vec;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("1048576 Frs", |b| b.iter(|| rand_fr_vec(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);