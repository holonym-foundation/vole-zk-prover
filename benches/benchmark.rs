use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::{Field, PrimeField};
use lazy_static::lazy_static;
use nalgebra::SMatrix;
use num_bigint::BigUint;
use rand::{rngs::ThreadRng, Rng};
use volonym::{vecccom::expand_seed_to_field_vec, smallvole::{VOLE, TestMOLE}, Fr, FrRepr, FVec, DotProduct, subspacevole::RAAACode, FMatrix, circom::{witness::wtns_from_reader, r1cs::R1CSFile}, actors::actors::Prover, zkp::R1CSWithMetadata, SparseVec, codeparams::{n_choose_k_square_matrix, calc_multi_transition_prob_matrix}};
use std::{fs::File, io::BufReader, str::FromStr};
// use num_modular::{ModularCoreOps};


lazy_static! {
    pub static ref WITNESS: FVec<Fr> = {
        let wtns_file = File::open("src/circom/examples/witness.wtns").unwrap();
        let mut wtns_reader = BufReader::new(wtns_file);
        wtns_from_reader(wtns_reader).unwrap()
    };
    pub static ref CIRCUIT: R1CSWithMetadata<Fr> = {
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

    // let fvx = FrVec::random(640);
    // let fvy = FrVec::random(640);

    let a_fr = Fr::random(&mut ThreadRng::default());
    let b_fr = Fr::random(&mut ThreadRng::default());

    // let a_f4096 = F4096::random(&mut ThreadRng::default());
    // let b_f4096 = F4096::random(&mut ThreadRng::default());

    let mut buf_a = [0u8; 512]; let mut buf_b = [0u8; 512]; 
    ThreadRng::default().fill(&mut buf_a); ThreadRng::default().fill(&mut buf_b); 
    let a_bn = &BigUint::from_bytes_be(&buf_a); let b_bn = &BigUint::from_bytes_be(&buf_b); 
    // 4096-bit prime:
    let m_bn = &BigUint::from_str("634855612763766422949715574974655080572035763098431131173542330951942241520024700573868828181875621819195714866311437400816437746351211821353414491923219000912234985725095974892889028816064244543361245372124850030069611467228117440723680827209043807052645248502562196494398621359970816592705372991289818789085381574821075762441634062960973119681297273928213955393528432295358851956847833931230825840619981691776819186538578209987035463631999172135178899211162484025988921269126198272856159210966145375942240775932278101156163916577767371288428479079115032264646619868163409215860522109374862139439813952153138830598467505836738004996067393286053093730557796288152122588452120725539979401515493820670371193659258133728269206579946230739134712941540113848399816418778005674468871822082743365729620057268465360056920036324997856481876137980138195559616047932205745558953520273772852725319481463187321548139333439002632998147862748563966354503022564820365369708839203425236575958538737699275909049963146375090216548591476247358626989468986537784638657350541405932010938179366613394571228824057256373154131296157692664851840428966735839389686098207300228898576134858803950752431317775566875031174212555605566723709772142593674181603978861");
    
    // let a_modint = MontgomeryInt::new(a_bn, &m_bn);
    // let x = a_bn.mulm(b_bn, m_bn);

    // let fr_matrix = FrMatrix(
    //     (0..1024).map(|_|
    //         FrVec((0..1024).map(|_|Fr::random(&mut ThreadRng::default())).collect::<Vec<_>>())
    //     ).collect()
    // );

    // let non_sparse_vec_1 = FrVec::random(640);
    // let non_sparse_vec_2 = FrVec::random(640);
    // let sparse_vec = SparseVec(vec![
    //     (5, Fr::random(&mut ThreadRng::default())), 
    //     (21, Fr::random(&mut ThreadRng::default())), 
    //     (615, Fr::random(&mut ThreadRng::default())), 
    // ]);
    // let mut repr = [0u8; 32];
    // ThreadRng::default().fill(&mut repr);
    group.sample_size(10);

    // group.bench_function("512x1024 MOLE", |b|b.iter(|| TestMOLE::init(
    //     black_box([123u8; 32]), 
    //     black_box(512), 
    //     black_box(1024)
    // )));
    
    // group.bench_function("Tc-1 times 1024-bit vector", |b|b.iter(|| RAAACode::repeat_extended_inverse(black_box(&FrVec(nx.clone())), 2)));
    // group.bench_function("add", |b|b.iter(|| black_box(a_) + black_box(b_)));
    group.bench_function("Fr mul", |b|b.iter(|| black_box(a_fr)* black_box(b_fr)));
    // group.bench_function("bn 4096 mul", |b|b.iter(|| black_box(a_bn).mulm(black_box(b_bn), black_box(m_bn))));
    // group.bench_function("bn 4096 add", |b|b.iter(|| black_box(a_bn) * black_box(b_bn)));

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
    // group.bench_function("Calculate 1024x1024 binomial coefficient (n choose k) mat", |b|b.iter(||n_choose_k_square_matrix(black_box(1024))));
    // group.bench_function("Calculate 1024x1024 3x probability matrix", |b|b.iter(||calc_multi_transition_prob_matrix(64, 3)));
    // group.bench_function("Calculate circuit id", |b|b.iter(||{ let c = black_box(&*CIRCUIT); c.circuit_id().unwrap() }));
    

}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);