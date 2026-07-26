//! Benchmark: building the rolling-hash k-mer index across sizes and k values.

mod common;

use _rusty_dot::kmer_hash::KmerIndex;
use codspeed_criterion_compat::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use common::random_dna;

fn bench_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_index_build");
    for &len in &[10_000usize, 50_000] {
        let seq = random_dna(len, 42);
        for &k in &[15usize, 51, 121] {
            group.bench_with_input(BenchmarkId::new(format!("k{k}"), len), &len, |b, _| {
                b.iter(|| KmerIndex::build(black_box(&seq), black_box(k)))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_build);
criterion_main!(benches);
