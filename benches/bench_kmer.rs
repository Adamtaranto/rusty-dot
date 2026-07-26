//! Benchmark: per-sequence index primitives.
//!
//! Contrasts the new rolling-hash index build against the legacy building
//! blocks (exact k-mer set + FM-index construction) on the same sequence.

mod common;

use _rusty_dot::kmer::{build_kmer_set, sequence_to_index_text, FmIdx};
use _rusty_dot::kmer_hash::KmerIndex;
use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, Criterion};
use common::random_dna;

fn bench_primitives(c: &mut Criterion) {
    let seq = random_dna(20_000, 3);
    let seq_str = std::str::from_utf8(&seq).unwrap();
    let k = 21usize;

    let mut group = c.benchmark_group("primitives");
    group.bench_function("build_kmer_set", |b| {
        b.iter(|| build_kmer_set(black_box(seq_str), k).unwrap())
    });
    group.bench_function("fm_index_build", |b| {
        b.iter(|| {
            let text = sequence_to_index_text(black_box(seq_str));
            FmIdx::new(text).unwrap()
        })
    });
    group.bench_function("kmer_index_build", |b| {
        b.iter(|| KmerIndex::build(black_box(&seq), k))
    });
    group.finish();
}

criterion_group!(benches, bench_primitives);
criterion_main!(benches);
