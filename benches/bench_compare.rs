//! Benchmark: forward and reverse-strand k-mer matching between two sequences.
//!
//! This is the headline path — computing shared-k-mer coordinates for a
//! homologous pair with planted inversions.  Index construction is done once
//! outside the measured loop so the benchmark isolates the matching cost.

mod common;

use _rusty_dot::kmer_hash::{shared_fwd_coords, shared_rev_coords, KmerIndex};
use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, Criterion};
use common::homologous_pair;

fn bench_compare(c: &mut Criterion) {
    let k = 21usize;
    let (query, target) = homologous_pair(50_000, 7, 0.02, 3);
    let q_index = KmerIndex::build(&query, k);
    let t_index = KmerIndex::build(&target, k);

    let mut group = c.benchmark_group("compare");
    group.bench_function("forward", |b| {
        b.iter(|| shared_fwd_coords(black_box(&query), &q_index, black_box(&target), &t_index, k))
    });
    group.bench_function("reverse", |b| {
        b.iter(|| shared_rev_coords(black_box(&query), &q_index, black_box(&target), &t_index, k))
    });
    group.finish();
}

criterion_group!(benches, bench_compare);
criterion_main!(benches);
