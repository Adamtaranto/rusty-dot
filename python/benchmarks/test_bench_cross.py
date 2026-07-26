"""CodSpeed benchmarks for the CrossIndex build + match path.

Run with ``pytest python/benchmarks --codspeed``.  These measure the two
performance-critical steps for a cross-genome comparison: building the
per-sequence rolling-hash index for every contig (parallelised across contigs)
and computing k-mer matches for every cross-group pair.
"""

from __future__ import annotations

from rusty_dot.paf_io import CrossIndex

from ._synth import multi_contig_group

# Kept modest so CodSpeed's instruction-count simulation stays quick while still
# exercising many-contig parallel builds and an N x N pairwise match grid.
_K = 15
_N_CONTIGS = 6
_CONTIG_LEN = 3000


def _build_cross() -> CrossIndex:
    """Build a CrossIndex with two multi-contig groups loaded."""
    cross = CrossIndex(k=_K)
    for name, seq in multi_contig_group(_N_CONTIGS, _CONTIG_LEN, seed=1):
        cross.add_sequence(name, seq, group='A')
    for name, seq in multi_contig_group(_N_CONTIGS, _CONTIG_LEN, seed=1000):
        cross.add_sequence(name, seq, group='B')
    return cross


def test_bench_build_cross_index(benchmark):
    """Benchmark building the per-contig k-mer index for both groups."""
    benchmark(_build_cross)


def test_bench_compute_matches(benchmark):
    """Benchmark the end-to-end build + all-pairs match computation."""

    def run():
        cross = _build_cross()
        cross.compute_matches(query_group='A', target_group='B')
        return cross

    benchmark(run)
