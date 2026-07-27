"""CodSpeed benchmarks for the collinearity contig-reordering path.

Run with ``pytest python/benchmarks --codspeed``.  These cover the d-genies
gravity sort (squared weighting + best-matching-chromosome argmax) and the
reverse-orientation detection added alongside it, on both the ``CrossIndex``
(Rust ordering + Python reversal) and the pure-Python ``compute_gravity_contigs``
paths.
"""

from __future__ import annotations

from _synth import multi_contig_group, revcomp

from rusty_dot.paf_io import CrossIndex, PafRecord, compute_gravity_contigs

# Modest sizes so CodSpeed's instruction-count simulation stays quick while
# still exercising a multi-contig, multi-target layout with inversions.
_K = 15
_N_CONTIGS = 5
_CONTIG_LEN = 2000


def _cross_with_matches() -> CrossIndex:
    """Build a CrossIndex whose group B mirrors A, with half the contigs RC'd."""
    contigs = multi_contig_group(_N_CONTIGS, _CONTIG_LEN, seed=7)
    cross = CrossIndex(k=_K)
    for name, seq in contigs:
        cross.add_sequence(name, seq, group='A')
    # Group B reuses A's sequences (so they align), reverse-complementing every
    # other one to exercise reverse-orientation detection.
    for i, (_name, seq) in enumerate(contigs):
        cross.add_sequence(f'ref_{i}', revcomp(seq) if i % 2 else seq, group='B')
    cross.compute_matches('A', 'B')
    return cross


def test_bench_reorder_for_colinearity(benchmark):
    """Benchmark build + match + gravity reorder + reversal detection."""

    def run():
        cross = _cross_with_matches()
        cross.reorder_for_colinearity('A', 'B')
        return cross

    benchmark(run)


def _synthetic_records() -> list[PafRecord]:
    """Fabricate a multi-target PAF record set (no k-mer engine involved)."""
    records: list[PafRecord] = []
    n_queries, n_targets, t_len = 40, 4, 5000
    for q in range(n_queries):
        best_t = q % n_targets
        # A few collinear blocks against the best target, plus one small hit to
        # a neighbouring target to exercise the argmax step.
        for b in range(4):
            qs = b * 200
            ts = (q // n_targets) * 400 + b * 200
            records.append(
                PafRecord.from_line(
                    f'q{q}\t2000\t{qs}\t{qs + 180}\t+\t'
                    f't{best_t}\t{t_len}\t{ts}\t{ts + 180}\t170\t180\t255'
                )
            )
        other_t = (best_t + 1) % n_targets
        records.append(
            PafRecord.from_line(
                f'q{q}\t2000\t900\t920\t+\tt{other_t}\t{t_len}\t50\t70\t18\t20\t255'
            )
        )
    return records


_RECORDS = _synthetic_records()
_QUERIES = [f'q{q}' for q in range(40)]
_TARGETS = [f't{t}' for t in range(4)]


def test_bench_compute_gravity_contigs(benchmark):
    """Benchmark the pure-Python gravity sort + reversal detection."""
    benchmark(lambda: compute_gravity_contigs(_RECORDS, _QUERIES, _TARGETS))
