"""CodSpeed benchmarks for dotplot rendering.

Run with ``pytest python/benchmarks --codspeed``.  Guards the vectorised
``LineCollection`` rendering path in ``DotPlotter._plot_panel``.
"""

from __future__ import annotations

import matplotlib

# Force a non-interactive backend before pyplot is imported anywhere.
matplotlib.use('Agg')

from _synth import homologous_pair

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter

_K = 15
_LEN = 4000


def _build_index() -> SequenceIndex:
    """Build a two-sequence index with forward + reverse-complement homology."""
    query, target = homologous_pair(_LEN, seed=7, n_inversions=3)
    idx = SequenceIndex(k=_K)
    idx.add_sequence('query', query)
    idx.add_sequence('target', target)
    return idx


def test_bench_plot_kmer(benchmark, tmp_path):
    """Benchmark rendering a dotplot from k-mer matches to a file (default auto)."""
    idx = _build_index()
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'bench.png')

    def run():
        fig = plotter.plot(output_path=output)
        # Free the figure so repeated iterations do not accumulate memory.
        import matplotlib.pyplot as plt

        plt.close(fig)

    benchmark(run)


def test_bench_plot_kmer_chained(benchmark, tmp_path):
    """Benchmark rendering with co-linear chaining enabled (chain_gap > 0)."""
    idx = _build_index()
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'bench_chained.svg')

    def run():
        fig = plotter.plot(output_path=output, chain_gap=500)
        import matplotlib.pyplot as plt

        plt.close(fig)

    benchmark(run)
