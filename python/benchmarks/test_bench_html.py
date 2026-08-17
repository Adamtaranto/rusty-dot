"""CodSpeed benchmarks for the interactive HTML report pipeline.

Covers the two costs specific to HTML output on top of a normal plot:
payload serialisation (:func:`build_panel_payload`, including sequence
slicing) and the full :meth:`DotPlotter.to_html` render.
"""

import matplotlib

matplotlib.use('Agg')

from _synth import homologous_pair
import matplotlib.pyplot as plt
import pytest

from rusty_dot import SequenceIndex
from rusty_dot._html.serialize import build_panel_payload
from rusty_dot.dotplot import DotPlotter

# Small enough to keep the figure render tractable under instrumentation,
# large enough to produce hundreds of match segments per panel.
_LENGTH = 20_000
_K = 15


@pytest.fixture(scope='module')
def html_bench_index() -> SequenceIndex:
    """Index over a homologous pair with planted inversions."""
    query, target = homologous_pair(_LENGTH, seed=42, divergence=0.02)
    idx = SequenceIndex(k=_K)
    idx.add_sequence('query', query)
    idx.add_sequence('target', target)
    return idx


@pytest.fixture(scope='module')
def html_capture(html_bench_index: SequenceIndex) -> dict:
    """Raw per-panel capture mirroring what plot() assembles for HTML output."""
    names = sorted(html_bench_index.sequence_names())
    panels = {}
    for row, q_name in enumerate(names):
        for col, t_name in enumerate(names):
            fwd = []
            rev = []
            for qs, qe, ts, te, strand in html_bench_index.compare_sequences_stranded(
                q_name, t_name, True
            ):
                (rev if strand == '-' else fwd).append([qs, qe, ts, te])
            panels[f'rd-panel-{row}-{col}'] = {
                'query': q_name,
                'target': t_name,
                'query_id': q_name,
                'qlen': html_bench_index.get_sequence_length(q_name),
                'tlen': html_bench_index.get_sequence_length(t_name),
                'segments': {'fwd': fwd, 'rev': rev, 'identity': []},
            }
    return {'panels': panels}


def test_bench_build_panel_payload(benchmark, html_capture, html_bench_index):
    """Serialise the captured panels, including sequence embedding."""
    result = benchmark(
        build_panel_payload,
        html_capture,
        get_sequence=html_bench_index.get_sequence,
    )
    assert result['has_sequences'] is True
    assert len(result['panels']) == 4


def test_bench_to_html(benchmark, html_bench_index, tmp_path):
    """Full interactive report render (plot + SVG + payload + assembly)."""
    out = tmp_path / 'bench_report.html'

    def run():
        fig = DotPlotter(html_bench_index).to_html(out)
        plt.close(fig)

    benchmark(run)
    assert out.exists()
