"""Tests for the browser app's plot-layout state (W2: interactive plot).

Covers the ``'colinearity_ref'`` ordering mode, the explicit-order
computation (target axis fixed, queries reordered/reoriented against it),
the panel (row, col) -> contig-name mapping, and the reordered-FASTA export
for both the k-mer and PAF-with-sequences paths.
"""

from pathlib import Path
import random
import sys

import matplotlib
import pytest

matplotlib.use('Agg')

APP_DIR = Path(__file__).resolve().parent.parent / 'app'
sys.path.insert(0, str(APP_DIR))

from core.align import paf_alignment_from_text  # noqa: E402
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache  # noqa: E402
from core.export import reordered_fasta_text  # noqa: E402
from core.fasta import parse_fasta_bytes  # noqa: E402
from core.panels import panel_pair, resolve_orders  # noqa: E402
from core.state import ORDER_CHOICES, PlotConfig  # noqa: E402

# ---------------------------------------------------------------- fixtures


def _rand_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(n))


@pytest.fixture(scope='module')
def kmer_setup():
    """CrossIndex with a known collinear layout.

    Target t1 = A + C (800 bp), t2 = B (400 bp).  Queries: qA matches the
    first half of t1, qRevA reverse-matches the second half of t1, qB
    matches t2.  Query upload order is deliberately scrambled.
    """
    from rusty_dot.paf_io import reverse_complement

    a = _rand_seq(400, seed=1)
    b = _rand_seq(400, seed=2)
    c = _rand_seq(400, seed=3)
    q_fa = parse_fasta_bytes(
        f'>qB\n{b}\n>qRevA\n{reverse_complement(c)}\n>qA\n{a}\n'.encode()
    )
    t_fa = parse_fasta_bytes(f'>t1\n{a + c}\n>t2\n{b}\n'.encode())
    idx = SessionCache().kmer_index(11, q_fa, t_fa)
    return idx, q_fa, t_fa


# ---------------------------------------------------------------- state


def test_order_choices_has_colinearity_ref():
    assert 'colinearity_ref' in ORDER_CHOICES
    assert ORDER_CHOICES['colinearity_ref'] == 'Colinearity vs fixed reference'


def test_plot_kwargs_colinearity_ref_maps_to_no_plot_reorder():
    # colinearity_ref ordering is computed app-side and passed as explicit
    # name lists, never as a DotPlotter.plot contig_order strategy.
    kwargs = PlotConfig(contig_order='colinearity_ref').plot_kwargs()
    assert kwargs['contig_order'] is None
    # Plot-level strategies still pass through unchanged.
    assert PlotConfig(contig_order='length').plot_kwargs()['contig_order'] == 'length'
    assert (
        PlotConfig(contig_order='colinearity').plot_kwargs()['contig_order']
        == 'colinearity'
    )


# ---------------------------------------------------------------- panels


def test_panel_pair_maps_row_major():
    q = ['q1', 'q2']
    t = ['t1', 't2', 't3']
    assert panel_pair(q, t, 0, 0) == ('q1', 't1')
    assert panel_pair(q, t, 1, 2) == ('q2', 't3')


def test_panel_pair_rejects_out_of_grid():
    with pytest.raises(IndexError):
        panel_pair(['q1'], ['t1'], 1, 0)
    with pytest.raises(IndexError):
        panel_pair(['q1'], ['t1'], 0, 1)
    with pytest.raises(IndexError):
        panel_pair(['q1'], ['t1'], -1, 0)


def test_resolve_orders_input_and_length():
    q_in, t_in = ['qs', 'ql'], ['tl', 'ts']
    lengths = {'qs': 10, 'ql': 100, 'tl': 90, 'ts': 5}
    q, t, rev = resolve_orders('input', [], q_in, t_in, lengths)
    assert (q, t, rev) == (q_in, t_in, set())
    assert q is not q_in  # copies, never aliases
    q, t, rev = resolve_orders('length', [], q_in, t_in, lengths)
    assert q == ['ql', 'qs']
    assert t == ['tl', 'ts']
    assert rev == set()


def test_resolve_orders_unknown_mode():
    with pytest.raises(ValueError, match='Unknown contig-order mode'):
        resolve_orders('bogus', [], [], [], {})


# ------------------------------------------- colinearity vs fixed reference


def test_colinearity_ref_keeps_target_fixed(kmer_setup):
    idx, q_fa, t_fa = kmer_setup
    records = idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
    assert records, 'fixture must produce k-mer matches'
    q_order, t_order, rev = resolve_orders(
        'colinearity_ref', records, q_fa.names, t_fa.names, {}
    )
    assert t_order == ['t1', 't2']  # reference axis untouched
    assert q_order == ['qA', 'qRevA', 'qB']  # gravity along t1 then t2
    assert rev == {'qRevA'}


def test_colinearity_ref_respects_arbitrary_target_order(kmer_setup):
    idx, q_fa, _t_fa = kmer_setup
    records = idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
    # Reference given in a deliberately 'suboptimal' order: still fixed.
    q_order, t_order, _rev = resolve_orders(
        'colinearity_ref', records, q_fa.names, ['t2', 't1'], {}
    )
    assert t_order == ['t2', 't1']
    assert q_order == ['qB', 'qA', 'qRevA']


def test_colinearity_reorders_both_axes(kmer_setup):
    idx, q_fa, t_fa = kmer_setup
    records = idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
    q_order, t_order, rev = resolve_orders(
        'colinearity', records, q_fa.names, t_fa.names, {}
    )
    assert sorted(q_order) == sorted(q_fa.names)
    assert sorted(t_order) == sorted(t_fa.names)
    assert rev == {'qRevA'}


def test_resolve_orders_does_not_mutate_index_or_inputs(kmer_setup):
    idx, q_fa, t_fa = kmer_setup
    before_groups = {g: idx.sequence_names(group=g) for g in idx.group_names}
    q_in, t_in = list(q_fa.names), list(t_fa.names)
    records = idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
    resolve_orders('colinearity_ref', records, q_in, t_in, {})
    resolve_orders('colinearity', records, q_in, t_in, {})
    assert q_in == q_fa.names
    assert t_in == t_fa.names
    assert {g: idx.sequence_names(group=g) for g in idx.group_names} == before_groups


# ---------------------------------------------------------------- export


def test_fasta_export_kmer_path(kmer_setup, tmp_path):
    from rusty_dot.paf_io import reverse_complement

    idx, q_fa, t_fa = kmer_setup
    records = idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
    order, _t, rev = resolve_orders(
        'colinearity_ref', records, q_fa.names, t_fa.names, {}
    )
    out = tmp_path / 'reordered.fasta'
    idx.write_fasta(out, group=QUERY_GROUP, order=order, reverse=rev)
    text = out.read_text()
    headers = [line for line in text.splitlines() if line.startswith('>')]
    assert headers == ['>qA', '>qRevA reverse_complement', '>qB']
    # The flipped contig now reads forward along t1's second half.
    reparsed = parse_fasta_bytes(text.encode())
    seq_map = dict(q_fa.records)
    assert dict(reparsed.records)['qRevA'] == reverse_complement(seq_map['qRevA'])
    assert dict(reparsed.records)['qA'] == seq_map['qA']


def test_fasta_export_paf_with_sequences_path():
    from rusty_dot.paf_io import reverse_complement

    paf = (
        'q1\t100\t0\t50\t+\tt1\t200\t100\t150\t50\t50\t60\n'
        'q2\t100\t0\t50\t-\tt1\t200\t0\t50\t50\t50\t60\n'
    )
    aln = paf_alignment_from_text(paf)
    q_fa = parse_fasta_bytes(b'>q1\nACGTACGTAC\n>q2\nTTGGCCAATT\n>q3\nGGGG\n')
    lengths = {n: aln.get_sequence_length(n) for n in aln.sequence_names()}
    order, t_order, rev = resolve_orders(
        'colinearity_ref', aln.records, aln.query_names, aln.target_names, lengths
    )
    assert t_order == ['t1']
    assert order == ['q2', 'q1']  # q2's match sits earlier on the fixed axis
    assert rev == {'q2'}  # single '-'-strand match => reverse-oriented
    text = reordered_fasta_text(q_fa.records, order, rev)
    headers = [line for line in text.splitlines() if line.startswith('>')]
    # q3 has no alignments: appended after the ordered contigs.
    assert headers == ['>q2 reverse_complement', '>q1', '>q3']
    reparsed = parse_fasta_bytes(text.encode())
    assert dict(reparsed.records)['q2'] == reverse_complement('TTGGCCAATT')
    assert dict(reparsed.records)['q1'] == 'ACGTACGTAC'
    assert dict(reparsed.records)['q3'] == 'GGGG'


def test_reordered_fasta_text_line_wrapping():
    text = reordered_fasta_text([('c1', 'A' * 130)], ['c1'], set(), line_width=60)
    lines = text.splitlines()
    assert lines[0] == '>c1'
    assert [len(line) for line in lines[1:]] == [60, 60, 10]
    assert reordered_fasta_text([('c1', 'ACGT')], ['c1'], set(), line_width=0) == (
        '>c1\nACGT\n'
    )


# ------------------------------------------------------- report injection


def test_inject_panel_bridge_places_script_before_body_close():
    pytest.importorskip('shiny')  # app.py imports shiny at module level
    import app as app_module

    html = '<html><body><svg></svg></body></html>'
    out = app_module.inject_panel_bridge(html)
    assert 'rd-panel-dblclick' in out
    assert out.index('<script>') < out.index('</body>')
    assert out.endswith('</body></html>')
    # No closing body tag: script is appended.
    assert app_module.inject_panel_bridge('<svg/>').startswith('<svg/>')


def test_generated_report_contains_bridge_and_panels(kmer_setup, tmp_path):
    """End-to-end: to_html + injection yields panels and a valid bridge."""
    pytest.importorskip('shiny')  # app.py imports shiny at module level
    import app as app_module
    from rusty_dot import DotPlotter
    from rusty_dot.paf_io import PafAlignment

    idx, q_fa, t_fa = kmer_setup
    paf = PafAlignment(idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP))
    plotter = DotPlotter(idx, paf_alignment=paf)
    out = tmp_path / 'report.html'
    q_int = [idx.make_internal_name(QUERY_GROUP, n) for n in q_fa.names]
    t_int = [idx.make_internal_name(TARGET_GROUP, n) for n in t_fa.names]
    fig = plotter.to_html(out, query_names=q_int, target_names=t_int)
    import matplotlib.pyplot as plt

    plt.close(fig)
    html = app_module.inject_panel_bridge(out.read_text())
    assert 'id="rd-panel-0-0"' in html
    assert f'id="rd-panel-{len(q_int) - 1}-{len(t_int) - 1}"' in html
    assert html.count('rd-panel-dblclick') == 1
    assert 'window.parent.postMessage' in html
    # Injected regex survived Python string escaping intact.
    assert r'/^rd-panel-(\d+)-(\d+)$/' in html


# ----------------------------------------------------- UX round 2 helpers


def test_strip_report_header_hides_report_hint_bar():
    pytest.importorskip('shiny')  # app.py imports shiny at module level
    import app as app_module

    html = '<html><head><title>t</title></head><body><svg/></body></html>'
    out = app_module.strip_report_header(html)
    assert '#rd-header{display:none}' in out
    assert out.index('#rd-header') < out.index('</head>')
    # No head tag: style is prepended, document body untouched.
    naked = app_module.strip_report_header('<svg/>')
    assert naked.startswith('<style>')
    assert naked.endswith('<svg/>')


def test_method_choices_exclude_paf_input_mode():
    pytest.importorskip('shiny')
    import app as app_module

    choices = app_module._method_choices()
    assert 'paf' not in choices
    assert {'kmer', 'minimap2', 'lastz', 'nucmer'} <= set(choices)


def test_plot_config_identity_kwargs():
    from core.state import PlotConfig

    kwargs = PlotConfig().plot_kwargs()
    assert kwargs['color_by_identity'] is False
    assert kwargs['identity_palette'] == 'viridis'
    kwargs = PlotConfig(color_by_identity=True, identity_palette='plasma')
    assert kwargs.plot_kwargs()['color_by_identity'] is True
    assert kwargs.plot_kwargs()['identity_palette'] == 'plasma'


def test_identity_colored_plot_from_mixed_identity_records():
    """color_by_identity renders per-record identity colours for PAF results."""
    import matplotlib.pyplot as plt

    from rusty_dot import DotPlotter
    from rusty_dot.paf_io import PafAlignment, PafRecord

    def rec(q, t, matches, block):
        return PafRecord(
            query_name=q,
            query_len=5000,
            query_start=0,
            query_end=block,
            strand='+',
            target_name=t,
            target_len=5000,
            target_start=0,
            target_end=block,
            residue_matches=matches,
            alignment_block_len=block,
            mapping_quality=255,
        )

    aln = PafAlignment([rec('q1', 't1', 4000, 4000), rec('q2', 't1', 2000, 4000)])
    fig = DotPlotter(aln).plot(
        query_names=['q1', 'q2'],
        target_names=['t1'],
        color_by_identity=True,
        identity_palette='plasma',
    )
    assert fig.axes
    plt.close(fig)


def test_self_align_kmer_index_and_plot():
    """Self-alignment reuses one FastaInput for both groups."""
    from core.cache import SessionCache
    from core.fasta import parse_fasta_bytes
    import matplotlib.pyplot as plt

    from rusty_dot import DotPlotter
    from rusty_dot.paf_io import PafAlignment

    seq = 'ACGTTGCAAGGCTTAACCGGTTAACGGCCAATT' * 8
    q = parse_fasta_bytes(f'>c1\n{seq}\n>c2\n{seq[::-1]}\n'.encode())
    cache = SessionCache()
    idx = cache.kmer_index(11, q, q)
    # Same digest on both sides: the cache key is still well-formed and a
    # second call hits the cache.
    assert cache.kmer_index(11, q, q) is idx
    paf = PafAlignment(idx.get_records_for_pair(QUERY_GROUP, TARGET_GROUP))
    # The self-vs-self diagonal must be present (c1 x c1 full-length match).
    diag = [r for r in paf.records if r.query_name == r.target_name]
    assert diag
    fig = DotPlotter(idx, paf_alignment=paf).plot(
        query_group=QUERY_GROUP, target_group=TARGET_GROUP
    )
    assert fig.axes
    plt.close(fig)
