"""Tests for the interactive HTML dotplot report (rusty_dot._html)."""

import json
import re

import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pytest

from rusty_dot._html.serialize import build_panel_payload
from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter
from rusty_dot.paf_io import PafAlignment, PafRecord


def _panel_segments(ax):
    """Return ``[(x0, y0, x1, y1), ...]`` for every match segment on ``ax``.

    Local copy of the helper used in test_dotplot.py: match layers are
    LineCollections, so per-segment inspection goes through ``ax.collections``.
    """
    from matplotlib.collections import LineCollection

    segments = []
    for coll in ax.collections:
        if isinstance(coll, LineCollection):
            for seg in coll.get_segments():
                (x0, y0), (x1, y1) = seg[0], seg[-1]
                segments.append((x0, y0, x1, y1))
    return segments


def _read_payload(html_path):
    """Extract and parse the embedded JSON payload from a report file."""
    html = html_path.read_text()
    m = re.search(
        r'<script type="application/json" id="rd-data">(.*?)</script>',
        html,
        re.S,
    )
    assert m is not None, 'payload script block missing'
    return json.loads(m.group(1).replace('<\\/', '</'))


def _segment_count(payload):
    """Total number of serialized segments across all panels and layers."""
    return sum(
        len(panel['segments'][layer])
        for panel in payload['panels'].values()
        for layer in ('fwd', 'rev', 'identity')
    )


@pytest.fixture
def html_index(fasta_file):
    """Sequence index built from the shared FASTA fixture."""
    idx = SequenceIndex(k=4)
    idx.load_fasta(fasta_file)
    return idx


def _make_paf_records():
    """Two-record PAF alignment (one forward, one reverse) for two contigs."""
    return [
        PafRecord(
            query_name='q1',
            query_len=1000,
            query_start=100,
            query_end=400,
            strand='+',
            target_name='t1',
            target_len=900,
            target_start=50,
            target_end=350,
            residue_matches=270,
            alignment_block_len=300,
            mapping_quality=60,
            tags={},
        ),
        PafRecord(
            query_name='q1',
            query_len=1000,
            query_start=500,
            query_end=800,
            strand='-',
            target_name='t1',
            target_len=900,
            target_start=400,
            target_end=700,
            residue_matches=240,
            alignment_block_len=300,
            mapping_quality=60,
            tags={},
        ),
    ]


# ---------------------------------------------------------------------------
# File structure
# ---------------------------------------------------------------------------


def test_to_html_creates_selfcontained_file(html_index, tmp_path):
    out = tmp_path / 'report.html'
    fig = DotPlotter(html_index).to_html(out, title='my report')
    plt.close(fig)

    assert out.exists()
    html = out.read_text()
    assert '<svg' in html
    assert 'id="rd-data"' in html
    assert '<title>my report</title>' in html
    # Assets are inlined, not referenced.
    assert '<link' not in html
    assert 'src=' not in html.split('<svg')[0]
    # No unreplaced template slots remain.
    assert '__RD_' not in html


def test_copy_buttons_start_in_fetch_state(html_index, tmp_path):
    """The copy buttons say what the *next* press does, not what they are for.

    Sequences are fetched lazily from the embedding app, so the first press
    only fetches; the label flips to 'Copy …' once one is cached.
    """
    out = tmp_path / 'report.html'
    plt.close(DotPlotter(html_index).to_html(out))
    html = out.read_text()

    assert re.search(r'id="rd-copy-query"[^>]*>Fetch query seq<', html)
    assert re.search(r'id="rd-copy-target"[^>]*>Fetch target seq<', html)
    # The label is derived from state, never captured as a constant.
    assert 'function labelFor(' in html
    assert 'Press again to copy' not in html


def test_plot_html_suffix_dispatch(html_index, tmp_path):
    """plot(output_path='x.html') routes to the HTML renderer."""
    out = tmp_path / 'grid.html'
    fig = DotPlotter(html_index).plot(output_path=out)
    plt.close(fig)
    assert '<!DOCTYPE html>' in out.read_text()


def test_plot_explicit_html_format(html_index, tmp_path):
    """format='html' forces HTML output regardless of suffix."""
    out = tmp_path / 'grid.dat'
    fig = DotPlotter(html_index).plot(output_path=out, format='html')
    plt.close(fig)
    assert '<!DOCTYPE html>' in out.read_text()


def test_html_suffix_with_other_format_saves_normally(html_index, tmp_path):
    """An explicit non-HTML format wins over a .html suffix."""
    out = tmp_path / 'grid.html'
    fig = DotPlotter(html_index).plot(output_path=out, format='svg')
    plt.close(fig)
    content = out.read_text()
    assert '<!DOCTYPE html>' not in content
    assert '<svg' in content


# ---------------------------------------------------------------------------
# Payload content
# ---------------------------------------------------------------------------


def test_payload_panel_count_and_gids(html_index, tmp_path):
    out = tmp_path / 'report.html'
    fig = DotPlotter(html_index).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    names = html_index.sequence_names()
    assert len(payload['panels']) == len(names) * len(names)

    html = out.read_text()
    for gid in payload['panels']:
        assert f'id="{gid}"' in html
        row, col = gid.rsplit('-', 2)[1:]
        panel = payload['panels'][gid]
        for layer in ('fwd', 'rev'):
            if panel['segments'][layer]:
                assert f'id="rd-matches-{row}-{col}-{layer}"' in html

    # Metadata sanity: names and lengths match the index.
    for panel in payload['panels'].values():
        assert panel['query'] in names
        assert panel['target'] in names
        assert panel['qlen'] == html_index.get_sequence_length(panel['query'])
        assert panel['tlen'] == html_index.get_sequence_length(panel['target'])


def test_payload_segment_counts_match_drawn_collections(html_index, tmp_path):
    out = tmp_path / 'report.html'
    fig = DotPlotter(html_index).to_html(out)

    payload = _read_payload(out)
    names = sorted(html_index.sequence_names())
    n = len(names)
    total_drawn = 0
    for row in range(n):
        for col in range(n):
            ax = fig.axes[row * n + col]
            drawn = len(_panel_segments(ax))
            panel = payload['panels'][f'rd-panel-{row}-{col}']
            serialized = sum(
                len(panel['segments'][layer]) for layer in ('fwd', 'rev', 'identity')
            )
            assert serialized == drawn
            total_drawn += drawn
    assert total_drawn > 0
    plt.close(fig)


def test_payload_svg_path_counts_match_segments(html_index, tmp_path):
    """The nth SVG element under a match group maps to the nth payload row."""
    out = tmp_path / 'report.html'
    fig = DotPlotter(html_index).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    html = out.read_text()
    for gid, panel in payload['panels'].items():
        row, col = gid.rsplit('-', 2)[1:]
        for layer in ('fwd', 'rev', 'identity'):
            n_segs = len(panel['segments'][layer])
            marker = f'id="rd-matches-{row}-{col}-{layer}"'
            if n_segs == 0:
                continue
            i = html.find(marker)
            assert i != -1
            group = html[i : html.find('</g>', i)]
            assert group.count('<path') + group.count('<use') == n_segs


def test_reverse_strand_segments_serialized(tmp_path):
    """RC-only homology lands in the 'rev' layer with fwd empty."""
    idx = SequenceIndex(k=6)
    idx.add_sequence('fwd_seq', 'ATGCTAGCTAGGATCGGATC')
    idx.add_sequence('rc_seq', 'GATCCGATCCTAGCTAGCAT')
    out = tmp_path / 'rc.html'
    fig = DotPlotter(idx).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    panel = payload['panels']['rd-panel-0-1']
    assert panel['segments']['rev']
    # Segment coordinate rows are [qs, qe, ts, te] ints within bounds.
    for qs, qe, ts, te in panel['segments']['rev']:
        assert all(isinstance(v, int) for v in (qs, qe, ts, te))
        assert 0 <= qs < qe <= panel['qlen']
        assert 0 <= ts < te <= panel['tlen']


def test_empty_match_panel_serialized(tmp_path):
    """A pair with no shared k-mers yields a panel with empty segment lists."""
    idx = SequenceIndex(k=10)
    idx.add_sequence('a', 'ACGTACGTACGTACGTACGTACGT')
    idx.add_sequence('b', 'GGGGGGGGGGGGGGGGGGGGGGGG')
    out = tmp_path / 'empty.html'
    fig = DotPlotter(idx).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    panel = payload['panels']['rd-panel-0-1']
    assert panel['segments'] == {'fwd': [], 'rev': [], 'identity': []}


# ---------------------------------------------------------------------------
# Sequence embedding
# ---------------------------------------------------------------------------


def test_sequences_embedded_for_sequence_index(html_index, tmp_path):
    out = tmp_path / 'report.html'
    fig = DotPlotter(html_index).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    assert payload['has_sequences'] is True
    for panel in payload['panels'].values():
        for layer in ('fwd', 'rev'):
            seqs = panel['seqs'][layer]
            assert len(seqs) == len(panel['segments'][layer])
            for seq, seg in zip(seqs, panel['segments'][layer]):
                assert len(seq) == seg[1] - seg[0]
                # Query-side slice of the stored sequence.
                full = html_index.get_sequence(panel['query'])
                assert seq == full[seg[0] : seg[1]]


def test_paf_only_report_has_no_sequences(tmp_path):
    """PafAlignment lacks get_sequence, so has_sequences is False."""
    aln = PafAlignment.from_records(_make_paf_records())
    out = tmp_path / 'paf.html'
    fig = DotPlotter(aln).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    assert payload['has_sequences'] is False
    assert all('seqs' not in p for p in payload['panels'].values())
    # But coordinates are still present for the q1 vs t1 panel.
    assert _segment_count(payload) > 0


def test_sequence_embedding_cap(html_index, tmp_path):
    """Above the residue cap, sequences are dropped from the payload."""
    out = tmp_path / 'report.html'
    plotter = DotPlotter(html_index)

    # Rebuild the payload from a real capture with a tiny cap.
    fig = plotter.to_html(out)
    plt.close(fig)
    assert _read_payload(out)['has_sequences'] is True

    # Direct serializer call with max_residues=0 must omit sequences.
    capture = {
        'panels': {
            'rd-panel-0-0': {
                'query': 'seq1',
                'target': 'seq1',
                'query_id': 'seq1',
                'qlen': 20,
                'tlen': 20,
                'segments': {'fwd': [[0, 20, 0, 20]], 'rev': [], 'identity': []},
            }
        }
    }
    payload = build_panel_payload(
        capture, get_sequence=html_index.get_sequence, max_residues=0
    )
    assert payload['has_sequences'] is False
    assert 'seqs' not in payload['panels']['rd-panel-0-0']


# ---------------------------------------------------------------------------
# Identity mode
# ---------------------------------------------------------------------------


def test_identity_mode_serializes_identity_layer(tmp_path):
    aln = PafAlignment.from_records(_make_paf_records())
    out = tmp_path / 'identity.html'
    fig = DotPlotter(aln).to_html(out, color_by_identity=True)
    plt.close(fig)

    payload = _read_payload(out)
    q_row = [
        p
        for p in payload['panels'].values()
        if p['query'] == 'q1' and p['target'] == 't1'
    ]
    assert len(q_row) == 1
    identity = q_row[0]['segments']['identity']
    assert len(identity) == 2
    for qs, qe, ts, te, ident, strand in identity:
        assert 0 <= qs < qe
        assert 0 <= ts < te
        assert 0.0 <= ident <= 1.0
        assert strand in ('+', '-')
    # fwd/rev layers stay empty in identity mode.
    assert q_row[0]['segments']['fwd'] == []
    assert q_row[0]['segments']['rev'] == []
    assert pytest.approx(identity[0][4], abs=1e-3) == 0.9


# ---------------------------------------------------------------------------
# Edge-case sequences
# ---------------------------------------------------------------------------


def test_n_run_sequences(tmp_path):
    """Sequences with N runs render and serialize without error."""
    idx = SequenceIndex(k=5)
    idx.add_sequence('n_seq', 'ACGTANNNNNACGTACGTANNNNNACGTA')
    idx.add_sequence('other', 'ACGTAACGTACGTAACGTA')
    out = tmp_path / 'nrun.html'
    fig = DotPlotter(idx).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    assert len(payload['panels']) == 4
    # Embedded sequences reflect the raw stored bases (may include N).
    if payload['has_sequences']:
        for panel in payload['panels'].values():
            for layer in ('fwd', 'rev'):
                for seq in panel['seqs'][layer]:
                    assert set(seq) <= set('ACGTN')


def test_palindromic_sequence(tmp_path):
    """A reverse-complement palindrome self-plot has both strand layers."""
    palindrome = 'ACGTACGTACGT' + 'ACGTACGTACGT'[::-1].translate(
        str.maketrans('ACGT', 'TGCA')
    )
    idx = SequenceIndex(k=6)
    idx.add_sequence('pal', palindrome)
    out = tmp_path / 'pal.html'
    fig = DotPlotter(idx).to_html(out)
    plt.close(fig)

    payload = _read_payload(out)
    panel = payload['panels']['rd-panel-0-0']
    assert panel['segments']['fwd']  # the self diagonal
    assert panel['segments']['rev']  # the palindromic anti-diagonal


# ---------------------------------------------------------------------------
# Unsupported paths + state hygiene
# ---------------------------------------------------------------------------


def test_plot_single_html_raises(html_index, tmp_path):
    plotter = DotPlotter(html_index)
    with pytest.raises(ValueError, match='HTML output'):
        plotter.plot_single('seq1', 'seq2', output_path=tmp_path / 'one.html')
    plt.close('all')


def test_plot_identity_colorbar_html_raises(html_index, tmp_path):
    plotter = DotPlotter(html_index)
    with pytest.raises(ValueError, match='HTML output'):
        plotter.plot_identity_colorbar(output_path=tmp_path / 'cbar.html')
    plt.close('all')


def test_capture_reset_between_plots(html_index, tmp_path):
    """HTML capture does not leak into subsequent non-HTML plots."""
    plotter = DotPlotter(html_index)
    fig1 = plotter.to_html(tmp_path / 'a.html')
    plt.close(fig1)
    assert plotter._html_capture is None

    fig2 = plotter.plot(output_path=tmp_path / 'b.png')
    plt.close(fig2)
    assert plotter._html_capture is None
    # Non-HTML figures carry no gid tags.
    assert all(ax.get_gid() is None for ax in fig2.axes)


# ---------------------------------------------------------------------------
# Review regression tests
# ---------------------------------------------------------------------------


def _revcomp(seq):
    """Reverse complement of an uppercase ACGTN string."""
    return seq.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def test_reverse_contigs_sequences_match_display_orientation(html_index, tmp_path):
    """Mirrored panels embed the revcomp of the original mirrored region."""
    out = tmp_path / 'rev.html'
    fig = DotPlotter(html_index).to_html(out, reverse_contigs={'seq1'})
    plt.close(fig)

    payload = _read_payload(out)
    assert payload['has_sequences'] is True
    for panel in payload['panels'].values():
        full = html_index.get_sequence(panel['query'])
        qlen = panel['qlen']
        mirrored = panel['query'] == 'seq1'
        for layer in ('fwd', 'rev'):
            for seq, seg in zip(panel['seqs'][layer], panel['segments'][layer]):
                qs, qe = seg[0], seg[1]
                if mirrored:
                    assert seq == _revcomp(full[qlen - qe : qlen - qs])
                else:
                    assert seq == full[qs:qe]


def test_rasterized_true_still_produces_clickable_paths(html_index, tmp_path):
    """HTML output forces vector match layers even with rasterized=True."""
    out = tmp_path / 'raster.html'
    fig = DotPlotter(html_index).to_html(out, rasterized=True)
    plt.close(fig)

    payload = _read_payload(out)
    html = out.read_text()
    for gid, panel in payload['panels'].items():
        row, col = gid.rsplit('-', 2)[1:]
        for layer in ('fwd', 'rev'):
            n_segs = len(panel['segments'][layer])
            if n_segs == 0:
                continue
            i = html.find(f'id="rd-matches-{row}-{col}-{layer}"')
            assert i != -1
            group = html[i : html.find('</g>', i)]
            assert '<image' not in group
            assert group.count('<path') + group.count('<use') == n_segs


def test_title_is_html_escaped(html_index, tmp_path):
    """Markup in the report title is escaped, not injected."""
    out = tmp_path / 'esc.html'
    fig = DotPlotter(html_index).to_html(out, title='<script>alert(1)</script>')
    plt.close(fig)

    html = out.read_text()
    head = html.split('</head>')[0]
    assert '<script>alert(1)</script>' not in head
    assert '&lt;script&gt;alert(1)&lt;/script&gt;' in html


def test_stale_capture_does_not_leak_into_colorbar(html_index, tmp_path, monkeypatch):
    """An aborted HTML plot() must not let later figures reuse its capture."""
    plotter = DotPlotter(html_index)

    # Abort plot() mid-grid, after capture is active and partially filled.
    original = DotPlotter._plot_panel
    calls = {'n': 0}

    def exploding_panel(self, *args, **kwargs):
        calls['n'] += 1
        if calls['n'] >= 2:
            raise RuntimeError('boom')
        return original(self, *args, **kwargs)

    monkeypatch.setattr(DotPlotter, '_plot_panel', exploding_panel)
    with pytest.raises(RuntimeError, match='boom'):
        plotter.plot(output_path=tmp_path / 'aborted.html')
    monkeypatch.setattr(DotPlotter, '_plot_panel', original)
    plt.close('all')
    assert plotter._html_capture is not None  # leaked by the abort...

    # ...but the colorbar clears it and refuses HTML output.
    with pytest.raises(ValueError, match='HTML output'):
        plotter.plot_identity_colorbar(output_path=tmp_path / 'cbar.html')
    plt.close('all')
    assert plotter._html_capture is None


def test_track_entry_carries_stamped_uid():
    from rusty_dot._html.serialize import _track_entry
    from rusty_dot.annotation import GffAnnotation

    ann = GffAnnotation.from_text('c1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n')
    feat = ann.records[0]
    # Standalone renders never stamp one; the payload still has the key.
    assert _track_entry('x', 0, 0, feat)['uid'] == ''
    # The app's override pass stamps the table's uid onto the record.
    feat.uid = 'query:0'
    assert _track_entry('x', 0, 0, feat)['uid'] == 'query:0'
