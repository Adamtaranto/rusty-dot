"""Tests for the dotplot visualization module."""

import os

import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter


def _panel_segments(ax):
    """Return ``[(x0, y0, x1, y1), ...]`` for every match segment on ``ax``.

    ``DotPlotter._plot_panel`` renders matches with one or more
    :class:`matplotlib.collections.LineCollection` (a single artist holding all
    segments) rather than one :class:`~matplotlib.lines.Line2D` per match, so
    per-segment inspection goes through ``ax.collections`` instead of
    ``ax.lines``.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The panel axes to inspect.

    Returns
    -------
    list of tuple of float
        One ``(x0, y0, x1, y1)`` tuple per drawn line segment.
    """
    from matplotlib.collections import LineCollection

    segments = []
    for coll in ax.collections:
        if isinstance(coll, LineCollection):
            for seg in coll.get_segments():
                (x0, y0), (x1, y1) = seg[0], seg[-1]
                segments.append((x0, y0, x1, y1))
    return segments


@pytest.fixture
def dotplot_index():
    """Create an index for dotplot testing."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGTACGTACGTACGTACGT')
    idx.add_sequence('seq2', 'TACGTACGTACGTACGTACG')
    idx.add_sequence('seq3', 'GCGCGCGCGCGCGCGCGCGC')
    return idx


def test_dotplotter_creation(dotplot_index):
    """Test DotPlotter can be created."""
    plotter = DotPlotter(dotplot_index)
    assert plotter.index is dotplot_index


def test_plot_all_vs_all(dotplot_index, tmp_path):
    """Test all-vs-all dotplot generation."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'dotplot.png')
    plotter.plot(output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single(dotplot_index, tmp_path):
    """Test single pair dotplot generation."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'single.png')
    plotter.plot_single('seq1', 'seq2', output_path=output)
    assert os.path.exists(output)


def test_plot_subset(dotplot_index, tmp_path):
    """Test plotting a subset of sequences."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'subset.png')
    plotter.plot(
        query_names=['seq1', 'seq2'],
        target_names=['seq2', 'seq3'],
        output_path=output,
    )
    assert os.path.exists(output)


def test_plot_empty_index_raises():
    """Test that plotting with empty index raises ValueError."""
    idx = SequenceIndex(k=4)
    plotter = DotPlotter(idx)
    with pytest.raises(ValueError):
        plotter.plot()


def test_plot_scale_sequences(dotplot_index, tmp_path):
    """Test that scale_sequences=True produces a valid plot file."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'scaled.png')
    plotter.plot(output_path=output, scale_sequences=True)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_rc_only_sequences(tmp_path):
    """Sequences with only reverse-complement matches produce a non-empty dotplot."""
    # 'q' k-mer AAAC has revcomp GTTT which is in 't'; no forward shared k-mers
    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'AAACAAACAAAC')
    idx.add_sequence('t', 'GTTTGTTTGTTT')
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'rc_only.png')
    plotter.plot_single('q', 't', output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_rc_lines_are_anti_diagonal(tmp_path):
    """RC matches are drawn as anti-diagonal lines in the plot axes."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # 'q' = AAAC repeated; 't' = GTTT repeated (only RC matches, no forward)
    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'AAACAAACAAAC')
    idx.add_sequence('t', 'GTTTGTTTGTTT')
    plotter = DotPlotter(idx)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'q', 't')
    plt.close(fig)

    # At least one line must have been drawn
    segments = _panel_segments(ax)
    assert len(segments) > 0

    # For every RC match line, the x-data must be decreasing (anti-diagonal)
    for x0, _y0, x1, _y1 in segments:
        # anti-diagonal: x goes from t_end to t_start (decreasing)
        assert x0 >= x1, (
            f'Expected anti-diagonal line (x decreasing) for RC match, '
            f'got x=({x0}, {x1})'
        )


def test_plot_fwd_lines_are_diagonal(tmp_path):
    """Forward matches are drawn as diagonal lines (x increasing) in the axes."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'ACGTACGTACGT')
    idx.add_sequence('t', 'ACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'q', 't')
    plt.close(fig)

    fwd_segments = [
        (x0, y0, x1, y1) for x0, y0, x1, y1 in _panel_segments(ax) if x0 <= x1
    ]
    assert len(fwd_segments) > 0, 'Expected at least one forward (diagonal) match line'


def test_plot_rc_color_parameter(tmp_path):
    """The rc_color parameter is accepted and produces a valid plot."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'AAACAAACAAAC')
    idx.add_sequence('t', 'GTTTGTTTGTTT')
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'rc_color.png')
    plotter.plot_single('q', 't', output_path=output, rc_color='green')
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_scale_sequences_false_matches_default(dotplot_index, tmp_path):
    """Test that scale_sequences=True (default) and not passing it produce same-size files."""
    plotter = DotPlotter(dotplot_index)
    out1 = str(tmp_path / 'scaled1.png')
    out2 = str(tmp_path / 'scaled2.png')
    plotter.plot(output_path=out1, scale_sequences=True)
    plotter.plot(output_path=out2)  # default (True)
    # Both should exist and be non-empty
    assert os.path.getsize(out1) > 0
    assert os.path.getsize(out2) > 0


def test_plot_svg_format_extension(dotplot_index, tmp_path):
    """Test that SVG output is produced when the output path has a .svg extension."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'dotplot.svg')
    plotter.plot(output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0
    # SVG files start with XML/SVG markup
    with open(output, 'r') as f:
        content = f.read(100)
    assert '<svg' in content or '<?xml' in content


def test_plot_single_svg_format_extension(dotplot_index, tmp_path):
    """Test that plot_single produces SVG output via file extension."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'single.svg')
    plotter.plot_single('seq1', 'seq2', output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0
    with open(output, 'r') as f:
        content = f.read(100)
    assert '<svg' in content or '<?xml' in content


def test_plot_svg_explicit_format_parameter(dotplot_index, tmp_path):
    """Test that format='svg' forces SVG output even with a .png extension."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'dotplot.png')
    plotter.plot(output_path=output, format='svg')
    assert os.path.exists(output)
    with open(output, 'r') as f:
        content = f.read(100)
    assert '<svg' in content or '<?xml' in content


def test_plot_single_svg_explicit_format_parameter(dotplot_index, tmp_path):
    """Test that format='svg' in plot_single forces SVG output."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'single.png')
    plotter.plot_single('seq1', 'seq2', output_path=output, format='svg')
    assert os.path.exists(output)
    with open(output, 'r') as f:
        content = f.read(100)
    assert '<svg' in content or '<?xml' in content


def test_plot_min_length_filters_short_matches(tmp_path):
    """Setting min_length filters out short matches from the dotplot."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # Sequences with only short (k=4) matches when unmerged
    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'ACGTACGTACGT')
    idx.add_sequence('t', 'ACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig, ax_all = plt.subplots()
    plotter._plot_panel(ax_all, 'q', 't', min_length=0)
    n_lines_all = len(_panel_segments(ax_all))
    plt.close(fig)

    fig, ax_filtered = plt.subplots()
    # A very large min_length should filter everything out
    plotter._plot_panel(ax_filtered, 'q', 't', min_length=10000)
    n_lines_filtered = len(_panel_segments(ax_filtered))
    plt.close(fig)

    assert n_lines_all > 0
    assert n_lines_filtered == 0


def test_plot_min_length_parameter_accepted(dotplot_index, tmp_path):
    """Test that min_length parameter is accepted in plot() without error."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'minlen.png')
    plotter.plot(output_path=output, min_length=5)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single_min_length_parameter_accepted(dotplot_index, tmp_path):
    """Test that min_length parameter is accepted in plot_single() without error."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'minlen_single.png')
    plotter.plot_single('seq1', 'seq2', output_path=output, min_length=5)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0

    """In a multi-panel grid plot, x-labels appear only on the bottom row
    and y-labels appear only on the leftmost column."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    idx = SequenceIndex(k=4)
    idx.add_sequence('q1', 'ACGTACGTACGT')
    idx.add_sequence('q2', 'TACGTACGTACG')
    idx.add_sequence('t1', 'ACGTACGTACGT')
    idx.add_sequence('t2', 'GCGCGCGCGCGC')
    plotter = DotPlotter(idx)

    query_names = ['q1', 'q2']
    target_names = ['t1', 't2']
    nrows, ncols = len(query_names), len(target_names)

    fig, axes = plt.subplots(nrows, ncols, squeeze=False)
    for row_idx, q_name in enumerate(query_names):
        for col_idx, t_name in enumerate(target_names):
            plotter._plot_panel(
                axes[row_idx][col_idx],
                q_name,
                t_name,
                show_xlabel=(row_idx == nrows - 1),
                show_ylabel=(col_idx == 0),
            )
    plt.close(fig)

    # Bottom row: xlabel present
    assert axes[nrows - 1][0].get_xlabel() != '', 'bottom-left should have x-label'
    assert axes[nrows - 1][ncols - 1].get_xlabel() != '', (
        'bottom-right should have x-label'
    )
    # Non-bottom rows: xlabel absent
    assert axes[0][0].get_xlabel() == '', 'top-left should not have x-label'
    assert axes[0][ncols - 1].get_xlabel() == '', 'top-right should not have x-label'
    # Leftmost column: ylabel present
    assert axes[0][0].get_ylabel() != '', 'top-left should have y-label'
    assert axes[nrows - 1][0].get_ylabel() != '', 'bottom-left should have y-label'
    # Non-leftmost columns: ylabel absent
    assert axes[0][ncols - 1].get_ylabel() == '', 'top-right should not have y-label'
    assert axes[nrows - 1][ncols - 1].get_ylabel() == '', (
        'bottom-right should not have y-label'
    )


# ---------------------------------------------------------------------------
# Inline / return-value tests (Jupyter notebook support)
# ---------------------------------------------------------------------------


def test_plot_returns_figure(dotplot_index, tmp_path):
    """plot() returns a matplotlib Figure instance."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot(output_path=str(tmp_path / 'out.png'))
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_single_returns_figure(dotplot_index, tmp_path):
    """plot_single() returns a matplotlib Figure instance."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_single('seq1', 'seq2', output_path=str(tmp_path / 'out.png'))
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_no_output_path_returns_figure(dotplot_index):
    """plot() with output_path=None returns a Figure and writes no file."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_single_no_output_path_returns_figure(dotplot_index):
    """plot_single() with output_path=None returns a Figure and writes no file."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_single('seq1', 'seq2')
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_no_output_path_creates_no_file(dotplot_index, tmp_path, monkeypatch):
    """plot() with output_path=None does not create any file."""
    monkeypatch.chdir(tmp_path)
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot()
    plt.close(fig)
    # No files should have been created in the working directory
    assert list(tmp_path.iterdir()) == []


def test_plot_single_no_output_path_creates_no_file(
    dotplot_index, tmp_path, monkeypatch
):
    """plot_single() with output_path=None does not create any file."""
    monkeypatch.chdir(tmp_path)
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_single('seq1', 'seq2')
    plt.close(fig)
    assert list(tmp_path.iterdir()) == []


# ---------------------------------------------------------------------------
# Identity colouring tests
# ---------------------------------------------------------------------------


def _make_paf_alignment(query_name='seq1', target_name='seq2'):
    """Build a small PafAlignment with varying identity records."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    records = [
        PafRecord(
            query_name=query_name,
            query_len=20,
            query_start=0,
            query_end=10,
            strand='+',
            target_name=target_name,
            target_len=20,
            target_start=0,
            target_end=10,
            residue_matches=9,  # 90% identity
            alignment_block_len=10,
            mapping_quality=255,
        ),
        PafRecord(
            query_name=query_name,
            query_len=20,
            query_start=12,
            query_end=20,
            strand='-',
            target_name=target_name,
            target_len=20,
            target_start=12,
            target_end=20,
            residue_matches=8,  # identity = 8/8 = 100%
            alignment_block_len=8,
            mapping_quality=255,
        ),
    ]
    return PafAlignment(records)


def test_dotplotter_accepts_paf_alignment(dotplot_index):
    """DotPlotter can be constructed with a paf_alignment argument."""
    paf = _make_paf_alignment()
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)
    assert plotter.paf_alignment is paf


def test_dotplotter_paf_alignment_defaults_none(dotplot_index):
    """paf_alignment defaults to None when not supplied."""
    plotter = DotPlotter(dotplot_index)
    assert plotter.paf_alignment is None


def test_plot_color_by_identity_with_paf(dotplot_index, tmp_path):
    """color_by_identity=True with a PafAlignment produces a valid plot file."""
    import matplotlib

    matplotlib.use('Agg')

    paf = _make_paf_alignment(query_name='seq1', target_name='seq2')
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)
    output = str(tmp_path / 'identity.png')
    fig = plotter.plot_single(
        'seq1', 'seq2', output_path=output, color_by_identity=True
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_color_by_identity_lines_drawn(dotplot_index):
    """color_by_identity=True draws lines from PAF records, not k-mer matches."""
    import matplotlib

    matplotlib.use('Agg')

    paf = _make_paf_alignment(query_name='seq1', target_name='seq2')
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'seq1', 'seq2', color_by_identity=True)
    plt.close(fig)

    # Two PAF records → two line segments drawn
    assert len(_panel_segments(ax)) == 2


def test_plot_color_by_identity_uses_de_tag(dotplot_index):
    """Records with equal cols-10/11 identity but different de tags differ."""
    import dataclasses

    import matplotlib

    from rusty_dot.paf_io import PafAlignment

    matplotlib.use('Agg')

    base = _make_paf_alignment(query_name='seq1', target_name='seq2').records[0]
    identical = dataclasses.replace(base)
    diverged = dataclasses.replace(
        base,
        query_start=12,
        query_end=20,
        target_start=12,
        target_end=20,
        tags={'de': 0.4},
    )
    paf = PafAlignment([identical, diverged])
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'seq1', 'seq2', color_by_identity=True)
    colors = [
        coll.get_colors() for coll in ax.collections if len(coll.get_segments()) == 2
    ][0]
    plt.close(fig)

    # Same 9/10 blast identity, but the de tag overrides for the second
    # record (identity 0.6) so the two segments get different colours.
    assert not (colors[0] == colors[1]).all()


def test_plot_color_by_identity_custom_palette(dotplot_index, tmp_path):
    """identity_palette parameter is accepted without error."""
    import matplotlib

    matplotlib.use('Agg')

    paf = _make_paf_alignment(query_name='seq1', target_name='seq2')
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)
    output = str(tmp_path / 'plasma.png')
    fig = plotter.plot_single(
        'seq1',
        'seq2',
        output_path=output,
        color_by_identity=True,
        identity_palette='plasma',
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_color_by_identity_warns_without_paf(dotplot_index, caplog):
    """color_by_identity=True without a PafAlignment logs a warning."""
    import logging

    import matplotlib

    matplotlib.use('Agg')

    plotter = DotPlotter(dotplot_index)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        fig, ax = plt.subplots()
        plotter._plot_panel(ax, 'seq1', 'seq2', color_by_identity=True)
        plt.close(fig)

    assert any('color_by_identity' in msg for msg in caplog.messages)


def test_plot_color_by_identity_fallback_uses_kmer_matches(dotplot_index):
    """After warning, the panel falls back to k-mer match lines (not zero lines)."""
    import logging

    import matplotlib

    matplotlib.use('Agg')

    plotter = DotPlotter(dotplot_index)

    fig_kmer, ax_kmer = plt.subplots()
    plotter._plot_panel(ax_kmer, 'seq1', 'seq2', color_by_identity=False)
    n_kmer = len(_panel_segments(ax_kmer))
    plt.close(fig_kmer)

    # Suppress the expected warning to keep test output clean.
    fig_warn, ax_warn = plt.subplots()
    logging.disable(logging.WARNING)
    try:
        plotter._plot_panel(ax_warn, 'seq1', 'seq2', color_by_identity=True)
    finally:
        logging.disable(logging.NOTSET)
    n_fallback = len(_panel_segments(ax_warn))
    plt.close(fig_warn)

    # fallback should produce the same number of segments as the normal k-mer plot
    assert n_fallback == n_kmer
    assert n_kmer > 0


def test_plot_identity_colorbar_returns_figure(dotplot_index):
    """plot_identity_colorbar() returns a matplotlib Figure."""
    import matplotlib

    matplotlib.use('Agg')

    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_identity_colorbar()
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_identity_colorbar_saves_file(dotplot_index, tmp_path):
    """plot_identity_colorbar() saves a file when output_path is given."""
    import matplotlib

    matplotlib.use('Agg')

    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'colorbar.png')
    fig = plotter.plot_identity_colorbar(output_path=output)
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_identity_colorbar_custom_palette(dotplot_index, tmp_path):
    """plot_identity_colorbar() accepts a custom palette name."""
    import matplotlib

    matplotlib.use('Agg')

    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'colorbar_plasma.png')
    fig = plotter.plot_identity_colorbar(palette='plasma', output_path=output)
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_identity_colorbar_no_file_without_output(
    dotplot_index, tmp_path, monkeypatch
):
    """plot_identity_colorbar() with output_path=None does not create any file."""
    import matplotlib

    matplotlib.use('Agg')

    monkeypatch.chdir(tmp_path)
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_identity_colorbar()
    plt.close(fig)
    assert list(tmp_path.iterdir()) == []


def test_plot_color_by_identity_grid(dotplot_index, tmp_path):
    """color_by_identity=True works in multi-panel grid plot()."""
    import matplotlib

    matplotlib.use('Agg')

    paf_records = []
    from rusty_dot.paf_io import PafAlignment, PafRecord

    for q in ['seq1', 'seq2']:
        for t in ['seq1', 'seq2']:
            paf_records.append(
                PafRecord(
                    query_name=q,
                    query_len=20,
                    query_start=0,
                    query_end=10,
                    strand='+',
                    target_name=t,
                    target_len=20,
                    target_start=0,
                    target_end=10,
                    residue_matches=10,
                    alignment_block_len=10,
                    mapping_quality=255,
                )
            )
    paf = PafAlignment(paf_records)
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)
    output = str(tmp_path / 'grid_identity.png')
    fig = plotter.plot(
        query_names=['seq1', 'seq2'],
        target_names=['seq1', 'seq2'],
        output_path=output,
        color_by_identity=True,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_color_by_identity_min_length_filters(dotplot_index):
    """min_length filtering works correctly with identity-coloured PAF records."""
    import matplotlib

    matplotlib.use('Agg')

    from rusty_dot.paf_io import PafAlignment, PafRecord

    records = [
        PafRecord(
            query_name='seq1',
            query_len=20,
            query_start=0,
            query_end=5,
            strand='+',
            target_name='seq2',
            target_len=20,
            target_start=0,
            target_end=5,
            residue_matches=5,
            alignment_block_len=5,
            mapping_quality=255,
        ),
        PafRecord(
            query_name='seq1',
            query_len=20,
            query_start=10,
            query_end=20,
            strand='+',
            target_name='seq2',
            target_len=20,
            target_start=10,
            target_end=20,
            residue_matches=10,
            alignment_block_len=10,
            mapping_quality=255,
        ),
    ]
    paf = PafAlignment(records)
    plotter = DotPlotter(dotplot_index, paf_alignment=paf)

    fig, ax_all = plt.subplots()
    plotter._plot_panel(ax_all, 'seq1', 'seq2', color_by_identity=True, min_length=0)
    n_all = len(_panel_segments(ax_all))
    plt.close(fig)

    fig, ax_filtered = plt.subplots()
    plotter._plot_panel(
        ax_filtered, 'seq1', 'seq2', color_by_identity=True, min_length=8
    )
    n_filtered = len(_panel_segments(ax_filtered))
    plt.close(fig)

    assert n_all == 2
    assert n_filtered == 1


# ---------------------------------------------------------------------------
# PafAlignment as DotPlotter index (no SequenceIndex required)
# ---------------------------------------------------------------------------


def _make_paf_alignment_index():
    """Build a PafAlignment suitable for use as a DotPlotter index."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    return PafAlignment(
        [
            PafRecord(
                query_name='contigA',
                query_len=100,
                query_start=0,
                query_end=80,
                strand='+',
                target_name='contigB',
                target_len=90,
                target_start=0,
                target_end=80,
                residue_matches=78,
                alignment_block_len=80,
                mapping_quality=60,
            ),
            PafRecord(
                query_name='contigA',
                query_len=100,
                query_start=5,
                query_end=40,
                strand='-',
                target_name='contigC',
                target_len=50,
                target_start=10,
                target_end=45,
                residue_matches=33,
                alignment_block_len=35,
                mapping_quality=60,
            ),
        ]
    )


def test_paf_alignment_as_index_sequence_names():
    """PafAlignment.sequence_names() returns all unique query+target names."""

    aln = _make_paf_alignment_index()
    names = aln.sequence_names()
    assert set(names) == {'contigA', 'contigB', 'contigC'}


def test_paf_alignment_as_index_get_sequence_length():
    """PafAlignment.get_sequence_length() returns lengths from PAF records."""
    aln = _make_paf_alignment_index()
    assert aln.get_sequence_length('contigA') == 100
    assert aln.get_sequence_length('contigB') == 90
    assert aln.get_sequence_length('contigC') == 50


def test_paf_alignment_as_index_get_sequence_length_missing():
    """PafAlignment.get_sequence_length() raises KeyError for unknown names."""
    aln = _make_paf_alignment_index()
    with pytest.raises(KeyError):
        aln.get_sequence_length('nonexistent')


def test_dotplotter_accepts_paf_alignment_as_index(tmp_path):
    """DotPlotter can be constructed with a PafAlignment as the index."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)
    # paf_alignment is automatically set to the PafAlignment index
    assert plotter.paf_alignment is aln


def test_dotplotter_paf_as_index_plot(tmp_path):
    """DotPlotter(PafAlignment) produces a valid plot file without SequenceIndex."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)
    output = str(tmp_path / 'paf_index_plot.png')
    fig = plotter.plot(
        query_names=['contigA'],
        target_names=['contigB'],
        output_path=output,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_dotplotter_paf_as_index_strand_colors(tmp_path):
    """DotPlotter(PafAlignment) draws forward and RC alignments with correct colours."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'contigA', 'contigB', dot_color='blue', rc_color='red')
    plt.close(fig)

    # One '+' alignment should produce one line segment
    assert len(_panel_segments(ax)) == 1


def test_dotplotter_paf_as_index_rc_alignment():
    """DotPlotter(PafAlignment) renders RC alignment as anti-diagonal line."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'contigA', 'contigC')
    plt.close(fig)

    segments = _panel_segments(ax)
    assert len(segments) == 1
    x0, _y0, x1, _y1 = segments[0]
    # RC strand: x goes from target_end to target_start (decreasing)
    assert x0 > x1


def test_dotplotter_paf_as_index_plot_grid(tmp_path):
    """DotPlotter(PafAlignment).plot() works for multi-panel grids."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    q_order, t_order = aln.reorder_contigs()
    plotter = DotPlotter(aln)
    output = str(tmp_path / 'paf_grid.png')
    fig = plotter.plot(
        query_names=q_order,
        target_names=t_order,
        output_path=output,
        scale_sequences=True,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_dotplotter_paf_as_index_identity_coloring(tmp_path):
    """DotPlotter(PafAlignment).plot(color_by_identity=True) works without explicit paf_alignment."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)
    output = str(tmp_path / 'paf_identity.png')
    fig = plotter.plot(
        query_names=['contigA'],
        target_names=['contigB'],
        output_path=output,
        color_by_identity=True,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


# ---------------------------------------------------------------------------
# DotPlotter + CrossIndex integration (query_group / target_group)
# ---------------------------------------------------------------------------


def _make_cross_index():
    """Build a small CrossIndex with two groups and pre-computed matches."""
    import matplotlib

    matplotlib.use('Agg')

    from rusty_dot.paf_io import CrossIndex

    cross = CrossIndex(k=4)
    cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='group_a')
    cross.add_sequence('q2', 'TACGTACGTACGTACG', group='group_a')
    cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='group_b')
    return cross


def test_dotplotter_cross_index_query_group_param(tmp_path):
    """DotPlotter.plot() with query_group/target_group resolves names from CrossIndex."""
    import matplotlib

    matplotlib.use('Agg')

    cross = _make_cross_index()
    cross.compute_matches()
    plotter = DotPlotter(cross)
    output = str(tmp_path / 'cross_group.png')
    fig = plotter.plot(
        query_group='group_a',
        target_group='group_b',
        output_path=output,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_dotplotter_cross_index_group_not_cross_raises():
    """query_group/target_group raises ValueError when index is not CrossIndex."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('s1', 'ACGT' * 4)
    plotter = DotPlotter(idx)
    with pytest.raises(ValueError, match='CrossIndex'):
        plotter.plot(query_group='a', target_group='b')


def test_dotplotter_cross_index_uses_precomputed_records():
    """When compute_matches() was called, DotPlotter uses the cached records."""
    import matplotlib

    matplotlib.use('Agg')

    from rusty_dot.paf_io import CrossIndex, PafAlignment

    cross = CrossIndex(k=4)
    cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
    cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
    cross.compute_matches()
    plotter = DotPlotter(cross)

    fig, ax_pre = plt.subplots()
    plotter._plot_panel(
        ax_pre,
        'a:q1',
        'b:t1',
        paf_alignment_override=PafAlignment(cross.get_records_for_pair('a', 'b')),
    )
    plt.close(fig)
    # The panel rendered without error and axes have correct limits
    assert ax_pre.get_xlim()[1] > 0


def test_dotplotter_cross_index_label_strips_prefix():
    """Axis labels should not contain the group prefix 'group:' for CrossIndex."""
    import matplotlib

    matplotlib.use('Agg')

    from rusty_dot.paf_io import CrossIndex

    cross = CrossIndex(k=4)
    cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='grp_a')
    cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='grp_b')
    cross.compute_matches()
    plotter = DotPlotter(cross)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'grp_a:q1', 'grp_b:t1', show_xlabel=True, show_ylabel=True)
    plt.close(fig)

    # Labels should be plain names, not prefixed names
    assert ax.get_xlabel() == 't1', f'Expected "t1", got {ax.get_xlabel()!r}'
    assert ax.get_ylabel() == 'q1', f'Expected "q1", got {ax.get_ylabel()!r}'


def test_dotplotter_cross_index_strip_group_prefix_static():
    """_strip_group_prefix helper correctly strips 'group:' prefix."""
    plotter = DotPlotter(SequenceIndex(k=4))
    assert plotter._strip_group_prefix('group_a:seq1') == 'seq1'
    assert plotter._strip_group_prefix('seq1') == 'seq1'
    # Only the first ':' is treated as the group separator; subsequent colons
    # are part of the unprefixed name.
    assert plotter._strip_group_prefix('a:b:c') == 'b:c'
    # A name with a colon but no group prefix still returns everything after
    # the first colon — callers should only pass CrossIndex internal names or
    # plain names to this helper.
    assert plotter._strip_group_prefix('name:with:colons') == 'with:colons'


def test_dotplotter_cross_index_plot_single_with_groups(tmp_path):
    """plot_single() with query_group/target_group works for CrossIndex."""
    import matplotlib

    matplotlib.use('Agg')

    from rusty_dot.paf_io import CrossIndex

    cross = CrossIndex(k=4)
    cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='grp_a')
    cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='grp_b')
    cross.compute_matches()
    plotter = DotPlotter(cross)
    output = str(tmp_path / 'cross_single.png')
    fig = plotter.plot_single(
        'q1',
        't1',
        query_group='grp_a',
        target_group='grp_b',
        output_path=output,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_dotplotter_cross_index_plot_single_group_not_cross_raises():
    """plot_single() with groups raises ValueError when index is not CrossIndex."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('s1', 'ACGT' * 4)
    plotter = DotPlotter(idx)
    with pytest.raises(ValueError, match='CrossIndex'):
        plotter.plot_single('s1', 's1', query_group='a', target_group='b')


def test_dotplotter_resolve_group_names_no_groups():
    """_resolve_group_names returns inputs unchanged when no groups provided."""
    plotter = DotPlotter(SequenceIndex(k=4))
    q, t, paf = plotter._resolve_group_names(None, None, ['a'], ['b'])
    assert q == ['a']
    assert t == ['b']
    assert paf is None


def test_dotplotter_resolve_group_names_raises_for_non_cross():
    """_resolve_group_names raises ValueError when groups given but index is not CrossIndex."""
    plotter = DotPlotter(SequenceIndex(k=4))
    with pytest.raises(ValueError, match='CrossIndex'):
        plotter._resolve_group_names('x', 'y', None, None)


# --- _resolve_rasterized -----------------------------------------------------


def test_resolve_rasterized_explicit_bool():
    """Explicit True/False are returned verbatim regardless of count."""
    from rusty_dot.dotplot import _resolve_rasterized

    assert _resolve_rasterized(1_000_000, True, 50) is True
    assert _resolve_rasterized(1_000_000, False, 50) is False


def test_resolve_rasterized_auto_threshold():
    """'auto' rasterizes only above the threshold."""
    from rusty_dot.dotplot import _resolve_rasterized

    assert _resolve_rasterized(10, 'auto', 50) is False
    assert _resolve_rasterized(51, 'auto', 50) is True
    assert _resolve_rasterized(50, 'auto', 50) is False  # boundary: not above


def test_resolve_rasterized_invalid_string():
    """An unknown string raises ValueError."""
    from rusty_dot.dotplot import _resolve_rasterized

    with pytest.raises(ValueError, match='auto'):
        _resolve_rasterized(1, 'sometimes', 50)


# --- _chain_blocks -----------------------------------------------------------


def test_chain_blocks_gap_zero_is_noop():
    """chain_gap=0 returns the blocks unchanged."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (100, 110, 100, 110, '+')]
    assert _chain_blocks(blocks, 0) == blocks


def test_chain_blocks_merges_same_diagonal_within_gap():
    """Collinear forward blocks within the gap merge into one."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (12, 20, 12, 20, '+')]  # diagonal 0, gap 2
    chained = _chain_blocks(blocks, 5)
    assert chained == [(0, 20, 0, 20, '+')]


def test_chain_blocks_keeps_blocks_beyond_gap_separate():
    """Collinear blocks farther apart than the gap stay separate."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (12, 20, 12, 20, '+')]  # gap 2
    chained = _chain_blocks(blocks, 1)  # gap tolerance 1 < 2
    assert len(chained) == 2


def test_chain_blocks_different_diagonals_stay_separate():
    """Blocks on different diagonals are never chained."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (0, 10, 50, 60, '+')]  # diagonals 0 and 50
    chained = _chain_blocks(blocks, 1000)
    assert len(chained) == 2


def test_chain_blocks_merges_overlapping():
    """Overlapping collinear blocks merge into their union."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (5, 15, 5, 15, '+')]
    chained = _chain_blocks(blocks, 1)
    assert chained == [(0, 15, 0, 15, '+')]


def test_chain_blocks_reverse_strand_antidiagonal():
    """Reverse-complement blocks on the same anti-diagonal chain correctly."""
    from rusty_dot.dotplot import _chain_blocks

    # Anti-diagonal invariant q_start + t_end == 100 for both.
    blocks = [(0, 10, 90, 100, '-'), (12, 20, 78, 88, '-')]  # gap 2
    chained = _chain_blocks(blocks, 5)
    assert chained == [(0, 20, 78, 100, '-')]


def test_chain_blocks_does_not_merge_across_strands():
    """A '+' and a '-' block are never chained together."""
    from rusty_dot.dotplot import _chain_blocks

    blocks = [(0, 10, 0, 10, '+'), (12, 20, 12, 20, '-')]
    chained = _chain_blocks(blocks, 1000)
    assert len(chained) == 2


# --- auto-rasterization end to end -------------------------------------------


def test_plot_panel_auto_vector_below_threshold():
    """With rasterized='auto' and few segments, the match layer stays vector."""
    import matplotlib

    matplotlib.use('Agg')

    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'ACGTACGTACGT')
    idx.add_sequence('t', 'ACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'q', 't', rasterized='auto', rasterization_threshold=50_000)
    from matplotlib.collections import LineCollection

    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    assert lcs, 'expected at least one match LineCollection'
    assert all(lc.get_rasterized() is False for lc in lcs)
    plt.close(fig)


def test_plot_panel_auto_rasterizes_above_threshold():
    """With a tiny threshold, the match layer is rasterized."""
    import matplotlib

    matplotlib.use('Agg')

    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'ACGTACGTACGT')
    idx.add_sequence('t', 'ACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'q', 't', rasterized='auto', rasterization_threshold=0)
    from matplotlib.collections import LineCollection

    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    assert lcs
    assert all(lc.get_rasterized() is True for lc in lcs)
    plt.close(fig)


def test_plot_panel_explicit_rasterized_false():
    """rasterized=False forces vector even with a zero threshold."""
    import matplotlib

    matplotlib.use('Agg')

    idx = SequenceIndex(k=4)
    idx.add_sequence('q', 'ACGTACGTACGT')
    idx.add_sequence('t', 'ACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'q', 't', rasterized=False, rasterization_threshold=0)
    from matplotlib.collections import LineCollection

    lcs = [c for c in ax.collections if isinstance(c, LineCollection)]
    assert all(lc.get_rasterized() is False for lc in lcs)
    plt.close(fig)


def test_plot_svg_match_layer_is_vector(dotplot_index, tmp_path):
    """Default SVG output keeps the match layer as vector paths (no <image>)."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'vector.svg')
    fig = plotter.plot(output_path=output)
    plt.close(fig)
    with open(output) as fh:
        content = fh.read()
    # A rasterized match layer would embed a raster <image>; vector output uses
    # <path>/<use> elements only.
    assert '<image' not in content
    assert '<path' in content


def test_plot_panel_chain_gap_reduces_segments():
    """Chaining never increases, and can decrease, the number of drawn segments."""
    import matplotlib

    matplotlib.use('Agg')

    # A sequence vs itself: a long conserved diagonal of k-mer matches.
    seq = 'ACGT' * 50
    idx = SequenceIndex(k=6)
    idx.add_sequence('q', seq)
    idx.add_sequence('t', seq)
    plotter = DotPlotter(idx)

    fig, ax0 = plt.subplots()
    plotter._plot_panel(ax0, 'q', 't', merge=False, chain_gap=0)
    n_unchained = len(_panel_segments(ax0))
    plt.close(fig)

    fig, ax1 = plt.subplots()
    plotter._plot_panel(ax1, 'q', 't', merge=False, chain_gap=500)
    n_chained = len(_panel_segments(ax1))
    plt.close(fig)

    assert n_unchained > 0
    assert n_chained <= n_unchained


# ---------------------------------------------------------------------------
# Reverse-oriented contig rendering (GAP 3)
# ---------------------------------------------------------------------------


def _forward_paf_alignment():
    """A PafAlignment with a single forward match: query 0-50 → target 0-50."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    rec = PafRecord.from_line('q\t100\t0\t50\t+\tt\t100\t0\t50\t48\t50\t255')
    return PafAlignment.from_records([rec])


def test_reverse_contig_mirrors_query_axis():
    """reverse_contigs mirrors the query coordinates (q → q_len - q)."""
    aln = _forward_paf_alignment()
    plotter = DotPlotter(aln)

    fig_fwd = plotter.plot(query_names=['q'], target_names=['t'])
    fwd_segs = _panel_segments(fig_fwd.axes[0])
    plt.close(fig_fwd)

    fig_rev = plotter.plot(query_names=['q'], target_names=['t'], reverse_contigs={'q'})
    rev_segs = _panel_segments(fig_rev.axes[0])
    plt.close(fig_rev)

    assert fwd_segs and rev_segs
    # Forward match occupies query 0-50; the reversed rendering mirrors it to
    # query 50-100 (q_len = 100), so the maximum y-coordinate jumps to ~100.
    fwd_max_y = max(max(s[1], s[3]) for s in fwd_segs)
    rev_max_y = max(max(s[1], s[3]) for s in rev_segs)
    assert fwd_max_y == 50
    assert rev_max_y == 100


def test_reverse_contigs_empty_leaves_unchanged():
    """Passing an empty set leaves the panel identical to the default."""
    aln = _forward_paf_alignment()
    plotter = DotPlotter(aln)

    fig_a = plotter.plot(query_names=['q'], target_names=['t'])
    segs_a = _panel_segments(fig_a.axes[0])
    plt.close(fig_a)

    fig_b = plotter.plot(query_names=['q'], target_names=['t'], reverse_contigs=set())
    segs_b = _panel_segments(fig_b.axes[0])
    plt.close(fig_b)

    assert segs_a == segs_b


def test_reverse_contigs_auto_from_pafalignment():
    """When reverse_contigs is None the set is pulled from the PafAlignment."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    recs = [
        PafRecord.from_line('qr\t300\t0\t150\t-\tref\t1000\t850\t1000\t140\t150\t255'),
        PafRecord.from_line('qr\t300\t150\t300\t-\tref\t1000\t700\t850\t140\t150\t255'),
    ]
    aln = PafAlignment.from_records(recs)
    aln.reorder_contigs(['qr'], ['ref'])  # populates reversed_contigs
    assert aln.reversed_contigs == {'qr'}

    plotter = DotPlotter(aln)
    # No reverse_contigs argument → auto-pull flags 'qr' as reversed.
    fig = plotter.plot(query_names=['qr'], target_names=['ref'])
    segs = _panel_segments(fig.axes[0])
    plt.close(fig)
    assert segs  # rendered without error using the auto-detected reversed set


def test_hide_internal_axes_removes_internal_boundaries(dotplot_index):
    """hide_internal_axes=True hides internal ticks/spines, keeps outer frame."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot(hide_internal_axes=True)
    fig.canvas.draw()

    nrows = ncols = 3
    axes = [fig.axes[i * ncols : (i + 1) * ncols] for i in range(nrows)]
    for row_idx in range(nrows):
        for col_idx in range(ncols):
            ax = axes[row_idx][col_idx]
            pos = f'panel ({row_idx}, {col_idx})'

            # Spines: internal edges hidden, outer frame intact.
            assert ax.spines['top'].get_visible() == (row_idx == 0), pos
            assert ax.spines['bottom'].get_visible() == (row_idx == nrows - 1), pos
            assert ax.spines['left'].get_visible() == (col_idx == 0), pos
            assert ax.spines['right'].get_visible() == (col_idx == ncols - 1), pos

            # Tick marks: only on the bottom row / left column.
            x_ticks = [t.tick1line.get_visible() for t in ax.xaxis.get_major_ticks()]
            y_ticks = [t.tick1line.get_visible() for t in ax.yaxis.get_major_ticks()]
            assert all(v == (row_idx == nrows - 1) for v in x_ticks), pos
            assert all(v == (col_idx == 0) for v in y_ticks), pos

    plt.close(fig)


def test_hide_internal_axes_zero_panel_spacing(dotplot_index):
    """hide_internal_axes=True collapses inter-panel gaps to zero."""
    plotter = DotPlotter(dotplot_index)
    for scale_sequences in (True, False):
        fig = plotter.plot(hide_internal_axes=True, scale_sequences=scale_sequences)
        gs = fig.axes[0].get_subplotspec().get_gridspec()
        params = gs.get_subplot_params(fig)
        assert params.wspace == 0.0
        assert params.hspace == 0.0
        plt.close(fig)


def test_default_keeps_internal_axes(dotplot_index):
    """Regression: the default path keeps internal spines/ticks and only
    suppresses redundant internal tick labels."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot()
    fig.canvas.draw()

    nrows = ncols = 3
    axes = [fig.axes[i * ncols : (i + 1) * ncols] for i in range(nrows)]
    for row_idx in range(nrows):
        for col_idx in range(ncols):
            ax = axes[row_idx][col_idx]
            pos = f'panel ({row_idx}, {col_idx})'

            # All spines and tick marks stay visible on every panel.
            for side in ('top', 'bottom', 'left', 'right'):
                assert ax.spines[side].get_visible(), pos
            assert all(t.tick1line.get_visible() for t in ax.xaxis.get_major_ticks()), (
                pos
            )
            assert all(t.tick1line.get_visible() for t in ax.yaxis.get_major_ticks()), (
                pos
            )

            # Tick labels only on the bottom row (x) and left column (y).
            x_labels = [t.label1.get_visible() for t in ax.xaxis.get_major_ticks()]
            y_labels = [t.label1.get_visible() for t in ax.yaxis.get_major_ticks()]
            assert all(v == (row_idx == nrows - 1) for v in x_labels), pos
            assert all(v == (col_idx == 0) for v in y_labels), pos

    # Default spacing is not collapsed to zero.
    gs = fig.axes[0].get_subplotspec().get_gridspec()
    params = gs.get_subplot_params(fig)
    assert params.wspace > 0.0
    assert params.hspace > 0.0

    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot-time contig ordering (contig_order / auto_reverse)
# ---------------------------------------------------------------------------


def _random_seq(rng, n):
    """Return a deterministic pseudo-random DNA sequence of length *n*."""
    return ''.join(rng.choice('ACGT') for _ in range(n))


def _revcomp(seq):
    """Return the reverse complement of *seq*."""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp[b] for b in reversed(seq))


def _shuffled_cross_index():
    """Build a two-group CrossIndex with the query contigs deliberately shuffled.

    The target group ``'ref'`` holds one chromosome made of three distinct
    60 bp blocks; the query group ``'qry'`` holds one contig per block, added
    out of colinear order (``c2, c3, c1``).

    Returns
    -------
    CrossIndex
        The populated index (matches not yet computed).
    """
    import random

    from rusty_dot.paf_io import CrossIndex

    rng = random.Random(42)
    blocks = [_random_seq(rng, 60) for _ in range(3)]
    cross = CrossIndex(k=11)
    cross.add_sequence('chr1', ''.join(blocks), group='ref')
    cross.add_sequence('c2', blocks[1], group='qry')
    cross.add_sequence('c3', blocks[2], group='qry')
    cross.add_sequence('c1', blocks[0], group='qry')
    return cross


def _grid_query_ylabels(fig, ncols):
    """Return the y-axis labels of the first column, top-to-bottom."""
    return [ax.get_ylabel() for i, ax in enumerate(fig.axes) if i % ncols == 0]


def _grid_target_titles(fig, ncols):
    """Return the column titles of the top row, left-to-right."""
    return [fig.axes[i].get_title() for i in range(ncols)]


def test_contig_order_invalid_value_raises(dotplot_index):
    """An unknown contig_order value raises ValueError listing the options."""
    plotter = DotPlotter(dotplot_index)
    with pytest.raises(ValueError, match='contig_order'):
        plotter.plot(contig_order='bogus')


def test_contig_order_default_none_unchanged(dotplot_index):
    """contig_order=None keeps the default alphabetical ordering."""
    plotter = DotPlotter(dotplot_index)
    fig_default = plotter.plot()
    fig_none = plotter.plot(contig_order=None)
    try:
        assert _grid_target_titles(fig_none, 3) == _grid_target_titles(fig_default, 3)
        assert _grid_target_titles(fig_none, 3) == ['seq1', 'seq2', 'seq3']
    finally:
        plt.close(fig_default)
        plt.close(fig_none)


def test_contig_order_length_sorts_descending():
    """contig_order='length' orders panels by descending sequence length."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('short', 'ACGTACGTACGT')  # 12 bp
    idx.add_sequence('long', 'ACGTACGTACGTACGTACGTACGTACGT')  # 28 bp
    idx.add_sequence('mid', 'ACGTACGTACGTACGTACGT')  # 20 bp
    plotter = DotPlotter(idx)

    fig = plotter.plot(contig_order='length')
    try:
        expected = ['long', 'mid', 'short']
        assert _grid_target_titles(fig, 3) == expected
        assert _grid_query_ylabels(fig, 3) == expected
    finally:
        plt.close(fig)


def test_contig_order_length_explicit_names_take_precedence():
    """Explicit query_names keep their order; the unset axis is length-sorted."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('short', 'ACGTACGTACGT')
    idx.add_sequence('long', 'ACGTACGTACGTACGTACGTACGTACGT')
    idx.add_sequence('mid', 'ACGTACGTACGTACGTACGT')
    plotter = DotPlotter(idx)

    fig = plotter.plot(contig_order='length', query_names=['short', 'long'])
    try:
        # Explicit query order preserved (rows), targets length-sorted (cols).
        assert _grid_query_ylabels(fig, 3) == ['short', 'long']
        assert _grid_target_titles(fig, 3) == ['long', 'mid', 'short']
    finally:
        plt.close(fig)


def test_contig_order_length_crossindex_groups():
    """contig_order='length' on a CrossIndex reorders the plotted groups."""
    cross = _shuffled_cross_index()
    # Make the contigs distinguishable by length.
    cross.add_sequence('tiny', 'ACGTACGTACGT', group='qry')
    plotter = DotPlotter(cross)

    fig = plotter.plot(query_group='qry', target_group='ref', contig_order='length')
    try:
        labels = _grid_query_ylabels(fig, 1)
        assert labels[0] in {'c1', 'c2', 'c3'}  # all 60 bp, ties allowed
        assert labels[-1] == 'tiny'
    finally:
        plt.close(fig)
    # The group order itself was updated in place.
    assert cross.contig_order['qry'][-1] == 'tiny'


def test_contig_order_colinearity_crossindex_restores_order():
    """'colinearity' on a shuffled two-group CrossIndex restores colinear order."""
    cross = _shuffled_cross_index()
    assert cross.contig_order['qry'] == ['c2', 'c3', 'c1']  # shuffled on entry
    plotter = DotPlotter(cross)

    fig = plotter.plot(
        query_group='qry', target_group='ref', contig_order='colinearity'
    )
    try:
        assert _grid_query_ylabels(fig, 1) == ['c1', 'c2', 'c3']
    finally:
        plt.close(fig)
    assert cross.contig_order['qry'] == ['c1', 'c2', 'c3']


def test_contig_order_colinearity_sequenceindex_uses_optimal_order():
    """'colinearity' on a bare SequenceIndex uses optimal_contig_order."""
    import random

    rng = random.Random(7)
    blocks = [_random_seq(rng, 60) for _ in range(3)]
    idx = SequenceIndex(k=11)
    idx.add_sequence('chr1', ''.join(blocks))
    idx.add_sequence('c2', blocks[1])
    idx.add_sequence('c3', blocks[2])
    idx.add_sequence('c1', blocks[0])
    plotter = DotPlotter(idx)

    fig = plotter.plot(
        target_names=['chr1'],
        contig_order='colinearity',
    )
    try:
        labels = _grid_query_ylabels(fig, 1)
        # The three block contigs appear in colinear order along chr1.
        block_rows = [n for n in labels if n in {'c1', 'c2', 'c3'}]
        assert block_rows == ['c1', 'c2', 'c3']
    finally:
        plt.close(fig)


def test_contig_order_auto_reverse_mirrors_reversed_contig():
    """auto_reverse=True mirrors a reverse-oriented contig like reverse_contigs."""
    import random

    from rusty_dot.paf_io import CrossIndex

    rng = random.Random(99)
    chrom = _random_seq(rng, 120)
    cross = CrossIndex(k=11)
    # 'qry' added first so group auto-detection uses it as the query axis.
    cross.add_sequence('q1', _revcomp(chrom), group='qry')
    cross.add_sequence('chr1', chrom, group='ref')
    plotter = DotPlotter(cross)

    def _segs(**kwargs):
        # No explicit groups: 'colinearity' auto-detects the two groups.
        fig = plotter.plot(contig_order='colinearity', **kwargs)
        segs = _panel_segments(fig.axes[0])
        plt.close(fig)
        return segs

    auto_segs = _segs(auto_reverse=True)
    explicit_segs = _segs(reverse_contigs={'q1'})
    plain_segs = _segs(reverse_contigs=set())

    assert cross.reversed_contigs('qry') == {'q1'}
    assert auto_segs  # something was drawn
    # auto_reverse reproduces the explicit reverse_contigs rendering ...
    assert auto_segs == explicit_segs
    # ... and differs from the unmirrored rendering.
    assert auto_segs != plain_segs


def test_contig_order_auto_reverse_explicit_reverse_contigs_wins():
    """An explicit reverse_contigs argument overrides auto_reverse."""
    import random

    from rusty_dot.paf_io import CrossIndex

    rng = random.Random(99)
    chrom = _random_seq(rng, 120)
    cross = CrossIndex(k=11)
    cross.add_sequence('chr1', chrom, group='ref')
    cross.add_sequence('q1', _revcomp(chrom), group='qry')
    plotter = DotPlotter(cross)

    fig_override = plotter.plot(
        query_group='qry',
        target_group='ref',
        contig_order='colinearity',
        auto_reverse=True,
        reverse_contigs=set(),  # explicitly disable mirroring
    )
    segs_override = _panel_segments(fig_override.axes[0])
    plt.close(fig_override)

    fig_auto = plotter.plot(
        query_group='qry',
        target_group='ref',
        contig_order='colinearity',
        auto_reverse=True,
    )
    segs_auto = _panel_segments(fig_auto.axes[0])
    plt.close(fig_auto)

    assert segs_override != segs_auto


def test_plot_from_computed_cache_renders_reverse_segments():
    """Plotting via groups from the compute_matches cache shows '-' matches.

    Regression test: the forward-strand-only cache used to blank panels for
    reverse-oriented contigs when the pre-computed records were used for
    rendering.
    """
    import random

    from rusty_dot.paf_io import CrossIndex

    rng = random.Random(21)
    chrom = ''.join(rng.choice('ACGT') for _ in range(200))
    cross = CrossIndex(k=15)
    cross.add_sequence('chr1', chrom, group='ref')
    cross.add_sequence('q_rc', _revcomp(chrom), group='qry')
    cross.compute_matches('qry', 'ref')
    assert ('qry', 'ref') in cross.computed_group_pairs

    plotter = DotPlotter(cross)
    fig = plotter.plot(query_group='qry', target_group='ref')
    segs = _panel_segments(fig.axes[0])
    plt.close(fig)

    assert segs, 'expected segments rendered from the cached records'
    # Reverse matches are drawn anti-diagonal (x decreasing).
    assert any(x0 > x1 for x0, _y0, x1, _y1 in segs), (
        'expected at least one reverse (anti-diagonal) segment from the cache'
    )


def test_identity_colorbar_appends_key_axes():
    """identity_colorbar=True adds one colorbar axes; off by default."""
    import matplotlib.pyplot as plt

    from rusty_dot.paf_io import PafAlignment, PafRecord

    records = [
        PafRecord('q1', 4000, 0, 3000, '+', 't1', 5000, 0, 3000, 2400, 3000, 255),
        PafRecord('q1', 4000, 3000, 4000, '+', 't1', 5000, 3000, 4000, 990, 1000, 255),
    ]
    aln = PafAlignment(records)
    fig_plain = DotPlotter(aln).plot(color_by_identity=True)
    fig_keyed = DotPlotter(aln).plot(color_by_identity=True, identity_colorbar=True)
    assert len(fig_keyed.axes) == len(fig_plain.axes) + 1
    cbar_ax = fig_keyed.axes[-1]
    assert cbar_ax.get_ylabel() == 'Identity (%)'
    # identity_colorbar without color_by_identity is a no-op.
    fig_off = DotPlotter(aln).plot(identity_colorbar=True)
    assert len(fig_off.axes) == len(fig_plain.axes)
    plt.close(fig_plain)
    plt.close(fig_keyed)
    plt.close(fig_off)


# ------------------------------------------------- focused-view axis labels


@pytest.mark.parametrize(
    ('span', 'divisor', 'unit'),
    [
        (500, 1.0, 'bp'),
        (9_999, 1.0, 'bp'),
        (10_000, 1e3, 'Kbp'),
        (999_999, 1e3, 'Kbp'),
        (1_000_000, 1e6, 'Mbp'),
        (37_000_000, 1e6, 'Mbp'),
    ],
)
def test_bp_unit_thresholds(span, divisor, unit):
    from rusty_dot.dotplot import _bp_unit

    assert _bp_unit(span) == (divisor, unit)


def test_focused_pair_axis_labels_and_units(dotplot_index):
    """A 1x1 plot labels the axes with contig names + units, no rotated title."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot(query_names=['seq1'], target_names=['seq2'])
    ax = fig.axes[0]
    assert ax.get_xlabel() == 'seq2 (bp)'
    assert ax.get_ylabel() == 'seq1 (bp)'
    assert ax.get_title() == ''
    fig.canvas.draw()
    tick_texts = [t.get_text() for t in ax.get_xticklabels()]
    assert tick_texts
    assert all('e' not in t and '×' not in t for t in tick_texts)
    plt.close(fig)


def test_focused_pair_ticks_scaled_to_unit():
    """Mbp-scale focused plots show scaled tick values, not raw bp."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    rec = PafRecord(
        'q1',
        3_000_000,
        0,
        2_500_000,
        '+',
        't1',
        4_000_000,
        0,
        2_500_000,
        2_000_000,
        2_500_000,
        255,
    )
    plotter = DotPlotter(PafAlignment([rec]))
    fig = plotter.plot(query_names=['q1'], target_names=['t1'])
    ax = fig.axes[0]
    assert ax.get_xlabel() == 't1 (Mbp)'
    assert ax.get_ylabel() == 'q1 (Mbp)'
    fig.canvas.draw()
    tick_texts = [t.get_text() for t in ax.get_xticklabels() if t.get_text()]
    # All tick values fit in the scaled range (no 7-digit raw bp text).
    assert all(len(t.replace('.', '').replace('-', '')) <= 4 for t in tick_texts)
    plt.close(fig)


def test_grid_keeps_rotated_titles(dotplot_index):
    """Multi-panel grids keep the rotated top titles and raw-bp ticks."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot(query_names=['seq1', 'seq2'], target_names=['seq1', 'seq2'])
    top_titles = [ax.get_title() for ax in fig.axes[:2]]
    assert top_titles == ['seq1', 'seq2']
    plt.close(fig)


def test_grid_shared_units_and_angled_x_ticks():
    """Grids share one bp unit across contigs and angle the x tick labels."""
    from rusty_dot.paf_io import PafAlignment, PafRecord

    records = [
        PafRecord(
            'q1',
            3_000_000,
            0,
            2_000_000,
            '+',
            't1',
            4_000_000,
            0,
            2_000_000,
            1_500_000,
            2_000_000,
            255,
        ),
        PafRecord(
            'q2', 40_000, 0, 30_000, '+', 't2', 50_000, 0, 30_000, 25_000, 30_000, 255
        ),
    ]
    plotter = DotPlotter(PafAlignment(records))
    fig = plotter.plot(query_names=['q1', 'q2'], target_names=['t1', 't2'])
    fig.canvas.draw()
    bottom_left, bottom_right = fig.axes[2], fig.axes[3]
    # The short q2/t2 contigs use the grid-wide Mbp unit, not their own Kbp.
    texts = [t.get_text() for t in bottom_right.get_xticklabels() if t.get_text()]
    assert texts
    assert all(float(t) < 1 for t in texts if t not in ('0',))
    # x tick labels are angled to avoid overlap.
    assert all(t.get_rotation() == 45 for t in bottom_left.get_xticklabels())
    # The shared unit is announced once per axis at figure level.
    fig_texts = [t.get_text() for t in fig.texts]
    assert fig_texts.count('Position (Mbp)') == 2
    plt.close(fig)


# ---------------------------------------------- focused-view canvas margins


def _extreme_pair_plotter():
    from rusty_dot.paf_io import PafAlignment, PafRecord

    rec = PafRecord(
        'tiny_query_contig',
        30_000,
        0,
        30_000,
        '+',
        'very_long_target_chromosome',
        3_000_000,
        0,
        30_000,
        28_000,
        30_000,
        255,
    )
    return DotPlotter(PafAlignment([rec]))


def test_focused_extreme_ratio_title_clear_of_panel():
    """1:100 pair: the title must not overlap the plot area."""
    plotter = _extreme_pair_plotter()
    fig = plotter.plot(
        query_names=['tiny_query_contig'],
        target_names=['very_long_target_chromosome'],
        title='tiny_query_contig vs very_long_target_chromosome',
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_bbox = fig.axes[0].get_window_extent(renderer)
    title_bbox = fig._suptitle.get_window_extent(renderer)
    assert not ax_bbox.overlaps(title_bbox)
    plt.close(fig)


def test_focused_extreme_ratio_labels_inside_canvas():
    """1:100 pair: axis name labels fit fully inside the figure canvas."""
    plotter = _extreme_pair_plotter()
    fig = plotter.plot(
        query_names=['tiny_query_contig'],
        target_names=['very_long_target_chromosome'],
        title='a title',
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_bbox = fig.bbox
    ax = fig.axes[0]
    for label in (ax.yaxis.label, ax.xaxis.label):
        bbox = label.get_window_extent(renderer)
        assert fig_bbox.contains(bbox.x0, bbox.y0)
        assert fig_bbox.contains(bbox.x1, bbox.y1)
    plt.close(fig)


def test_focused_extreme_ratio_panel_aspect_preserved():
    """Margins must not distort the panel's bp-per-inch parity."""
    plotter = _extreme_pair_plotter()
    fig = plotter.plot(
        query_names=['tiny_query_contig'],
        target_names=['very_long_target_chromosome'],
        title='a title',
    )
    fig.canvas.draw()
    ax = fig.axes[0]
    bbox = ax.get_window_extent(fig.canvas.get_renderer())
    x_bp_per_px = 3_000_000 / bbox.width
    y_bp_per_px = 30_000 / bbox.height
    assert abs(x_bp_per_px - y_bp_per_px) / x_bp_per_px < 0.01
    plt.close(fig)


def test_focused_thin_axis_single_end_tick():
    """An ultra-thin axis shows one end tick (the sequence length)."""
    plotter = _extreme_pair_plotter()
    fig = plotter.plot(
        query_names=['tiny_query_contig'],
        target_names=['very_long_target_chromosome'],
    )
    fig.canvas.draw()
    ax = fig.axes[0]
    texts = [t.get_text() for t in ax.get_yticklabels() if t.get_text()]
    assert texts == ['30']
    plt.close(fig)


# ----------------------------------------------- multi-row contig labels


def _grid_plotter(lengths, prefix='query_contig_with_a_long_name_'):
    import random

    from rusty_dot import SequenceIndex

    rng = random.Random(7)
    idx = SequenceIndex(k=11)
    names = []
    for i, length in enumerate(lengths, 1):
        name = f'{prefix}{i:03d}'
        names.append(name)
        idx.add_sequence(name, ''.join(rng.choice('ACGT') for _ in range(length)))
    idx.add_sequence('target_contig', ''.join(rng.choice('ACGT') for _ in range(3000)))
    return DotPlotter(idx), names


def test_multi_row_contig_labels_are_angled():
    """Vertical row labels collide when contigs differ in length.

    A 90-degree label is as tall as the name is long, so a short contig's
    row cannot contain it and neighbouring labels smear together.
    """
    pl, names = _grid_plotter([3000, 2500, 2000])
    fig = pl.plot(query_names=names, target_names=['target_contig'])
    left_col = [ax for ax in fig.axes if ax.get_ylabel()]
    assert len(left_col) == 3
    for ax in left_col:
        rot = ax.yaxis.label.get_rotation()
        assert 0 < rot < 90, f'row label should be angled, got {rot}'
    plt.close(fig)


def test_single_row_keeps_the_conventional_vertical_label():
    """One row cannot collide with anything, so leave it alone."""
    pl, names = _grid_plotter([3000])
    fig = pl.plot(query_names=names, target_names=['target_contig'])
    ax = next(ax for ax in fig.axes if ax.get_ylabel())
    assert ax.yaxis.label.get_rotation() == 90
    plt.close(fig)


def test_row_labels_shrink_to_fit_a_short_row():
    """A thin row gets a smaller label rather than one that overruns it."""
    pl, names = _grid_plotter([30000, 900])
    fig = pl.plot(query_names=names, target_names=['target_contig'])
    labelled = [ax for ax in fig.axes if ax.get_ylabel()]
    tall, short = labelled[0], labelled[1]
    assert short.yaxis.label.get_fontsize() < tall.yaxis.label.get_fontsize(), (
        'the short row should use a smaller label'
    )
    plt.close(fig)


def test_row_label_never_exceeds_its_row_height():
    """The guarantee the rotation exists to provide: no overlap."""
    import math

    pl, names = _grid_plotter([30000, 1200, 900, 800])
    fig = pl.plot(query_names=names, target_names=['target_contig'])
    fig_h = fig.get_size_inches()[1]
    for ax in [ax for ax in fig.axes if ax.get_ylabel()]:
        label = ax.yaxis.label
        extent_in = (
            len(ax.get_ylabel())
            * label.get_fontsize()
            * 0.62
            / 72.0
            * math.sin(math.radians(label.get_rotation()))
        )
        row_in = ax.get_position().height * fig_h
        assert extent_in <= row_in * 1.05, (
            f'label {ax.get_ylabel()!r} spans {extent_in:.2f}in of a {row_in:.2f}in row'
        )
    plt.close(fig)
