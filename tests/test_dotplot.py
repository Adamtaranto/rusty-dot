"""Tests for the dotplot visualization module."""

import os

import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter


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
    assert len(ax.lines) > 0

    # For every RC match line, the x-data must be decreasing (anti-diagonal)
    for line in ax.lines:
        xdata = line.get_xdata()
        if len(xdata) == 2:
            # anti-diagonal: x goes from t_end to t_start (decreasing)
            assert xdata[0] >= xdata[1], (
                f'Expected anti-diagonal line (x decreasing) for RC match, '
                f'got x={xdata}'
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

    fwd_lines = [
        line for line in ax.lines if list(line.get_xdata()) == sorted(line.get_xdata())
    ]
    assert len(fwd_lines) > 0, 'Expected at least one forward (diagonal) match line'


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
    n_lines_all = len(ax_all.lines)
    plt.close(fig)

    fig, ax_filtered = plt.subplots()
    # A very large min_length should filter everything out
    plotter._plot_panel(ax_filtered, 'q', 't', min_length=10000)
    n_lines_filtered = len(ax_filtered.lines)
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

    # Two PAF records → two lines drawn
    assert len(ax.lines) == 2


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
    n_kmer = len(ax_kmer.lines)
    plt.close(fig_kmer)

    # Suppress the expected warning to keep test output clean.
    fig_warn, ax_warn = plt.subplots()
    logging.disable(logging.WARNING)
    try:
        plotter._plot_panel(ax_warn, 'seq1', 'seq2', color_by_identity=True)
    finally:
        logging.disable(logging.NOTSET)
    n_fallback = len(ax_warn.lines)
    plt.close(fig_warn)

    # fallback should produce the same number of lines as the normal k-mer plot
    assert n_fallback == n_kmer


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
    n_all = len(ax_all.lines)
    plt.close(fig)

    fig, ax_filtered = plt.subplots()
    plotter._plot_panel(
        ax_filtered, 'seq1', 'seq2', color_by_identity=True, min_length=8
    )
    n_filtered = len(ax_filtered.lines)
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

    # One '+' alignment should produce one line
    assert len(ax.lines) == 1


def test_dotplotter_paf_as_index_rc_alignment():
    """DotPlotter(PafAlignment) renders RC alignment as anti-diagonal line."""
    import matplotlib

    matplotlib.use('Agg')

    aln = _make_paf_alignment_index()
    plotter = DotPlotter(aln)

    fig, ax = plt.subplots()
    plotter._plot_panel(ax, 'contigA', 'contigC')
    plt.close(fig)

    assert len(ax.lines) == 1
    xdata = ax.lines[0].get_xdata()
    # RC strand: x goes from target_end to target_start (decreasing)
    assert xdata[0] > xdata[1]


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
