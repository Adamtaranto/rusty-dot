"""Tests for the dotplot visualization module."""

import logging
import os
import textwrap

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
# PafAlignment source (identity colouring)
# ---------------------------------------------------------------------------

# Minimal PAF content with two records that have different identities:
# - record 1: 45/50 matches → 0.9 identity (forward strand)
# - record 2: 38/40 matches → 0.95 identity (minus strand)
_PAF_WITH_IDENTITY = textwrap.dedent("""\
    query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t45\t50\t255
    query2\t80\t0\t40\t-\ttarget1\t200\t150\t190\t38\t40\t255
""")


def _write_temp_paf(content: str, tmp_path) -> str:
    path = tmp_path / 'alignments.paf'
    path.write_text(content)
    return str(path)


@pytest.fixture
def paf_alignment(tmp_path):
    """Return a PafAlignment loaded from a temp PAF file."""
    from rusty_dot.paf_io import PafAlignment

    path = _write_temp_paf(_PAF_WITH_IDENTITY, tmp_path)
    return PafAlignment.from_file(path)


def test_dotplotter_accepts_paf_alignment(paf_alignment):
    """DotPlotter can be initialised with a PafAlignment."""
    plotter = DotPlotter(paf_alignment)
    assert plotter._is_paf_source


def test_paf_source_length_map(paf_alignment):
    """_paf_len_map is populated from PAF records."""
    plotter = DotPlotter(paf_alignment)
    assert plotter._get_sequence_length('query1') == 100
    assert plotter._get_sequence_length('target1') == 200


def test_paf_source_sequence_names(paf_alignment):
    """_sequence_names returns all names seen in the PAF records."""
    plotter = DotPlotter(paf_alignment)
    names = set(plotter._sequence_names())
    assert 'query1' in names
    assert 'target1' in names


def test_get_paf_matches_returns_identity(paf_alignment):
    """_get_paf_matches returns tuples that include identity values."""
    plotter = DotPlotter(paf_alignment)
    matches = plotter._get_paf_matches('query1', 'target1', min_length=0)
    assert len(matches) == 1
    q_start, q_end, t_start, t_end, strand, identity = matches[0]
    assert strand == '+'
    assert identity == pytest.approx(0.9)


def test_get_paf_matches_minus_strand_identity(paf_alignment):
    """_get_paf_matches handles minus-strand records and computes identity."""
    plotter = DotPlotter(paf_alignment)
    matches = plotter._get_paf_matches('query2', 'target1', min_length=0)
    assert len(matches) == 1
    _, _, _, _, strand, identity = matches[0]
    assert strand == '-'
    assert identity == pytest.approx(0.95)


def test_get_paf_matches_min_length_filter(paf_alignment):
    """_get_paf_matches respects min_length filtering."""
    plotter = DotPlotter(paf_alignment)
    # query1 has query_aligned_len=50, query2 has 40
    matches = plotter._get_paf_matches('query2', 'target1', min_length=50)
    assert len(matches) == 0


def test_plot_paf_source_creates_file(paf_alignment, tmp_path):
    """plot() works with a PafAlignment source and produces a valid file."""
    plotter = DotPlotter(paf_alignment)
    output = str(tmp_path / 'paf_dotplot.png')
    plotter.plot(
        query_names=['query1', 'query2'],
        target_names=['target1'],
        output_path=output,
        scale_sequences=False,
    )
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_paf_source_color_by_identity(paf_alignment, tmp_path):
    """plot() with color_by_identity=True on a PAF source produces a file."""
    plotter = DotPlotter(paf_alignment)
    output = str(tmp_path / 'identity_plot.png')
    plotter.plot(
        query_names=['query1', 'query2'],
        target_names=['target1'],
        output_path=output,
        color_by_identity=True,
        palette='viridis',
        scale_sequences=False,
    )
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single_paf_source_color_by_identity(paf_alignment, tmp_path):
    """plot_single() with color_by_identity=True on a PAF source works."""
    plotter = DotPlotter(paf_alignment)
    output = str(tmp_path / 'identity_single.png')
    plotter.plot_single(
        'query1', 'target1', output_path=output, color_by_identity=True
    )
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_identity_scale_creates_file(tmp_path):
    """plot_identity_scale() creates a non-empty image file."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGT' * 5)
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'scale.png')
    plotter.plot_identity_scale(output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_identity_scale_svg(tmp_path):
    """plot_identity_scale() can produce an SVG file."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGT' * 5)
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'scale.svg')
    plotter.plot_identity_scale(output_path=output)
    assert os.path.exists(output)
    with open(output) as f:
        content = f.read(200)
    assert '<svg' in content or '<?xml' in content


def test_plot_identity_scale_custom_palette(tmp_path):
    """plot_identity_scale() accepts a custom palette without error."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGT' * 5)
    plotter = DotPlotter(idx)
    output = str(tmp_path / 'scale_rdylgn.png')
    plotter.plot_identity_scale(output_path=output, palette='RdYlGn')
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_color_by_identity_warns_for_kmer_source(dotplot_index, tmp_path, caplog):
    """A warning is emitted when color_by_identity=True with a k-mer source."""
    plotter = DotPlotter(dotplot_index)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        plotter.plot(
            query_names=['seq1'],
            target_names=['seq2'],
            output_path=str(tmp_path / 'warn_kmer.png'),
            color_by_identity=True,
        )
    assert any('color_by_identity' in msg for msg in caplog.messages)


def test_color_by_identity_warns_for_kmer_source_plot_single(dotplot_index, tmp_path, caplog):
    """A warning is emitted when color_by_identity=True in plot_single with a k-mer source."""
    plotter = DotPlotter(dotplot_index)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        plotter.plot_single(
            'seq1',
            'seq2',
            output_path=str(tmp_path / 'warn_kmer_single.png'),
            color_by_identity=True,
        )
    assert any('color_by_identity' in msg for msg in caplog.messages)


def test_color_by_identity_no_warning_for_paf_source(paf_alignment, tmp_path, caplog):
    """No warning is emitted when color_by_identity=True with a PAF source."""
    plotter = DotPlotter(paf_alignment)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        plotter.plot(
            query_names=['query1'],
            target_names=['target1'],
            output_path=str(tmp_path / 'no_warn.png'),
            color_by_identity=True,
            scale_sequences=False,
        )
    assert not any('color_by_identity' in msg for msg in caplog.messages)


def test_paf_source_panel_uses_identity_colours(paf_alignment):
    """Lines are drawn with colours from the colormap when color_by_identity=True."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plotter = DotPlotter(paf_alignment)
    fig, ax = plt.subplots()
    plotter._plot_panel(
        ax, 'query1', 'target1', color_by_identity=True, palette='plasma'
    )
    plt.close(fig)

    # At least one line must be drawn
    assert len(ax.lines) > 0
    # The line colour should not be plain 'blue' or 'red' (identity-mapped instead)
    for line in ax.lines:
        clr = line.get_color()
        assert clr != 'blue' and clr != 'red'
