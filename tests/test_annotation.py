"""Tests for the GffAnnotation class."""

import os

import matplotlib.figure
import matplotlib.pyplot as plt
import pytest

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.annotation import (
    SUPPORTED_PALETTES,
    GffAnnotation,
    GffFeature,
    _parse_gff,
)
from rusty_dot.dotplot import DotPlotter

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_gff_file(tmp_path, lines):
    """Write a GFF file and return its path."""
    path = tmp_path / 'test.gff'
    path.write_text('\n'.join(lines) + '\n')
    return str(path)


@pytest.fixture
def simple_gff(tmp_path):
    """A small GFF3 file with two feature types."""
    lines = [
        '##gff-version 3',
        'seq1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1',
        'seq1\t.\tCDS\t10\t90\t.\t+\t0\tID=cds1',
        'seq2\t.\tgene\t200\t400\t.\t-\t.\tID=gene2',
        'seq2\t.\trepeat_region\t500\t600\t0.5\t.\t.\tID=rep1',
    ]
    return _make_gff_file(tmp_path, lines)


@pytest.fixture
def simple_annotation(simple_gff):
    """GffAnnotation loaded from simple_gff."""
    return GffAnnotation.from_file(simple_gff)


@pytest.fixture
def dotplot_index():
    """Minimal index for DotPlotter tests."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGTACGTACGTACGTACGT')
    idx.add_sequence('seq2', 'TACGTACGTACGTACGTACG')
    return idx


# ---------------------------------------------------------------------------
# _parse_gff
# ---------------------------------------------------------------------------


def test_parse_gff_yields_features(simple_gff):
    records = list(_parse_gff(simple_gff))
    assert len(records) == 4


def test_parse_gff_skips_comments(simple_gff):
    records = list(_parse_gff(simple_gff))
    for r in records:
        assert not r.seqname.startswith('#')


def test_parse_gff_coordinate_conversion(simple_gff):
    records = list(_parse_gff(simple_gff))
    gene1 = next(r for r in records if r.attributes == 'ID=gene1')
    # GFF 1-based → 0-based: start becomes 0, end stays 100
    assert gene1.start == 0
    assert gene1.end == 100


def test_parse_gff_score_parsing(simple_gff):
    records = list(_parse_gff(simple_gff))
    rep1 = next(r for r in records if 'rep1' in r.attributes)
    assert rep1.score == pytest.approx(0.5)


def test_parse_gff_dot_score_is_none(simple_gff):
    records = list(_parse_gff(simple_gff))
    gene1 = next(r for r in records if r.attributes == 'ID=gene1')
    assert gene1.score is None


def test_parse_gff_dot_frame_is_none(simple_gff):
    records = list(_parse_gff(simple_gff))
    gene1 = next(r for r in records if r.attributes == 'ID=gene1')
    assert gene1.frame is None


def test_parse_gff_integer_frame(simple_gff):
    records = list(_parse_gff(simple_gff))
    cds1 = next(r for r in records if 'cds1' in r.attributes)
    assert cds1.frame == 0


def test_parse_gff_skips_short_lines(tmp_path):
    path = _make_gff_file(tmp_path, ['seq1\tonly\tthree'])
    records = list(_parse_gff(path))
    assert records == []


def test_parse_gff_empty_file(tmp_path):
    path = _make_gff_file(tmp_path, [])
    records = list(_parse_gff(path))
    assert records == []


# ---------------------------------------------------------------------------
# GffAnnotation.from_file
# ---------------------------------------------------------------------------


def test_from_file_returns_annotation(simple_gff):
    ann = GffAnnotation.from_file(simple_gff)
    assert isinstance(ann, GffAnnotation)


def test_from_file_record_count(simple_annotation):
    assert len(simple_annotation) == 4


# ---------------------------------------------------------------------------
# GffAnnotation.feature_types
# ---------------------------------------------------------------------------


def test_feature_types_sorted(simple_annotation):
    types = simple_annotation.feature_types()
    assert types == sorted(types)


def test_feature_types_unique(simple_annotation):
    types = simple_annotation.feature_types()
    assert len(types) == len(set(types))


def test_feature_types_content(simple_annotation):
    types = simple_annotation.feature_types()
    assert 'gene' in types
    assert 'CDS' in types
    assert 'repeat_region' in types


# ---------------------------------------------------------------------------
# GffAnnotation.sequence_names
# ---------------------------------------------------------------------------


def test_sequence_names_sorted(simple_annotation):
    names = simple_annotation.sequence_names()
    assert names == sorted(names)


def test_sequence_names_content(simple_annotation):
    names = simple_annotation.sequence_names()
    assert 'seq1' in names
    assert 'seq2' in names


# ---------------------------------------------------------------------------
# Color assignment
# ---------------------------------------------------------------------------


def test_colors_assigned_for_all_types(simple_annotation):
    for ft in simple_annotation.feature_types():
        color = simple_annotation.get_color(ft)
        assert color.startswith('#')


def test_unknown_type_returns_fallback(simple_annotation):
    assert simple_annotation.get_color('nonexistent') == '#888888'


def test_custom_color_in_constructor(simple_gff):
    ann = GffAnnotation.from_file(simple_gff, colors={'gene': '#ff0000'})
    assert ann.get_color('gene') == '#ff0000'


def test_set_colors_overrides(simple_annotation):
    simple_annotation.set_colors({'CDS': '#00ff00'})
    assert simple_annotation.get_color('CDS') == '#00ff00'


def test_set_colors_does_not_affect_other_types(simple_annotation):
    original_gene_color = simple_annotation.get_color('gene')
    simple_annotation.set_colors({'CDS': '#00ff00'})
    assert simple_annotation.get_color('gene') == original_gene_color


# ---------------------------------------------------------------------------
# Palette selection
# ---------------------------------------------------------------------------


def test_default_palette_is_tab10(simple_gff):
    ann = GffAnnotation.from_file(simple_gff)
    assert ann._palette == 'tab10'


def test_valid_palette_accepted(simple_gff):
    for palette in SUPPORTED_PALETTES:
        ann = GffAnnotation.from_file(simple_gff, palette=palette)
        assert ann._palette == palette


def test_invalid_palette_raises(simple_gff):
    with pytest.raises(ValueError, match='Unknown palette'):
        GffAnnotation.from_file(simple_gff, palette='notapalette')


# ---------------------------------------------------------------------------
# keep_feature_types
# ---------------------------------------------------------------------------


def test_keep_feature_types_filters(simple_annotation):
    kept = simple_annotation.keep_feature_types(['gene'])
    assert kept.feature_types() == ['gene']


def test_keep_feature_types_returns_new_instance(simple_annotation):
    kept = simple_annotation.keep_feature_types(['gene'])
    assert kept is not simple_annotation


def test_keep_feature_types_preserves_colors(simple_annotation):
    original_color = simple_annotation.get_color('gene')
    kept = simple_annotation.keep_feature_types(['gene'])
    assert kept.get_color('gene') == original_color


def test_keep_feature_types_empty_list(simple_annotation):
    kept = simple_annotation.keep_feature_types([])
    assert len(kept) == 0


def test_keep_feature_types_multiple(simple_annotation):
    kept = simple_annotation.keep_feature_types(['gene', 'CDS'])
    assert set(kept.feature_types()) == {'gene', 'CDS'}


# ---------------------------------------------------------------------------
# drop_feature_types
# ---------------------------------------------------------------------------


def test_drop_feature_types_removes(simple_annotation):
    dropped = simple_annotation.drop_feature_types(['gene'])
    assert 'gene' not in dropped.feature_types()


def test_drop_feature_types_returns_new_instance(simple_annotation):
    dropped = simple_annotation.drop_feature_types(['gene'])
    assert dropped is not simple_annotation


def test_drop_feature_types_keeps_others(simple_annotation):
    dropped = simple_annotation.drop_feature_types(['gene'])
    assert 'CDS' in dropped.feature_types()
    assert 'repeat_region' in dropped.feature_types()


def test_drop_feature_types_unknown_type_noop(simple_annotation):
    original_count = len(simple_annotation)
    dropped = simple_annotation.drop_feature_types(['nonexistent'])
    assert len(dropped) == original_count


# ---------------------------------------------------------------------------
# filter_by_sequence
# ---------------------------------------------------------------------------


def test_filter_by_sequence(simple_annotation):
    filtered = simple_annotation.filter_by_sequence(['seq1'])
    seqs = filtered.sequence_names()
    assert seqs == ['seq1']


def test_filter_by_sequence_returns_new_instance(simple_annotation):
    filtered = simple_annotation.filter_by_sequence(['seq1'])
    assert filtered is not simple_annotation


def test_filter_by_sequence_count(simple_annotation):
    # seq1 has 2 features (gene + CDS)
    filtered = simple_annotation.filter_by_sequence(['seq1'])
    assert len(filtered) == 2


def test_filter_by_sequence_multiple(simple_annotation):
    filtered = simple_annotation.filter_by_sequence(['seq1', 'seq2'])
    assert len(filtered) == 4


def test_filter_by_sequence_missing_seq(simple_annotation):
    filtered = simple_annotation.filter_by_sequence(['nonexistent'])
    assert len(filtered) == 0


# ---------------------------------------------------------------------------
# get_features_for_sequence
# ---------------------------------------------------------------------------


def test_get_features_for_sequence(simple_annotation):
    feats = simple_annotation.get_features_for_sequence('seq1')
    assert len(feats) == 2
    assert all(f.seqname == 'seq1' for f in feats)


def test_get_features_for_sequence_missing(simple_annotation):
    feats = simple_annotation.get_features_for_sequence('nonexistent')
    assert feats == []


# ---------------------------------------------------------------------------
# records property
# ---------------------------------------------------------------------------


def test_records_returns_list(simple_annotation):
    records = simple_annotation.records
    assert isinstance(records, list)
    assert len(records) == 4


def test_records_returns_copy(simple_annotation):
    records = simple_annotation.records
    records.clear()
    # Internal state should be unchanged
    assert len(simple_annotation) == 4


# ---------------------------------------------------------------------------
# __len__ and __repr__
# ---------------------------------------------------------------------------


def test_len(simple_annotation):
    assert len(simple_annotation) == 4


def test_repr(simple_annotation):
    r = repr(simple_annotation)
    assert 'GffAnnotation' in r
    assert '4' in r


# ---------------------------------------------------------------------------
# DotPlotter integration: plot() with annotation
# ---------------------------------------------------------------------------


def test_plot_annotation_grid(dotplot_index, simple_annotation, tmp_path):
    """plot() with annotation parameter produces a valid output file."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'annotated_grid.png')
    fig = plotter.plot(
        query_names=['seq1', 'seq2'],
        target_names=['seq1', 'seq2'],
        output_path=output,
        annotation=simple_annotation,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_annotation_squares_drawn(dotplot_index, simple_annotation):
    """Annotation squares are added to self-vs-self diagonal panels."""
    plotter = DotPlotter(dotplot_index)

    fig, axes = plt.subplots(2, 2, squeeze=False)
    query_names = ['seq1', 'seq2']
    target_names = ['seq1', 'seq2']
    for row_idx, q in enumerate(query_names):
        for col_idx, t in enumerate(target_names):
            ax = axes[row_idx][col_idx]
            plotter._plot_panel(ax, q, t)
            if q == t:
                plotter._draw_annotation_squares(ax, q, simple_annotation)

    plt.close(fig)
    # seq1 has 2 features → 2 patches on diagonal panel (0,0)
    n_patches_diag = len(axes[0][0].patches)
    n_patches_off = len(axes[0][1].patches)
    assert n_patches_diag == 2  # two seq1 features
    assert n_patches_off == 0   # off-diagonal: no patches


def test_plot_warns_when_annotation_contains_unknown_sequences(dotplot_index, caplog):
    """A warning is emitted when annotation has features for sequences not in the index."""
    import logging

    # Build annotation with a sequence that is NOT in dotplot_index
    ann = GffAnnotation(
        [GffFeature('unknown_seq', '.', 'gene', 0, 100, None, '+', None, '')],
    )
    plotter = DotPlotter(dotplot_index)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        fig = plotter.plot(annotation=ann)
    plt.close(fig)
    assert any('unknown_seq' in msg for msg in caplog.messages)


def test_plot_grid_column_labels_on_top(dotplot_index):
    """In a grid plot, column (target) names appear as titles on the top row."""
    plotter = DotPlotter(dotplot_index)
    query_names = ['seq1', 'seq2']
    target_names = ['seq1', 'seq2']
    fig = plotter.plot(
        query_names=query_names,
        target_names=target_names,
    )
    # The first len(target_names) axes in fig.axes correspond to the top row.
    # Each should have a title set to the target (column) name.
    titles = [fig.axes[i].get_title() for i in range(len(target_names))]
    assert titles == target_names, f'Expected titles {target_names}, got {titles}'
    plt.close(fig)


def test_plot_grid_returns_figure(dotplot_index):
    """plot() with annotation returns a matplotlib Figure."""
    ann = GffAnnotation(
        [GffFeature('seq1', '.', 'gene', 0, 10, None, '+', None, '')],
    )
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot(annotation=ann)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


# ---------------------------------------------------------------------------
# DotPlotter integration: plot_single() with annotation
# ---------------------------------------------------------------------------


def test_plot_single_with_annotation(dotplot_index, simple_annotation, tmp_path):
    """plot_single() with annotation produces a valid output file."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'annotated_single.png')
    fig = plotter.plot_single(
        'seq1',
        'seq2',
        output_path=output,
        annotation=simple_annotation,
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single_annotation_tracks_drawn(dotplot_index, simple_annotation):
    """When annotation is provided, annotation track axes are created."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_single('seq1', 'seq2', annotation=simple_annotation)
    # Figure should have more than 1 axes (main + tracks + corner)
    assert len(fig.axes) > 1
    plt.close(fig)


def test_plot_single_no_annotation_single_ax(dotplot_index):
    """Without annotation, plot_single produces a single-axes figure."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_single('seq1', 'seq2')
    assert len(fig.axes) == 1
    plt.close(fig)


def test_plot_single_self_vs_self_with_annotation(dotplot_index, simple_annotation, tmp_path):
    """plot_single() works for self-vs-self with annotation tracks."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'self_annotated.png')
    fig = plotter.plot_single(
        'seq1', 'seq1', output_path=output, annotation=simple_annotation
    )
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single_annotation_warns_missing(dotplot_index, caplog, simple_gff):
    """plot_single() warns when annotation has features for unknown sequences."""
    import logging

    ann = GffAnnotation(
        [GffFeature('missing_seq', '.', 'gene', 0, 100, None, '+', None, '')],
    )
    plotter = DotPlotter(dotplot_index)
    with caplog.at_level(logging.WARNING, logger='rusty_dot.dotplot'):
        fig = plotter.plot_single('seq1', 'seq2', annotation=ann)
    plt.close(fig)
    assert any('missing_seq' in msg for msg in caplog.messages)


# ---------------------------------------------------------------------------
# DotPlotter.plot_annotation_legend
# ---------------------------------------------------------------------------


def test_plot_annotation_legend_returns_figure(dotplot_index, simple_annotation):
    """plot_annotation_legend() returns a matplotlib Figure."""
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_annotation_legend(simple_annotation)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)


def test_plot_annotation_legend_saves_file(dotplot_index, simple_annotation, tmp_path):
    """plot_annotation_legend() saves a file when output_path is provided."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / 'legend.png')
    fig = plotter.plot_annotation_legend(simple_annotation, output_path=output)
    plt.close(fig)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_annotation_legend_no_file_without_output(
    dotplot_index, tmp_path, monkeypatch
):
    """plot_annotation_legend() with output_path=None creates no file."""
    monkeypatch.chdir(tmp_path)
    # Build annotation directly without writing to tmp_path
    ann = GffAnnotation(
        [GffFeature('seq1', '.', 'gene', 0, 100, None, '+', None, '')],
    )
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_annotation_legend(ann)
    plt.close(fig)
    assert list(tmp_path.iterdir()) == []


def test_plot_annotation_legend_empty_annotation(dotplot_index):
    """plot_annotation_legend() works with an empty annotation."""
    ann = GffAnnotation([])
    plotter = DotPlotter(dotplot_index)
    fig = plotter.plot_annotation_legend(ann)
    assert isinstance(fig, matplotlib.figure.Figure)
    plt.close(fig)
