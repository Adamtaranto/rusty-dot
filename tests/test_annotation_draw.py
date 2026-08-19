"""Tests for GFF parser upgrades and annotation track drawing."""

import gzip
import json
import re

import matplotlib

matplotlib.use('Agg')

from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon
import matplotlib.pyplot as plt
import pytest

from rusty_dot import DotPlotter, SequenceIndex
from rusty_dot._annotation_draw import (
    annotation_legend_handles,
    assign_lanes,
    draw_track,
)
from rusty_dot.annotation import GffAnnotation, parse_attributes

GFF = '\n'.join(
    [
        '##gff-version 3',
        'c1\ttest\tgene\t51\t600\t.\t+\t.\tID=g1;Name=GeneA',
        'c1\ttest\tCDS\t51\t200\t.\t+\t0\tID=cds1;Parent=g1',
        'c1\ttest\tCDS\t351\t600\t.\t+\t0\tID=cds1;Parent=g1',
        'c1\ttest\tgene\t401\t900\t.\t-\t.\tID=g1b;Name=Gene%3BB',
        'c1\ttest\trepeat_region\t1001\t1500\t.\t.\t.\tID=rep1',
        'c2\ttest\tgene\t101\t800\t.\t-\t.\tID=g2',
    ]
)


@pytest.fixture(autouse=True)
def _close_figures():
    """Close every figure a test opened (matplotlib warns past 20 open)."""
    yield
    plt.close('all')


# ------------------------------------------------------------ parser


def test_parse_attributes_basic_and_percent_decoding():
    attrs = parse_attributes('ID=g1;Name=Gene%3BB;Note=a%2Cb;')
    assert attrs == {'ID': 'g1', 'Name': 'Gene;B', 'Note': 'a,b'}


def test_parse_attributes_gff2_fallback_and_empty():
    assert parse_attributes('gene_id "abc"; note') == {'gene_id': 'abc', 'note': ''}
    assert parse_attributes('') == {}


def test_feature_id_parent_name_properties():
    ann = GffAnnotation.from_text(GFF)
    gene = next(r for r in ann.records if r.feature_id == 'g1')
    assert gene.parent is None
    assert gene.name == 'GeneA'
    cds = next(r for r in ann.records if r.feature_type == 'CDS')
    assert cds.feature_id == 'cds1'
    assert cds.parent == 'g1'
    assert cds.name == 'cds1'  # Name falls back to ID
    escaped = next(r for r in ann.records if r.feature_id == 'g1b')
    assert escaped.name == 'Gene;B'


def test_parent_takes_first_of_comma_list():
    ann = GffAnnotation.from_text('c1\tt\texon\t1\t10\t.\t+\t.\tID=e;Parent=m1,m2')
    assert ann.records[0].parent == 'm1'


def test_from_bytes_plain_gzip_and_errors():
    plain = GffAnnotation.from_bytes(GFF.encode())
    assert len(plain) == 6
    zipped = GffAnnotation.from_bytes(gzip.compress(GFF.encode()))
    assert [r.feature_id for r in zipped.records] == [
        r.feature_id for r in plain.records
    ]
    with pytest.raises(ValueError, match='gzip'):
        GffAnnotation.from_bytes(b'\x1f\x8bnot-really-gzip')
    with pytest.raises(ValueError, match='UTF-8'):
        GffAnnotation.from_bytes(b'\xff\xfe\x00bad')


def test_fasta_section_terminates_parsing():
    text = GFF + '\n##FASTA\n>c1\nACGT\n'
    ann = GffAnnotation.from_text(text)
    assert len(ann) == 6  # nothing after ##FASTA is parsed


def test_get_features_for_sequence_uses_lazy_index():
    ann = GffAnnotation.from_text(GFF)
    first = ann.get_features_for_sequence('c1')
    second = ann.get_features_for_sequence('c1')
    assert [f.feature_id for f in first] == [f.feature_id for f in second]
    assert ann.get_features_for_sequence('missing') == []


def test_iter_groups_multi_part_and_singletons():
    ann = GffAnnotation.from_text(GFF)
    groups = ann.iter_groups('c1')
    keys = [key for key, _ in groups]
    cds_parts = next(parts for key, parts in groups if key == ('cds1', 'g1'))
    assert len(cds_parts) == 2
    assert [p.start for p in cds_parts] == sorted(p.start for p in cds_parts)
    # Genes without a Parent are singleton groups.
    singleton_sizes = [len(parts) for key, parts in groups if key != ('cds1', 'g1')]
    assert singleton_sizes == [1, 1, 1]
    assert ('cds1', 'g1') in keys


# ------------------------------------------------------------ lanes


def test_assign_lanes_disjoint_share_lane():
    assert assign_lanes([(0, 10), (20, 30), (40, 50)]) == [0, 0, 0]


def test_assign_lanes_overlaps_stack():
    lanes = assign_lanes([(0, 100), (50, 150), (60, 90)])
    assert lanes[0] == 0
    assert lanes[1] == 1
    assert lanes[2] == 2


def test_assign_lanes_touching_endpoints_share_lane():
    assert assign_lanes([(0, 10), (10, 20)]) == [0, 0]
    # ...unless a gap is required.
    assert assign_lanes([(0, 10), (10, 20)], gap=5) == [0, 1]


def test_assign_lanes_preserves_input_order():
    # Input deliberately unsorted; lanes are reported in input order.
    lanes = assign_lanes([(50, 150), (0, 100)])
    assert lanes == [1, 0]


# ------------------------------------------------------------ tracks


def _track_ax():
    _fig, ax = plt.subplots()
    ax.set_xlim(0, 2000)
    return ax


def test_draw_track_arrow_head_follows_strand():
    ann = GffAnnotation.from_text(
        'c1\tt\tgene\t101\t500\t.\t+\t.\tID=a\nc1\tt\tgene\t701\t1100\t.\t-\t.\tID=b\n'
    )
    ax = _track_ax()
    lanes = draw_track(ax, ann, 'c1', 2000, orientation='x')
    assert lanes == 1  # disjoint features share a lane
    polys = [p for p in ax.patches if isinstance(p, Polygon)]
    assert len(polys) == 2
    fwd, rev = polys
    # '+' head vertex (unique x at lane centre y) sits at the feature end.
    fwd_xs = fwd.get_xy()[:, 0]
    assert fwd_xs.max() == 500
    assert (fwd_xs == 500).sum() == 1  # single tip vertex
    rev_xs = rev.get_xy()[:, 0]
    assert rev_xs.min() == 700
    assert (rev_xs == 700).sum() == 1


def test_draw_track_reverse_mirrors_and_flips():
    ann = GffAnnotation.from_text('c1\tt\tgene\t101\t500\t.\t+\t.\tID=a')
    ax = _track_ax()
    draw_track(ax, ann, 'c1', 2000, orientation='x', reverse=True)
    poly = [p for p in ax.patches if isinstance(p, Polygon)][0]
    xs = poly.get_xy()[:, 0]
    # Interval mirrors to [1500, 1900]; '+' flips to '-' so the tip is at
    # the mirrored start.
    assert xs.min() == 1500
    assert xs.max() == 1900
    assert (xs == 1500).sum() == 1


def test_draw_track_y_orientation_head_at_start_for_minus():
    ann = GffAnnotation.from_text('c2\tt\tgene\t101\t800\t.\t-\t.\tID=g2')
    _fig, ax = plt.subplots()
    ax.set_ylim(0, 2320)
    ax.invert_yaxis()
    draw_track(ax, ann, 'c2', 2320, orientation='y')
    poly = [p for p in ax.patches if isinstance(p, Polygon)][0]
    ys = poly.get_xy()[:, 1]
    assert ys.min() == 100
    assert (ys == 100).sum() == 1  # '-' strand: tip at the start coordinate


def test_draw_track_unstranded_rounded_and_connectors():
    ann = GffAnnotation.from_text(
        'c1\tt\trepeat_region\t201\t900\t.\t.\t.\tID=r1\n'
        'c1\tt\tCDS\t1001\t1200\t.\t+\t0\tID=c;Parent=m\n'
        'c1\tt\tCDS\t1501\t1800\t.\t+\t0\tID=c;Parent=m\n'
    )
    ax = _track_ax()
    draw_track(ax, ann, 'c1', 2000, orientation='x')
    assert any(isinstance(p, FancyBboxPatch) for p in ax.patches)
    connectors = [ln for ln in ax.lines if isinstance(ln, Line2D)]
    assert len(connectors) == 1
    xs = connectors[0].get_xdata()
    assert list(xs) == sorted(xs)  # midpoints in ascending order


def test_draw_track_empty_sequence_hides_axis():
    ann = GffAnnotation.from_text('c1\tt\tgene\t1\t10\t.\t+\t.\tID=a')
    ax = _track_ax()
    assert draw_track(ax, ann, 'other', 2000, orientation='x') == 0
    assert not ax.axison


def test_annotation_legend_handles():
    ann = GffAnnotation.from_text(GFF)
    handles = annotation_legend_handles(ann)
    assert [h.get_label() for h in handles] == ann.feature_types()
    subset = annotation_legend_handles(ann, ['gene'])
    assert [h.get_label() for h in subset] == ['gene']


# ------------------------------------------------------------ plot()


def _plotter():
    seq = 'ACGTTGCAAGGCCTTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCG' * 30
    idx = SequenceIndex(k=11)
    idx.add_sequence('c1', seq)
    idx.add_sequence('c2', seq[::-1] + seq[:400])
    return DotPlotter(idx)


def test_plot_single_pair_with_tracks_adds_axes_and_legend():
    pl = _plotter()
    ann = GffAnnotation.from_text(GFF)
    fig = pl.plot(
        query_names=['c2'],
        target_names=['c1'],
        annotation_query=ann,
        annotation_target=ann,
        annotation_tracks=True,
    )
    assert len(fig.axes) == 3  # main + y-track + x-track
    assert len(fig.legends) == 1
    plt.close(fig)


def test_plot_grid_ignores_tracks_and_draws_diagonal_behind():
    pl = _plotter()
    ann = GffAnnotation.from_text(GFF)
    fig = pl.plot(annotation=ann, annotation_tracks=True)
    assert len(fig.axes) == 4  # 2x2 grid, no track axes
    diag = fig.axes[0]
    coll = [c for c in diag.collections if c.get_gid() is None and c.get_zorder() < 1]
    # Squares are behind alignment segments.
    assert any(c.get_zorder() == 0.5 for c in diag.collections)
    line_z = [c.get_zorder() for c in diag.collections if c.get_zorder() >= 1]
    assert all(z > 0.5 for z in line_z)
    assert coll is not None
    plt.close(fig)


def test_diagonal_squares_mirror_on_reversed_contig():
    pl = _plotter()
    ann = GffAnnotation.from_text('c1\tt\tgene\t101\t300\t.\t+\t.\tID=a')
    seq_len = pl.index.get_sequence_length('c1')
    fig = pl.plot(
        query_names=['c1'],
        target_names=['c1'],
        annotation=ann,
        reverse_contigs={'c1'},
    )
    ax = fig.axes[0]
    patch_colls = [c for c in ax.collections if c.get_zorder() == 0.5]
    assert patch_colls
    rect = patch_colls[0].get_paths()[0]
    xs = rect.vertices[:, 0]
    ys = rect.vertices[:, 1]
    # x (target axis) keeps forward coords; y mirrors to [len-end, len-start].
    assert xs.min() == 100 and xs.max() == 300
    assert ys.min() == seq_len - 300 and ys.max() == seq_len - 100
    plt.close(fig)


# ------------------------------------------------------------ HTML report


def test_html_report_annotations_payload_and_click_wiring(tmp_path):
    pl = _plotter()
    ann = GffAnnotation.from_text(GFF)
    out = tmp_path / 'report.html'
    pl.to_html(out, annotation=ann)
    html = out.read_text()

    m = re.search(
        r'<script type="application/json" id="rd-data">(.*?)</script>', html, re.S
    )
    payload = json.loads(m.group(1))
    p00 = payload['panels']['rd-panel-0-0']
    annots = p00['annotations']
    # One entry per c1 feature, in draw (file) order.
    assert [a['type'] for a in annots] == [
        'gene',
        'CDS',
        'CDS',
        'gene',
        'repeat_region',
    ]
    first = annots[0]
    assert first['name'] == 'GeneA'
    assert first['strand'] == '+'
    assert first['start'] == 50 and first['end'] == 600

    # The SVG group exists with exactly one child path per feature.
    gidx = html.index('id="rd-annot-0-0"')
    group = html[gidx : html.index('</g>', gidx)]
    assert group.count('<path') == len(annots)
    # Off-diagonal panels carry no annotations.
    assert 'annotations' not in payload['panels']['rd-panel-0-1']
    # Click wiring shipped with the report.
    assert 'showAnnotationDetail' in html
    assert 'rd-annot-' in html


def test_diagonal_squares_on_cross_index_self_alignment():
    """Same contig under two group prefixes counts as a self panel."""
    from rusty_dot.paf_io import CrossIndex

    seq = 'ACGTTGCAAGGCCTTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCG' * 10
    cross = CrossIndex(k=11)
    cross.add_sequence('c1', seq, group='query')
    cross.add_sequence('c1', seq, group='target')
    cross.compute_matches('query', 'target', True)
    ann = GffAnnotation.from_text('c1\tt\tgene\t11\t100\t.\t+\t.\tID=g')
    fig = DotPlotter(cross).plot(
        query_group='query', target_group='target', annotation=ann
    )
    squares = [c for ax in fig.axes for c in ax.collections if c.get_zorder() == 0.5]
    assert squares, 'self panel across group prefixes must draw squares'
    plt.close(fig)
