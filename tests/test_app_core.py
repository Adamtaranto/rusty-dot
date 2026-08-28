"""Tests for the browser app's pure-Python core (app/core)."""

import gzip
from pathlib import Path
import sys

import matplotlib
import pytest

matplotlib.use('Agg')

APP_DIR = Path(__file__).resolve().parent.parent / 'app'
sys.path.insert(0, str(APP_DIR))

from core.align import (  # noqa: E402
    AVAILABLE_METHODS,
    METHOD_LABELS,
    paf_alignment_from_text,
    paf_text_from_alignment,
)
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache  # noqa: E402
from core.fasta import content_digest, parse_fasta_bytes  # noqa: E402
from core.state import ORDER_CHOICES, PlotConfig  # noqa: E402

FASTA = b'>chr1 first contig\nACGTACGTAC\nGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAA\n'
PAF_LINE = 'q1\t100\t0\t50\t+\tt1\t200\t10\t60\t50\t50\t60'


# ---------------------------------------------------------------- fasta


def test_parse_basic_multiline():
    fi = parse_fasta_bytes(FASTA)
    assert fi.names == ['chr1', 'chr2']
    assert fi.records[0][1] == 'ACGTACGTACGTACGTACGT'
    assert fi.records[1][1] == 'TTTTGGGGCCCCAAAA'
    assert fi.total_length == 36


def test_parse_crlf_and_blank_lines():
    data = b'>a\r\nACGT\r\n\r\nACGT\r\n>b\r\nGGGG\r\n'
    fi = parse_fasta_bytes(data)
    assert fi.records == (('a', 'ACGTACGT'), ('b', 'GGGG'))


def test_parse_gzip_roundtrip():
    fi_plain = parse_fasta_bytes(FASTA)
    fi_gz = parse_fasta_bytes(gzip.compress(FASTA))
    assert fi_gz.records == fi_plain.records
    # Digest is of the *raw* bytes, so compressed and plain differ.
    assert fi_gz.digest != fi_plain.digest


def test_parse_semicolon_comment_lines_skipped():
    fi = parse_fasta_bytes(b';comment\n>a\nACGT\n;mid comment\nACGT\n')
    assert fi.records == (('a', 'ACGTACGT'),)


def test_parse_header_name_is_first_token():
    fi = parse_fasta_bytes(b'>contig_1 length=4 cov=12\nACGT\n')
    assert fi.names == ['contig_1']


def test_parse_internal_whitespace_stripped_from_sequence():
    # Spaces and tabs inside sequence lines are removed, not preserved.
    fi = parse_fasta_bytes(b'>a\nAC GT\tAC\n>b\n GG GG \n')
    assert fi.records == (('a', 'ACGTAC'), ('b', 'GGGG'))


def test_parse_no_trailing_newline():
    fi = parse_fasta_bytes(b'>a\nACGT\n>b\nTTTT')
    assert fi.records == (('a', 'ACGT'), ('b', 'TTTT'))


def test_parse_blank_line_before_header():
    fi = parse_fasta_bytes(b'>a\nACGT\n\n\n>b\nTTTT\n')
    assert fi.names == ['a', 'b']


def test_parse_error_line_numbers_reported():
    with pytest.raises(ValueError, match='line 3'):
        parse_fasta_bytes(b'>a\nACGT\n>a\nGGGG\n')
    with pytest.raises(ValueError, match='line 2'):
        parse_fasta_bytes(b';comment\nACGT\n>a\nACGT\n')


@pytest.mark.parametrize(
    ('data', 'match'),
    [
        (b'', 'No FASTA records'),
        (b'ACGT\n>a\nACGT\n', 'before first'),
        (b'>a\nACGT\n>a\nGGGG\n', 'Duplicate contig name'),
        (b'>\nACGT\n', 'Empty contig name'),
        (b'> \nACGT\n', 'Empty contig name'),
        (b'>a\n>b\nACGT\n', 'no sequence'),
        (b'\xff\xfe\x00bad', 'not UTF-8'),
        (b'\x1f\x8btruncated-gzip', 'gzip'),
    ],
)
def test_parse_errors(data, match):
    with pytest.raises(ValueError, match=match):
        parse_fasta_bytes(data)


def test_content_digest_stable_and_distinct():
    assert content_digest(FASTA) == content_digest(FASTA)
    assert content_digest(FASTA) != content_digest(FASTA + b'x')
    assert len(content_digest(FASTA)) == 16


# ---------------------------------------------------------------- align


def test_paf_from_text_skips_comments_and_blanks():
    text = f'# comment\n\n{PAF_LINE}\n\n'
    aln = paf_alignment_from_text(text)
    assert len(aln) == 1
    rec = aln.records[0]
    assert (rec.query_name, rec.target_name) == ('q1', 't1')
    assert rec.strand == '+'


def test_paf_from_text_reports_line_number():
    with pytest.raises(ValueError, match='line 2'):
        paf_alignment_from_text(f'{PAF_LINE}\nnot\tpaf\n')


def test_paf_from_text_empty():
    with pytest.raises(ValueError, match='No PAF records'):
        paf_alignment_from_text('# only a comment\n')


def test_paf_text_roundtrip():
    aln = paf_alignment_from_text(PAF_LINE)
    text = paf_text_from_alignment(aln)
    assert text.endswith('\n')
    again = paf_alignment_from_text(text)
    assert len(again) == len(aln)
    assert again.records[0].query_name == 'q1'


def test_method_registry_consistent():
    assert AVAILABLE_METHODS <= set(METHOD_LABELS)
    assert {'kmer', 'paf'} <= AVAILABLE_METHODS
    # BLAST is deliberately absent (no WASM build exists).
    assert 'blast' not in METHOD_LABELS


# ---------------------------------------------------------------- cache


def _two_assemblies():
    seq = 'ACGTTGCAAGGCTTAACCGGTTAACGGCCAATT' * 8
    q = parse_fasta_bytes(f'>q1\n{seq}\n>q2\n{seq[::-1]}\n'.encode())
    t = parse_fasta_bytes(f'>t1\n{seq}\n'.encode())
    return q, t


def test_kmer_index_builds_and_caches():
    cache = SessionCache()
    q, t = _two_assemblies()
    idx1 = cache.kmer_index(11, q, t)
    idx2 = cache.kmer_index(11, q, t)
    assert idx1 is idx2  # cache hit returns the same object
    idx3 = cache.kmer_index(13, q, t)
    assert idx3 is not idx1  # different k -> rebuild
    # Matches were computed at build time: PAF export works immediately.
    # (get_paf_all is hardcoded to legacy groups 'a'/'b', so the app uses
    # get_paf with explicit group pairs.)
    assert idx1.get_paf(group_pairs=[(QUERY_GROUP, TARGET_GROUP)], merge=True)


def test_kmer_index_group_assignment():
    cache = SessionCache()
    q, t = _two_assemblies()
    idx = cache.kmer_index(11, q, t)
    fig_input_names = idx.sequence_names()
    assert any(n.startswith(f'{QUERY_GROUP}:') for n in fig_input_names)
    assert any(n.startswith(f'{TARGET_GROUP}:') for n in fig_input_names)


def test_paf_cache_keying():
    cache = SessionCache()
    aln = paf_alignment_from_text(PAF_LINE)
    cache.put_paf('paf', {'x': 1}, aln, 'digest-a')
    assert cache.get_paf('paf', {'x': 1}, 'digest-a') is aln
    assert cache.get_paf('paf', {'x': 2}, 'digest-a') is None
    assert cache.get_paf('paf', {'x': 1}, 'digest-b') is None
    # Param order must not matter.
    cache.put_paf('m', {'a': 1, 'b': 2}, aln, 'd')
    assert cache.get_paf('m', {'b': 2, 'a': 1}, 'd') is aln


# ---------------------------------------------------------------- state


def test_plot_config_kwargs_mapping():
    cfg = PlotConfig()
    kwargs = cfg.plot_kwargs()
    assert kwargs['contig_order'] is None  # 'input' maps to None
    assert kwargs['auto_reverse'] is False
    cfg = PlotConfig(contig_order='colinearity', auto_reverse=True)
    kwargs = cfg.plot_kwargs()
    assert kwargs['contig_order'] == 'colinearity'
    assert kwargs['auto_reverse'] is True
    assert {'input', 'length', 'colinearity', 'colinearity_ref'} == set(ORDER_CHOICES)


# ------------------------------------------------------- app-level smoke


def test_cached_index_plots_like_the_app():
    """The exact call path the app's make_figure uses must produce a Figure."""
    from rusty_dot import DotPlotter

    cache = SessionCache()
    q, t = _two_assemblies()
    idx = cache.kmer_index(11, q, t)
    cfg = PlotConfig(contig_order='colinearity', auto_reverse=True)
    fig = DotPlotter(idx).plot(
        query_group=QUERY_GROUP, target_group=TARGET_GROUP, **cfg.plot_kwargs()
    )
    assert fig.axes
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_write_fasta_reordered_like_the_app(tmp_path):
    """FASTA export with explicit order/reverse matches the app download path."""
    cache = SessionCache()
    q, t = _two_assemblies()
    idx = cache.kmer_index(11, q, t)
    order, _t_order = idx.reorder_contigs(
        query_group=QUERY_GROUP, target_group=TARGET_GROUP
    )
    reverse = idx.reversed_contigs(QUERY_GROUP)
    out = tmp_path / 'reordered.fasta'
    idx.write_fasta(out, group=QUERY_GROUP, order=order, reverse=reverse)
    text = out.read_text()
    assert text.startswith('>')
    assert set(order) == {'q1', 'q2'}
    for name in order:
        assert f'>{name}' in text


# ------------------------------------------------------- wheel selection


def test_pick_wasm_wheel_matches_platform():
    from pathlib import Path

    from core.wheels import pick_wasm_wheel

    tag = 'emscripten_3_1_58_wasm32'
    good = Path('wheels/rusty_dot-0.1.0-cp312-cp312-emscripten_3_1_58_wasm32.whl')
    stale = Path('wheels/rusty_dot-0.1.0-cp312-cp312-emscripten_3_1_65_wasm32.whl')
    # A stale wheel from another Emscripten version sorts last lexically —
    # the picker must select by platform tag, not sort order.
    assert pick_wasm_wheel([good, stale], tag) == good
    assert pick_wasm_wheel([stale, good], tag) == good
    # Multiple matches: lexically last (highest version) wins.
    newer = Path('wheels/rusty_dot-0.2.0-cp312-cp312-emscripten_3_1_58_wasm32.whl')
    assert pick_wasm_wheel([good, newer, stale], tag) == newer


def test_pick_wasm_wheel_errors():
    from pathlib import Path

    from core.wheels import pick_wasm_wheel
    import pytest as _pytest

    with _pytest.raises(RuntimeError, match='No rusty-dot wasm wheel'):
        pick_wasm_wheel([], 'emscripten_3_1_58_wasm32')
    stale = Path('wheels/rusty_dot-0.1.0-cp312-cp312-emscripten_3_1_65_wasm32.whl')
    with _pytest.raises(RuntimeError, match='different Pyodide'):
        pick_wasm_wheel([stale], 'emscripten_3_1_58_wasm32')


def test_runtime_platform_tag_shape():
    from core.wheels import runtime_platform_tag

    tag = runtime_platform_tag()
    assert tag and '-' not in tag and '.' not in tag


# --------------------------------------------- PAF/assembly name validation


def test_validate_query_names_all_match():
    from core.validate import validate_query_names

    assert validate_query_names(['q1', 'q2'], ['q1', 'q2'], ['t1']) == []


def test_validate_query_names_none_match():
    from core.validate import validate_query_names

    warnings = validate_query_names(['a', 'b'], ['q1'], ['t1'])
    assert len(warnings) == 1
    assert 'None of the uploaded assembly contig names' in warnings[0]
    assert 'target column' not in warnings[0]


def test_validate_query_names_swapped_inputs_hint():
    from core.validate import validate_query_names

    warnings = validate_query_names(['t1', 't2'], ['q1'], ['t1', 't2'])
    assert len(warnings) == 1
    assert 'DO match the PAF target column' in warnings[0]


def test_validate_query_names_partial_overlap_both_directions():
    from core.validate import validate_query_names

    warnings = validate_query_names(['q1', 'extra'], ['q1', 'ghost'], ['t1'])
    assert len(warnings) == 2
    assert any('no alignments in the PAF' in w and 'extra' in w for w in warnings)
    assert any('not in the uploaded assembly' in w and 'ghost' in w for w in warnings)


def test_validate_query_names_ambiguous_duplicates():
    from core.validate import validate_query_names

    warnings = validate_query_names(['q1'], ['q1', 'shared'], ['t1', 'shared'])
    assert any('BOTH the PAF query and target columns' in w for w in warnings)
    assert any('shared' in w for w in warnings)


def test_validate_query_names_preview_truncates():
    from core.validate import validate_query_names

    missing = [f'c{i}' for i in range(10)]
    warnings = validate_query_names(['q1', *missing], ['q1'], ['t1'])
    assert len(warnings) == 1
    assert '10 assembly contig(s)' in warnings[0]
    assert '…' in warnings[0]


def test_kmer_index_min_block_len_in_cache_key():
    cache = SessionCache()
    q, t = _two_assemblies()
    unfiltered = cache.kmer_index(11, q, t)
    filtered = cache.kmer_index(11, q, t, min_block_len=1000)
    assert filtered is not unfiltered  # distinct cache entries
    # An absurd threshold removes every record.
    assert not filtered.get_paf(group_pairs=[(QUERY_GROUP, TARGET_GROUP)], merge=True)


def test_kmer_cache_evicts_least_recently_used():
    cache = SessionCache()
    q, t = _two_assemblies()
    first = cache.kmer_index(11, q, t)
    for k in (12, 13, 14):  # fill past KMER_MAX=3, evicting k=11
        cache.kmer_index(k, q, t)
    assert len(cache._kmer) == SessionCache.KMER_MAX
    rebuilt = cache.kmer_index(11, q, t)
    assert rebuilt is not first  # k=11 was evicted and rebuilt


def test_paf_cache_evicts_and_refreshes_on_hit():
    from rusty_dot.paf_io import PafAlignment, PafRecord

    cache = SessionCache()
    line = 'q\t100\t0\t50\t+\tt\t100\t0\t50\t50\t50\t255'
    aln = PafAlignment.from_records([PafRecord.from_line(line)])
    for i in range(SessionCache.PAF_MAX):
        cache.put_paf('minimap2', {'i': i}, aln, 'dq', 'dt')
    # Touch the oldest entry, then overflow: the second-oldest is evicted.
    assert cache.get_paf('minimap2', {'i': 0}, 'dq', 'dt') is not None
    cache.put_paf('minimap2', {'i': 99}, aln, 'dq', 'dt')
    assert cache.get_paf('minimap2', {'i': 0}, 'dq', 'dt') is not None
    assert cache.get_paf('minimap2', {'i': 1}, 'dq', 'dt') is None
    assert len(cache._paf) == SessionCache.PAF_MAX


# ------------------------------------------------- annotation controls


def test_slugify_and_type_slug_map():
    from core.annotation_state import slugify, type_slug_map

    assert slugify('five_prime_UTR') == 'five_prime_utr'
    assert slugify('repeat region!') == 'repeat_region'
    assert slugify('***') == 'type'
    # Colliding slugs get deterministic numeric suffixes.
    mapping = type_slug_map(['Gene', 'gene', 'ge ne'])
    assert mapping['Gene'] == 'gene'
    assert mapping['gene'] == 'gene_2'
    assert mapping['ge ne'] == 'ge_ne'
    assert len(set(mapping.values())) == 3


# --------------------------------------------------------------- colours


def test_normalise_type_folds_case_and_aliases():
    from core.annotation_colors import normalise_type

    assert normalise_type('CDS') == normalise_type('cds') == 'cds'
    assert normalise_type('  mRNA ') == 'mrna'
    assert normalise_type('transcript') == 'mrna'
    for alias in ('repeat_region', 'TE', 'LTR', 'tandem_repeat'):
        assert normalise_type(alias) == 'repeat'
    # Unknown types fold case only, they are not invented into aliases.
    assert normalise_type('Effector') == 'effector'


def test_reserved_types_get_their_conventional_colors():
    from core.annotation_colors import RESERVED, assign_shared_colors

    shared = assign_shared_colors(
        {'query': ['gene', 'CDS', 'mRNA', 'exon', 'repeat_region']}
    )
    assert shared['gene'] == RESERVED['gene']
    assert shared['cds'] == RESERVED['cds']
    assert shared['mrna'] == RESERVED['mrna']
    assert shared['exon'] == RESERVED['exon']
    assert shared['repeat'] == RESERVED['repeat']


def test_same_type_gets_one_colour_across_roles():
    """The whole point: query and target must agree, whatever the spelling.

    GffAnnotation._assign_colors indexes into each upload's own sorted type
    list, so 'gene' would otherwise land on a different palette entry per
    file -- and 'CDS' vs 'cds' on different entries entirely.
    """
    from core.annotation_colors import assign_shared_colors, color_map_for

    q_types = ['CDS', 'Effector', 'TE', 'gene']
    t_types = ['cds', 'Starship', 'intron', 'gene', 'repeat_region']
    shared = assign_shared_colors({'query': q_types, 'target': t_types})

    q_map = color_map_for(q_types, shared)
    t_map = color_map_for(t_types, shared)
    assert q_map['gene'] == t_map['gene']
    assert q_map['CDS'] == t_map['cds']  # case-insensitive
    assert q_map['TE'] == t_map['repeat_region']  # aliased


def test_shared_colors_are_order_independent_and_deterministic():
    from core.annotation_colors import assign_shared_colors

    a = assign_shared_colors({'query': ['b', 'a'], 'target': ['c']})
    b = assign_shared_colors({'target': ['c', 'a'], 'query': ['b']})
    assert a == b == assign_shared_colors({'query': ['a', 'b', 'c']})


def test_many_types_never_repeat_a_colour():
    from core.annotation_colors import PALETTE, assign_shared_colors

    types = [f'type_{i:02d}' for i in range(len(PALETTE) + 16)]
    shared = assign_shared_colors({'query': types})
    assert len(shared) == len(types)
    assert len(set(shared.values())) == len(types)


def test_color_map_for_omits_unknown_types():
    from core.annotation_colors import color_map_for

    assert color_map_for(['gene', 'ghost'], {'gene': '#111111'}) == {'gene': '#111111'}


def test_display_name_prefers_the_common_spelling():
    from core.annotation_colors import display_name

    assert display_name(['CDS', 'CDS', 'cds']) == 'CDS'
    assert display_name(['cds', 'CDS']) == 'CDS'  # tie -> alphabetical


def test_apply_annotation_config_filters_and_recolours():
    from core.annotation_state import apply_annotation_config

    from rusty_dot.annotation import GffAnnotation

    ann = GffAnnotation.from_text(
        'c1\tt\tgene\t1\t100\t.\t+\t.\tID=g\nc1\tt\tCDS\t1\t50\t.\t+\t0\tID=c\n'
    )
    out = apply_annotation_config(ann, {'CDS': False}, {'gene': '#112233', 'CDS': ''})
    assert out.feature_types() == ['gene']
    assert out.get_color('gene') == '#112233'
    # The original is untouched.
    assert ann.feature_types() == ['CDS', 'gene']
    # Everything off -> nothing to draw.
    assert apply_annotation_config(ann, {'CDS': False, 'gene': False}, {}) is None


def test_make_figure_call_path_accepts_annotations():
    """plot() accepts the exact annotation kwargs the app forwards."""
    from core.state import PlotConfig

    from rusty_dot import DotPlotter
    from rusty_dot.annotation import GffAnnotation

    cache = SessionCache()
    q, t = _two_assemblies()
    idx = cache.kmer_index(11, q, t)
    ann = GffAnnotation.from_text('q1\tt\tgene\t1\t100\t.\t+\t.\tID=g')
    cfg = PlotConfig()
    kwargs = cfg.plot_kwargs()
    kwargs.update(
        contig_order=None,
        annotation=ann,
        annotation_query=ann,
        annotation_target=None,
        annotation_tracks=True,
    )
    fig = DotPlotter(idx).plot(
        query_group=QUERY_GROUP, target_group=TARGET_GROUP, **kwargs
    )
    assert fig.axes


def test_has_self_pair_detects_shared_contigs():
    from core.panels import has_self_pair

    assert has_self_pair(['c1', 'c2'], ['c2', 'c3'])
    assert not has_self_pair(['c1', 'c2'], ['c3'])
    assert not has_self_pair([], [])


def test_merge_annotations_concatenates_and_keeps_colours():
    from core.annotation_state import merge_annotations

    from rusty_dot.annotation import GffAnnotation

    a = GffAnnotation.from_text('c1\tt\tgene\t1\t10\t.\t+\t.\tID=a')
    b = GffAnnotation.from_text('c1\tt\texon\t20\t30\t.\t+\t.\tID=b')
    a.set_colors({'gene': '#111111'})
    b.set_colors({'exon': '#222222'})

    merged = merge_annotations([a, b])
    assert len(merged) == 2
    assert merged.get_color('gene') == '#111111'
    assert merged.get_color('exon') == '#222222'
    # Degenerate inputs.
    assert merge_annotations([None, None]) is None
    assert merge_annotations([a, None]) is a  # single input passes through


def test_merge_annotations_needs_no_rusty_dot_when_there_is_nothing_to_merge():
    """The empty case must not import rusty_dot.

    Under Shinylive the app boots while the wasm wheel is still installing,
    and the annotation controls are primed at startup with no uploads yet.
    An import there raised ModuleNotFoundError, which a ``reactive.calc``
    caches and re-raises until its dependencies change -- so a transient
    startup failure disabled annotations for the whole session.
    """
    import builtins

    from core.annotation_state import merge_annotations

    real_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.startswith('rusty_dot'):
            raise ModuleNotFoundError("No module named 'rusty_dot'")
        return real_import(name, *args, **kwargs)

    builtins.__import__ = blocked
    try:
        assert merge_annotations([]) is None
        assert merge_annotations([None, None]) is None
    finally:
        builtins.__import__ = real_import


# ------------------------------------------------- per-feature overrides

FEATURE_GFF = (
    'c1\tsrcA\tgene\t1\t100\t.\t+\t.\tID=g1;Name=Alpha\n'
    'c1\tsrcA\tCDS\t10\t50\t.\t+\t0\tID=c1;Parent=g1;product=thing\n'
    'c2\tsrcA\tgene\t1\t80\t.\t-\t.\tID=g2\n'
)


def _feature_ann():
    from rusty_dot.annotation import GffAnnotation

    return GffAnnotation.from_text(FEATURE_GFF)


def test_build_feature_rows_filters_by_sequence_and_uses_1_based_coords():
    from core.annotation_state import build_feature_rows

    rows = build_feature_rows(_feature_ann(), 'c1', 'query')
    assert [r['type'] for r in rows] == ['gene', 'CDS']
    gene = rows[0]
    # Displayed 1-based inclusive, matching the source file.
    assert (gene['start'], gene['end'], gene['length']) == (1, 100, 100)
    assert gene['uid'] == 'query:0'
    assert gene['name'] == 'Alpha'
    assert gene['source'] == 'srcA'
    # uids index the whole record list, not the filtered subset.
    assert rows[1]['uid'] == 'query:1'
    assert build_feature_rows(None, 'c1', 'query') == []


def test_build_feature_rows_surfaces_extra_attributes():
    from core.annotation_state import build_feature_rows

    rows = build_feature_rows(_feature_ann(), 'c1', 'query')
    assert rows[1]['attributes'] == {'product': 'thing'}
    assert 'ID' not in rows[1]['attributes']  # promoted to its own column


def test_build_feature_rows_carries_source_file():
    from core.annotation_state import build_feature_rows

    ann = _feature_ann()
    for rec in ann.records:
        rec.source_file = 'genes.gff'
    assert build_feature_rows(ann, 'c1', 'query')[0]['source_file'] == 'genes.gff'


def test_apply_feature_overrides_hides_exactly_the_named_features():
    from core.annotation_state import apply_feature_overrides

    ann = _feature_ann()
    out = apply_feature_overrides(ann, frozenset({'query:1'}), {}, 'query')
    assert [r.feature_id for r in out.records] == ['g1', 'g2']


def test_apply_feature_overrides_sets_and_resets_colour():
    from core.annotation_state import apply_feature_overrides

    ann = _feature_ann()
    out = apply_feature_overrides(ann, frozenset(), {'query:0': '#0f0f0f'}, 'query')
    assert out.records[0].color == '#0f0f0f'
    assert out.records[1].color is None
    # Dropping the override must fall back to the type colour.
    out = apply_feature_overrides(out, frozenset(), {}, 'query')
    assert out.records[0].color is None


def test_apply_feature_overrides_hiding_everything_returns_none():
    from core.annotation_state import apply_feature_overrides

    ann = _feature_ann()
    uids = frozenset(f'query:{i}' for i in range(len(ann.records)))
    assert apply_feature_overrides(ann, uids, {}, 'query') is None
    assert apply_feature_overrides(None, frozenset(), {}, 'query') is None


def test_type_toggle_beats_feature_toggle():
    """apply_annotation_config runs last, so a type switched off wins."""
    from core.annotation_state import apply_annotation_config, apply_feature_overrides

    ann = _feature_ann()
    kept = apply_feature_overrides(ann, frozenset(), {}, 'query')
    out = apply_annotation_config(kept, {'gene': True, 'CDS': False}, {})
    assert {r.feature_type for r in out.records} == {'gene'}


def test_feature_colour_survives_the_type_filter():
    """keep_feature_types reuses record objects, so overrides must persist."""
    from core.annotation_state import apply_annotation_config, apply_feature_overrides

    ann = _feature_ann()
    kept = apply_feature_overrides(ann, frozenset(), {'query:0': '#abcabc'}, 'query')
    out = apply_annotation_config(kept, {'gene': True, 'CDS': True}, {})
    gene = next(r for r in out.records if r.feature_id == 'g1')
    assert gene.color == '#abcabc'


def test_feature_table_js_contract():
    """Pin the app <-> feature-table.js bridge at the asset level.

    There is no headless-browser harness, so a rename on either side would
    otherwise silently leave the annotations table inert.
    """
    js = (APP_DIR / 'www' / 'feature-table.js').read_text()
    app_py = (APP_DIR / 'app.py').read_text()

    assert "setInputValue('feature_table_change'" in js
    for kind in ("'vis'", "'color'", "'bulk'"):
        assert kind in js
    assert 'rd-feature-table' in js and 'rd-feature-table' in app_py
    assert 'input.feature_table_change' in app_py
    # The table must stay unbound: one Shiny input per feature would mean
    # ~1200 bindings on a real GFF and a visible Pyodide freeze.
    assert 'input_checkbox' not in app_py.split('def annotation_table')[1]


def test_report_frame_is_not_suspended_when_hidden():
    """Switching to the Annotations tab must not re-render the report.

    Shiny suspends hidden outputs by default, so returning to the Plot tab
    would trigger a fresh seconds-long matplotlib pass.
    """
    app_py = (APP_DIR / 'app.py').read_text()
    head = app_py.split('def report_frame')[0]
    assert head.rstrip().endswith('@render.ui')
    assert 'suspend_when_hidden=False' in head.rsplit('@output', 1)[1]


def test_new_result_clears_the_paf_pair_index():
    """A re-run from the overview must drop the previous run's records.

    ``paf_pair_index`` is otherwise cleared only on a ``focus`` change, and
    ``focus.set(None)`` is a no-op when the user is already in the overview
    -- sequence previews then resolved CIGARs (or their absence) from the
    stale records of the previous aligner run.
    """
    app_py = (APP_DIR / 'app.py').read_text()
    body = app_py.split('def _reset_focus_on_new_result')[1]
    body = body.split('@reactive.effect')[0]
    assert 'paf_pair_index.clear()' in body


# ------------------------------------------------------ annotation sources


def _ann(text='c1\tt\tgene\t1\t10\t.\t+\t.\tID=a'):
    from rusty_dot.annotation import GffAnnotation

    return GffAnnotation.from_text(text)


def test_replace_source_reports_no_change_for_an_identical_upload():
    """Re-running an alignment re-parses the same GenBank file every time.

    The app's reactive values invalidate on identity, so re-setting an
    equal-but-new tuple would rebuild the feature-type controls and discard
    the user's toggles and colours on every run.
    """
    from core.annotation_state import replace_source

    key = ('digest1', 'a.gbk')
    entries = replace_source((), 'genbank', 'a.gbk', _ann(), key)
    assert entries is not None and len(entries) == 1
    assert replace_source(entries, 'genbank', 'a.gbk', _ann(), key) is None


def test_replace_source_detects_new_content_under_the_same_name():
    from core.annotation_state import replace_source

    entries = replace_source((), 'genbank', 'a.gbk', _ann(), ('d1', 'a.gbk'))
    updated = replace_source(entries, 'genbank', 'a.gbk', _ann(), ('d2', 'a.gbk'))
    assert updated is not None
    assert updated[0]['key'] == ('d2', 'a.gbk')


def test_replace_source_detects_the_same_content_under_a_new_name():
    """source_file is displayed in the drill-down, so the name must refresh."""
    from core.annotation_state import replace_source

    entries = replace_source((), 'gff', 'a.gff', _ann(), ('d1', 'a.gff'))
    updated = replace_source(entries, 'gff', 'b.gff', _ann(), ('d1', 'b.gff'))
    assert updated is not None
    assert updated[0]['filename'] == 'b.gff'


def test_replace_source_clearing_an_empty_slot_is_a_no_op():
    """Otherwise an empty GFF upload wipes every per-feature override."""
    from core.annotation_state import replace_source

    assert replace_source((), 'gff', '', None, None) is None
    entries = replace_source((), 'gff', 'a.gff', _ann(), ('d1', 'a.gff'))
    cleared = replace_source(entries, 'gff', '', None, None)
    assert cleared == ()


def test_replace_source_leaves_the_other_kind_alone():
    """A role can hold GenBank-derived features and an uploaded GFF."""
    from core.annotation_state import replace_source

    entries = replace_source((), 'genbank', 'a.gbk', _ann(), ('d1', 'a.gbk'))
    both = replace_source(entries, 'gff', 'b.gff', _ann(), ('d2', 'b.gff'))
    assert {e['kind'] for e in both} == {'genbank', 'gff'}
    # Replacing the gff must not disturb the genbank entry...
    again = replace_source(both, 'gff', 'c.gff', _ann(), ('d3', 'c.gff'))
    assert {e['kind'] for e in again} == {'genbank', 'gff'}
    gb = next(e for e in again if e['kind'] == 'genbank')
    assert gb['key'] == ('d1', 'a.gbk')
    # ...and an unchanged genbank still short-circuits alongside a gff.
    assert replace_source(again, 'genbank', 'a.gbk', _ann(), ('d1', 'a.gbk')) is None


# --------------------------------------------------------------- nav tips


def test_nav_tips_hides_click_to_focus_on_a_single_panel():
    """Click-to-focus is disabled in the report when there is one panel."""
    from core.panels import nav_tips

    actions = lambda f, m: [a for a, _ in nav_tips(f, m)]  # noqa: E731

    assert 'click panel' in actions(False, True)
    assert 'click panel' not in actions(False, False)
    assert 'click panel' not in actions(True, True)  # already focused
    # Double-click drill-down works even on a single-panel overview.
    assert 'double-click panel' in actions(False, False)
    assert 'double-click panel' not in actions(True, False)
    # The pan/zoom tips are unconditional.
    for combo in ((True, True), (True, False), (False, True), (False, False)):
        assert {'scroll', 'drag', 'Esc'} <= set(actions(*combo))


# ------------------------------------------------ asset-level UI contracts


def test_report_js_matches_panel_ids_strictly():
    """Prefix-only matching let the axes background pose as a panel."""
    import rusty_dot._html as html_pkg

    js = (Path(html_pkg.__file__).parent / 'report.js').read_text()
    assert 'PANEL_ID_RE' in js
    assert 'function closestPanel(' in js
    # panelGroups must be filtered, not taken raw from the prefix selector.
    head = js.split('var panelGroups')[1].split('var selectedPanel')[0]
    assert 'PANEL_ID_RE.test' in head
    # The colliding id is gone everywhere.
    css = (Path(html_pkg.__file__).parent / 'report.css').read_text()
    dotplot = (Path(html_pkg.__file__).parent.parent / 'dotplot.py').read_text()
    for src in (js, css, dotplot):
        assert 'rd-panel-0-0-bg' not in src
        assert "f'{gid}-bg'" not in src


def test_panel_dblclick_bridge_walks_ancestors():
    """closest() with a prefix selector stopped at the background group."""
    app_py = (APP_DIR / 'app.py').read_text()
    bridge = app_py.split('_PANEL_DBLCLICK_JS = ')[1].split('"""', 2)[1]
    assert '.closest(' not in bridge  # the call, not the word in a comment
    assert 'parentNode' in bridge
    # The exact-shape regex is what makes the walk terminate correctly, and
    # tests/test_app_state.py pins that it survives Python escaping.
    assert r'/^rd-panel-(\\d+)-(\\d+)$/' in bridge


def test_annotation_toggles_are_static_so_a_rerun_cannot_reset_them():
    """gff_controls must depend on the type index alone.

    Anything that re-renders it recreates every gtyp_/gcol_ input at its
    default, silently discarding the user's annotation choices.
    """
    app_py = (APP_DIR / 'app.py').read_text()
    before, body = app_py.split('def gff_controls')
    body = body.split('@reactive.calc')[0]
    # Comments in this file legitimately name the things being asserted
    # against, so test the code only.
    code = '\n'.join(
        line for line in body.splitlines() if not line.lstrip().startswith('#')
    )

    assert 'self_panels()' not in code
    for toggle in ("'gff_diagonal'", "'gff_tracks'"):
        assert toggle not in code, f'{toggle} must live outside gff_controls'
        assert toggle in before
    # Visibility is driven by a hidden output, which must be exempt from
    # Shiny's suspend-when-hidden behaviour to update at all.
    assert 'def gff_mode' in app_py
    gff_mode_head = app_py.split('def gff_mode')[0].rsplit('@output', 1)[1]
    assert 'suspend_when_hidden=False' in gff_mode_head
    # Conditions test positively: output.gff_mode is undefined before the
    # first report, and `undefined !== ''` would flash the controls.
    assert "output.gff_mode === 'self'" in before
    assert 'output.gff_mode !==' not in before


def test_panel_click_is_delegated_so_any_click_can_reset():
    """Focus must be escapable by clicking anywhere, not only the same panel.

    Per-panel listeners each called stopPropagation, so no click outside the
    focused panel could ever reach a reset.
    """
    import rusty_dot._html as html_pkg

    js = (Path(html_pkg.__file__).parent / 'report.js').read_text()
    block = js.split('if (panelGroups.length > 1)')[1].split('// ---')[0]
    assert "svg.addEventListener('click'" in block
    assert 'panelGroups.forEach' not in block, 'still one listener per panel'
    assert 'selectedPanel !== null' in block, 'no click-anywhere reset branch'
    assert 'closestPanel(evt.target)' in block


# ----------------------------------------------- annotation name validation


def test_annotation_names_all_present_is_silent():
    from core.validate import validate_annotation_names

    assert validate_annotation_names(['c1', 'c2'], ['c1'], 'query') == []
    # Nothing to say when either side is empty.
    assert validate_annotation_names([], ['c1'], 'query') == []
    assert validate_annotation_names(['c1'], [], 'query') == []


def test_annotation_names_partial_mismatch_lists_the_missing():
    from core.validate import validate_annotation_names

    (msg,) = validate_annotation_names(['c1', 'c2'], ['c1', 'c9'], 'query')
    assert 'c9' in msg and 'query' in msg
    assert 'c1' not in msg  # only the offenders are named


def test_annotation_names_total_mismatch_suggests_the_likely_cause():
    """Wrong file or wrong role is the usual reason nothing lines up."""
    from core.validate import validate_annotation_names

    (msg,) = validate_annotation_names(['c1', 'c2'], ['zz', 'yy'], 'target')
    assert 'None of the target GFF' in msg
    assert 'right role' in msg


def test_annotation_names_preview_is_capped():
    from core.validate import validate_annotation_names

    missing = [f'x{i}' for i in range(12)]
    (msg,) = validate_annotation_names(['c1'], ['c1'] + missing, 'query')
    assert '12 sequence name(s)' in msg
    assert msg.count(',') <= 5  # 5 names shown then an ellipsis
    assert '…' in msg


def test_clear_annotations_resets_the_widget_not_just_the_server():
    """Shiny's file-input binding has a no-op setValue.

    Without the custom message the filename would linger after clearing,
    and re-picking the same file would fire no change event, so it could
    never be reloaded.
    """
    app_py = (APP_DIR / 'app.py').read_text()
    bridge = (APP_DIR / 'www' / 'bridge.js').read_text()

    assert "'clear_gff'" in app_py
    assert "send_custom_message(\n            'rd_clear_file_inputs'" in app_py
    assert "addCustomMessageHandler(\n        'rd_clear_file_inputs'" in bridge
    # The three pieces of DOM the binding writes and we must undo.
    for reset in ("el.value = ''", "text.value = ''", 'progress-bar'):
        assert reset in bridge


def test_sidebar_scroll_is_captured_before_the_picker_opens():
    """The jump happens during the picker interaction, not the re-render.

    By the time `change` fires the container is already at the top, so
    capturing there would restore 0.  The position has to be taken when
    the input is clicked.
    """
    js = (APP_DIR / 'www' / 'sidebar-scroll.js').read_text()
    click_block = js.split("'click'")[1].split('document.addEventListener')[0]
    assert 'mark = {' in click_block, 'position not captured on click'
    change_block = js.split("'change'")[1]
    assert 'mark && mark.el === scroller' in change_block, (
        'change handler must prefer the pre-picker position'
    )
    # A real gesture has to be able to cancel the hold.
    assert "'wheel'" in js and 'release' in js


def test_self_align_drops_target_annotations_but_spares_paf():
    """Self-alignment puts the query assembly on both axes.

    A target annotation set then has nothing of its own to annotate, and
    its upload is hidden, so leaving it loaded would keep it in the merged
    set with no visible control to remove it.  PAF input is exempt: it has
    both roles and no self-alignment, and `self_align` keeps its last value
    while its panel is hidden.
    """
    app_py = (APP_DIR / 'app.py').read_text()
    body = app_py.split('def _drop_target_annotations_when_self_aligning')[1]
    body = body.split('\n    _MIN_LEN_NOTIF_ID')[0]

    assert "input.input_mode() == 'paf'" in body, 'PAF must be exempt'
    assert 'not input.self_align()' in body
    # Both source kinds go, and the widget is reset too -- Shiny cannot
    # clear a file input, so the filename would otherwise linger.
    assert "_set_ann_source('target', kind, '', None)" in body
    assert "'rd_clear_file_inputs'" in body
    assert "'target_gff'" in body
    # The query role is never touched.
    assert "'query'" not in body

    # The upload is hidden on the matching condition.
    assert '"input.input_mode === \'paf\' || !input.self_align"' in app_py, (
        'target_gff visibility must match the drop condition'
    )


def test_feature_type_toggles_are_held_until_apply():
    """Each toggle otherwise re-renders the whole figure.

    Choices made before the first plot apply themselves when it is
    generated; afterwards they wait for the sidebar's Apply button.
    """
    app_py = (APP_DIR / 'app.py').read_text()

    # annotations() must go through the applied values, not the raw inputs.
    body = app_py.split('def annotations')[1].split('# --- end GFF')[0]
    assert 'gff_type_applied()' in body
    for raw in ('gtyp_', 'gcol_', 'gff_type_choices_live'):
        assert raw not in body, f'annotations() still reads {raw}* directly'

    # The commit runs on the button press and nothing else.
    apply_fn = app_py.split('def _apply_gff_types')[1].split('@reactive.effect')[0]
    assert 'gff_type_applied.set(live)' in apply_fn
    head = app_py.split('def _apply_gff_types')[0].rsplit('@reactive.effect', 1)[1]
    assert '@reactive.event(input.apply_gff_types' in head

    # Before the first plot the live values flow straight through; the
    # mirror backs off as soon as a result exists.
    auto = app_py.split('def _sync_types_before_first_plot')[1]
    auto = auto.split('@reactive.effect')[0]
    assert 'if result() is not None:' in auto and 'return' in auto
    assert 'gff_type_applied.set(live)' in auto


# ------------------------------------------------ pending feature overrides


def test_count_pending_overrides_counts_distinct_features():
    """The button says how many features change, not how many edits."""
    from core.annotation_state import count_pending_overrides

    # One feature both hidden and recoloured still counts once.
    assert (
        count_pending_overrides(
            frozenset({'query:0'}), {'query:0': '#fff'}, frozenset(), {}
        )
        == 1
    )
    assert (
        count_pending_overrides(
            frozenset({'query:0'}), {'query:1': '#fff'}, frozenset(), {}
        )
        == 2
    )


def test_count_pending_overrides_is_zero_when_settled():
    from core.annotation_state import count_pending_overrides

    hidden, colors = frozenset({'query:3'}), {'query:4': '#abc'}
    assert count_pending_overrides(hidden, colors, hidden, dict(colors)) == 0
    assert count_pending_overrides(frozenset(), {}, frozenset(), {}) == 0


def test_count_pending_overrides_sees_reverted_and_changed_edits():
    from core.annotation_state import count_pending_overrides

    # Un-hiding something that is currently hidden is a pending change.
    assert count_pending_overrides(frozenset(), {}, frozenset({'q:1'}), {}) == 1
    # So is dropping a colour override, or changing one.
    assert count_pending_overrides(frozenset(), {}, frozenset(), {'q:1': '#a'}) == 1
    assert (
        count_pending_overrides(frozenset(), {'q:1': '#b'}, frozenset(), {'q:1': '#a'})
        == 1
    )


def test_feature_edits_are_held_until_apply():
    """Per-feature edits must not redraw the figure on every click.

    ~600 rows per sequence, and each toggle previously re-rendered both the
    figure and the whole table.
    """
    app_py = (APP_DIR / 'app.py').read_text()

    handler = app_py.split('def _on_feature_table_change')[1]
    handler = handler.split('@reactive.effect')[0]
    assert 'feature_hidden_pending.set' in handler
    assert 'feature_colors_pending.set' in handler
    assert 'feature_hidden.set' not in handler, 'handler must not commit'

    apply_fn = app_py.split('def _apply_feature_overrides_now')[1]
    apply_fn = apply_fn.split('@reactive.effect')[0]
    assert 'feature_hidden.set(feature_hidden_pending())' in apply_fn
    # Equal-but-new objects would still invalidate and cost a render.
    assert 'return' in apply_fn

    # annotations() draws from the committed values only.
    ann = app_py.split('def annotations')[1].split('# --- end GFF')[0]
    assert 'feature_hidden()' in ann and '_pending' not in ann


def test_annotation_table_does_not_rerender_on_every_edit():
    """The row payload reads pending values under isolate() for this reason."""
    app_py = (APP_DIR / 'app.py').read_text()
    body = app_py.split('def _feature_table_payload')[1].split('\n    @')[0]
    assert 'with reactive.isolate():' in body
    assert 'feature_hidden_pending()' in body

    # The rows themselves are built client-side; the shell must not loop
    # over them (that render cost is what stalled the tab switch).
    shell = app_py.split('def annotation_table')[1].split('\n    @')[0]
    assert 'for r in rows' not in shell
    assert '<tbody></tbody>' in shell

    # With no re-render, the reset button has to repaint the swatches, so
    # each one carries the type colour it falls back to.
    js = (APP_DIR / 'www' / 'feature-table.js').read_text()
    assert 'data-type-color=' in js
    reset = js.split("action === 'reset'")[1]
    assert 'dataset.typeColor' in reset


# ------------------------------------------------ annotations table sorting


def test_feature_table_headers_declare_how_they_sort():
    """The client sorts, so every column has to say how it compares."""
    # Read the constant out of the source: importing app.py needs shiny,
    # which the test environment deliberately does not carry.
    import ast

    tree = ast.parse((APP_DIR / 'app.py').read_text())
    spec = next(
        node.value
        for node in tree.body
        if isinstance(node, ast.AnnAssign)
        and getattr(node.target, 'id', None) == '_FEATURE_COLUMNS'
    )
    kinds = dict(ast.literal_eval(spec))
    assert kinds['Start'] == kinds['End'] == kinds['Length'] == 'num'
    assert kinds['Show'] == 'check'
    assert kinds['Colour'] == 'color'
    assert kinds['Type'] == kinds['Name'] == kinds['Contig'] == 'text'

    app_py = (APP_DIR / 'app.py').read_text()
    body = app_py.split('def annotation_table')[1].split('\n    @')[0]
    assert 'data-sort=' in body and 'aria-sort=' in body
    # Rows remember their file position so the sort can be undone; they
    # are built client-side, so the marker lives in the row template.
    js = (APP_DIR / 'www' / 'feature-table.js').read_text()
    assert 'data-idx=' in js


def test_feature_table_sorts_on_the_client():
    """Sorting must not go through the server.

    Re-rendering hundreds of rows is the cost this table is built to
    avoid, and it would discard edits that have not been applied yet.
    """
    js = (APP_DIR / 'www' / 'feature-table.js').read_text()
    sort = js.split('// --- column sorting')[1].split('// --- client-side filter')[0]

    assert 'setInputValue' not in sort, 'sorting must stay client-side'
    # Rows are moved, not rebuilt, so their inputs keep their state.
    assert 'appendChild' in sort and 'createDocumentFragment' in sort
    # Numeric columns are formatted with thousands separators.
    assert "replace(/,/g, '')" in sort
    # asc -> desc -> file order, and ties keep file order so repeated
    # sorts do not depend on the previous one.
    assert 'data-idx' in sort
    assert 'aria-sort' in sort


def test_drilldown_plot_pane_is_a_flex_column_like_the_overview():
    """The aligner log must never be laid out over the report frame.

    In the drill-down `plot_area` wraps the frame in `ui.navset_tab`, whose
    `.tab-content` / `.tab-pane` are plain block boxes -- and Shiny's output
    wrappers are `display: contents`, so the frame is a flex item of
    whatever box encloses it.  Under block layout `.rd-plot-area` measured
    ~90px shorter than the frame it contained, and `aligner_log_ui` (the
    next sibling in the page body, with no z-index on either side) was laid
    out inside that overflow and painted on top of the frame -- hiding the
    report's fixed bottom detail bar, which is where a clicked annotation's
    metadata appears.

    These three declarations are what keeps the pane measuring its own
    contents; dropping any of them brings the overlap straight back.
    """
    css = (APP_DIR / 'www' / 'app.css').read_text()

    def rule(selector: str) -> str:
        head, _, rest = css.partition(selector + ' {')
        assert rest, f'no rule for {selector}'
        assert head.rstrip().endswith(('}', '*/')), f'{selector} matched mid-rule'
        return rest.split('}')[0]

    # The frame itself: an inline replaced element is sized through a line
    # box, which leaves its wrapper measuring 0 tall.
    assert 'display: block;' in rule('.rd-report-frame')

    for selector in (
        '.rd-plot-area .tab-content',
        '.rd-plot-area .tab-content > .tab-pane.active',
    ):
        body = rule(selector)
        assert 'display: flex;' in body, f'{selector} is not a flex container'
        assert 'flex-direction: column;' in body, f'{selector} is not a column'
        assert 'min-height: 0;' in body, f'{selector} cannot shrink'

    # The annotations tab is the same story: its scroller has to shrink with
    # the pane rather than hold a fixed 62vh and overflow it.
    for selector in ('.rd-ft-panel', '.rd-ft-scroll'):
        body = rule(selector)
        assert 'flex: 1 1 auto;' in body, f'{selector} does not grow/shrink'
        assert 'min-height: 0;' in body, f'{selector} cannot shrink'


# ------------------------------------------------- feature-row dedupe / uids


def test_apply_feature_overrides_stamps_uids_on_records():
    from core.annotation_state import apply_feature_overrides

    ann = _feature_ann()
    out = apply_feature_overrides(ann, frozenset(), {}, 'query')
    assert [r.uid for r in out.records] == ['query:0', 'query:1', 'query:2']
    # keep_feature_types reuses the record objects, so the uid survives the
    # type filter that runs after the override pass.
    filtered = out.keep_feature_types(['gene'])
    assert [r.uid for r in filtered.records] == ['query:0', 'query:2']


def test_dedupe_feature_rows_drops_cross_role_copies():
    from core.annotation_state import build_feature_rows, dedupe_feature_rows

    q = build_feature_rows(_feature_ann(), 'c1', 'query')
    t = build_feature_rows(_feature_ann(), 'c1', 'target')
    # Self panel with the same upload under both roles: every feature
    # would list twice; the first (query-role) copy wins.
    assert dedupe_feature_rows(q + t) == q


def test_dedupe_feature_rows_keeps_distinct_sources():
    from core.annotation_state import build_feature_rows, dedupe_feature_rows

    a = _feature_ann()
    b = _feature_ann()
    for rec in a.records:
        rec.source_file = 'genes.gb'
    for rec in b.records:
        rec.source_file = 'genes.gff'
    rows = build_feature_rows(a, 'c1', 'query') + build_feature_rows(b, 'c1', 'target')
    # Same coordinates from different files is deliberate (GenBank + GFF
    # both loaded); both copies stay.
    assert dedupe_feature_rows(rows) == rows


def test_count_pending_type_changes():
    from core.annotation_state import count_pending_type_changes

    assert count_pending_type_changes({}, {}) == 0
    live = {'gene': (True, ''), 'cds': (False, '#112233')}
    assert count_pending_type_changes(live, dict(live)) == 0
    # A missing key means the default (visible, type colour), so all-default
    # controls against a freshly reset applied side show nothing pending.
    assert count_pending_type_changes({'gene': (True, '')}, {}) == 0
    assert count_pending_type_changes(live, {}) == 1
    assert count_pending_type_changes(live, {'gene': (False, '')}) == 2
    # A colour-only edit counts as a change.
    assert (
        count_pending_type_changes({'gene': (True, '#ffffff')}, {'gene': (True, '')})
        == 1
    )
