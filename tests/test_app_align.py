"""Tests for the browser app's biowasm aligner support (app/core/align.py)."""

from pathlib import Path
import sys

import pytest

APP_DIR = Path(__file__).resolve().parent.parent / 'app'
sys.path.insert(0, str(APP_DIR))

from core.align import (  # noqa: E402
    AVAILABLE_METHODS,
    BIOWASM_TOOLS,
    LASTZ_GENERAL_FIELDS,
    METHOD_LABELS,
    MINIMAP2_PRESETS,
    QUERY_FILENAME,
    TARGET_FILENAME,
    alignment_from_tool_output,
    build_tool_args,
    fasta_text,
    lastz_general_to_records,
    nucmer_delta_to_records,
)
from core.cache import SessionCache  # noqa: E402
from core.fasta import parse_fasta_bytes  # noqa: E402

# ---------------------------------------------------------------- registry


def test_biowasm_methods_available():
    assert BIOWASM_TOOLS == {'minimap2', 'lastz', 'nucmer'}
    assert BIOWASM_TOOLS <= AVAILABLE_METHODS
    assert BIOWASM_TOOLS <= set(METHOD_LABELS)
    # BLAST must not appear: no production WASM build exists.
    assert not any('blast' in m.lower() for m in METHOD_LABELS)


# ---------------------------------------------------------------- fasta_text


def test_fasta_text_roundtrip():
    raw = b'>chr1 description here\nACGT\nACGT\n>chr2\nTTTT\n'
    fi = parse_fasta_bytes(raw)
    text = fasta_text(fi.records)
    assert text == '>chr1\nACGTACGT\n>chr2\nTTTT\n'
    # Re-parsing the reconstructed text yields identical records.
    assert parse_fasta_bytes(text.encode()).records == fi.records


def test_fasta_text_empty_records_raises():
    with pytest.raises(ValueError, match='No FASTA records'):
        fasta_text([])


@pytest.mark.parametrize('bad', [[('', 'ACGT')], [('chr1', '')]])
def test_fasta_text_empty_name_or_seq_raises(bad):
    with pytest.raises(ValueError, match='empty name or sequence'):
        fasta_text(bad)


# ---------------------------------------------------------------- build_tool_args


def test_minimap2_args_defaults():
    args = build_tool_args('minimap2', {'preset': 'asm20', 'k': 0, 'w': 0})
    assert args == ['-x', 'asm20', TARGET_FILENAME, QUERY_FILENAME]


@pytest.mark.parametrize('preset', MINIMAP2_PRESETS)
def test_minimap2_args_presets(preset):
    args = build_tool_args('minimap2', {'preset': preset})
    assert args[:2] == ['-x', preset]


def test_minimap2_args_k_and_w():
    args = build_tool_args('minimap2', {'preset': 'asm5', 'k': 19, 'w': 10})
    assert args == ['-x', 'asm5', '-k', '19', '-w', '10', 'target.fa', 'query.fa']


def test_minimap2_args_bad_preset():
    with pytest.raises(ValueError, match='preset'):
        build_tool_args('minimap2', {'preset': 'map-ont'})


def test_minimap2_args_negative_k():
    with pytest.raises(ValueError, match='-k'):
        build_tool_args('minimap2', {'preset': 'asm5', 'k': -3})


def test_lastz_args_defaults():
    args = build_tool_args('lastz', {'step': 1, 'gapped': True, 'notransition': False})
    assert args == [
        'target.fa[multiple]',
        'query.fa',
        '--format=general:' + ','.join(LASTZ_GENERAL_FIELDS),
    ]
    # Target precedes query and is marked multi-sequence.
    assert args[0].startswith(TARGET_FILENAME)


def test_lastz_args_all_options():
    args = build_tool_args('lastz', {'step': 10, 'gapped': False, 'notransition': True})
    assert '--step=10' in args
    assert '--notransition' in args
    assert '--nogapped' in args


def test_lastz_args_bad_step():
    with pytest.raises(ValueError, match='step'):
        build_tool_args('lastz', {'step': 0})


def test_lastz_default_step_is_assembly_scale():
    # Omitting step falls back to 20 (large-genome seeding), not the tool's
    # step=1 default that is unusably slow on whole assemblies in wasm.
    args = build_tool_args('lastz', {})
    assert '--step=20' in args


def test_nucmer_default_params_are_assembly_scale():
    # Omitting -l/-c falls back to whole-genome settings (100/200) rather
    # than mummer's small-region defaults (20/65).
    args = build_tool_args('nucmer', {})
    assert args[:4] == ['-l', '100', '-c', '200']


def test_nucmer_args_defaults():
    args = build_tool_args('nucmer', {'l': 20, 'c': 65, 'maxmatch': False})
    assert args == ['-l', '20', '-c', '65', TARGET_FILENAME, QUERY_FILENAME]


def test_nucmer_args_maxmatch():
    args = build_tool_args('nucmer', {'l': 15, 'c': 40, 'maxmatch': True})
    assert args[0] == '--maxmatch'
    assert args[1:] == ['-l', '15', '-c', '40', 'target.fa', 'query.fa']


@pytest.mark.parametrize('params', [{'l': 0}, {'c': -1}])
def test_nucmer_args_bad_values(params):
    with pytest.raises(ValueError, match='nucmer'):
        build_tool_args('nucmer', params)


def test_build_tool_args_unknown_tool():
    with pytest.raises(ValueError, match='Unknown biowasm tool'):
        build_tool_args('blast', {})


# ---------------------------------------------------------------- LASTZ parser

LASTZ_HEADER = '#' + '\t'.join(LASTZ_GENERAL_FIELDS)


def _lastz_line(
    name1='t1',
    zstart1=10,
    end1=60,
    length1=200,
    name2='q1',
    strand2='+',
    zstart2=0,
    end2=50,
    length2=100,
    nmatch=48,
):
    return '\t'.join(
        str(v)
        for v in [
            name1,
            zstart1,
            end1,
            length1,
            name2,
            strand2,
            zstart2,
            end2,
            length2,
            nmatch,
        ]
    )


def test_lastz_parser_forward():
    text = LASTZ_HEADER + '\n' + _lastz_line() + '\n'
    (rec,) = lastz_general_to_records(text)
    assert rec.query_name == 'q1'
    assert rec.target_name == 't1'
    assert (rec.query_start, rec.query_end, rec.query_len) == (0, 50, 100)
    assert (rec.target_start, rec.target_end, rec.target_len) == (10, 60, 200)
    assert rec.strand == '+'
    assert rec.residue_matches == 48
    assert rec.alignment_block_len == 50
    assert rec.mapping_quality == 255


def test_lastz_parser_reverse_strand_forward_coords():
    # zstart2+/end2+ are already forward-strand coords; the parser must
    # keep them as-is (PAF requires forward query coords for '-' records).
    text = _lastz_line(strand2='-', zstart2=20, end2=70)
    (rec,) = lastz_general_to_records(text)
    assert rec.strand == '-'
    assert rec.query_start == 20
    assert rec.query_end == 70
    assert rec.query_start < rec.query_end


def test_lastz_parser_multi_contig():
    text = '\n'.join(
        [
            LASTZ_HEADER,
            _lastz_line(name1='t1', name2='q1'),
            _lastz_line(name1='t2', name2='q2', strand2='-'),
            _lastz_line(name1='t1', name2='q2'),
        ]
    )
    recs = lastz_general_to_records(text)
    assert [(r.query_name, r.target_name) for r in recs] == [
        ('q1', 't1'),
        ('q2', 't2'),
        ('q2', 't1'),
    ]


def test_lastz_parser_blank_lines_and_header_skipped():
    text = '\n\n' + LASTZ_HEADER + '\n' + _lastz_line() + '\n\n'
    assert len(lastz_general_to_records(text)) == 1


def test_lastz_parser_wrong_column_count():
    with pytest.raises(ValueError, match='columns'):
        lastz_general_to_records('t1\t10\t60\n')


def test_lastz_parser_non_numeric():
    with pytest.raises(ValueError, match='line 1'):
        lastz_general_to_records(_lastz_line(zstart1='abc'))


def test_lastz_parser_bad_strand():
    with pytest.raises(ValueError, match='strand'):
        lastz_general_to_records(_lastz_line(strand2='.'))


def test_lastz_parser_empty():
    with pytest.raises(ValueError, match='no alignments'):
        lastz_general_to_records(LASTZ_HEADER + '\n')


# ---------------------------------------------------------------- delta parser

DELTA_FWD = """\
/data/target.fa /data/query.fa
NUCMER
>t1 q1 1000 500
101 300 1 200 5 5 0
10
-3
0
"""

DELTA_REV = """\
/data/target.fa /data/query.fa
NUCMER
>t1 q1 1000 500
101 300 200 1 5 5 0
0
"""


def test_nucmer_parser_forward():
    (rec,) = nucmer_delta_to_records(DELTA_FWD)
    assert (rec.target_name, rec.query_name) == ('t1', 'q1')
    assert (rec.target_len, rec.query_len) == (1000, 500)
    # 1-based inclusive -> 0-based half-open.
    assert (rec.target_start, rec.target_end) == (100, 300)
    assert (rec.query_start, rec.query_end) == (0, 200)
    assert rec.strand == '+'
    assert rec.alignment_block_len == 200
    assert rec.residue_matches == 195  # block - errors
    assert rec.mapping_quality == 255


def test_nucmer_parser_reverse():
    (rec,) = nucmer_delta_to_records(DELTA_REV)
    assert rec.strand == '-'
    # Query coords normalised to forward strand, 0-based half-open.
    assert (rec.query_start, rec.query_end) == (0, 200)
    assert rec.query_start < rec.query_end


def test_nucmer_parser_multi_contig_and_multi_alignment():
    text = (
        'ref.fa qry.fa\nNUCMER\n'
        '>t1 q1 1000 500\n'
        '1 100 1 100 0 0 0\n0\n'
        '201 400 301 500 2 2 0\n-5\n0\n'
        '>t2 q2 800 300\n'
        '50 149 300 201 1 1 0\n0\n'
    )
    recs = nucmer_delta_to_records(text)
    assert len(recs) == 3
    assert [(r.target_name, r.query_name) for r in recs] == [
        ('t1', 'q1'),
        ('t1', 'q1'),
        ('t2', 'q2'),
    ]
    assert [r.strand for r in recs] == ['+', '+', '-']
    assert (recs[2].query_start, recs[2].query_end) == (200, 300)
    assert (recs[2].target_start, recs[2].target_end) == (49, 149)


def test_nucmer_parser_indel_offsets_not_treated_as_alignments():
    # Large indel offsets (single ints) must be skipped, not parsed.
    recs = nucmer_delta_to_records(
        'ref qry\nNUCMER\n>t q 100 100\n1 50 1 50 0 0 0\n12\n-7\n3\n0\n'
    )
    assert len(recs) == 1


def test_nucmer_parser_bad_header():
    with pytest.raises(ValueError, match='sequence header'):
        nucmer_delta_to_records('ref qry\nNUCMER\n>t1 q1 1000\n1 10 1 10 0 0 0\n0\n')


def test_nucmer_parser_bad_coordinate_line():
    with pytest.raises(ValueError, match='coordinate line'):
        nucmer_delta_to_records('ref qry\nNUCMER\n>t1 q1 10 10\n1 10 1\n')


def test_nucmer_parser_non_numeric_coordinates():
    with pytest.raises(ValueError, match='coordinate line'):
        nucmer_delta_to_records('ref qry\nNUCMER\n>t1 q1 10 10\n1 x 1 10 0 0 0\n')


def test_nucmer_parser_empty():
    with pytest.raises(ValueError, match='no alignments'):
        nucmer_delta_to_records('ref qry\nNUCMER\n')


# ------------------------------------------------------ alignment_from_tool_output


def test_alignment_from_minimap2_paf():
    paf = 'q1\t100\t0\t50\t+\tt1\t200\t10\t60\t50\t50\t60\ttp:A:P\n'
    aln = alignment_from_tool_output('minimap2', paf)
    assert len(aln) == 1
    assert aln.records[0].query_name == 'q1'


def test_alignment_from_lastz_output():
    aln = alignment_from_tool_output('lastz', LASTZ_HEADER + '\n' + _lastz_line())
    assert len(aln) == 1
    assert aln.query_names == ['q1']
    assert aln.target_names == ['t1']


def test_alignment_from_nucmer_delta():
    aln = alignment_from_tool_output('nucmer', DELTA_FWD)
    assert len(aln) == 1


def test_alignment_from_unknown_tool():
    with pytest.raises(ValueError, match='Unknown biowasm tool'):
        alignment_from_tool_output('bowtie2', '')


# ---------------------------------------------------------------- cache keying


def test_cache_keying_by_method_params_and_digests():
    cache = SessionCache()
    query = parse_fasta_bytes(b'>q1\nACGTACGT\n')
    target = parse_fasta_bytes(b'>t1\nACGTACGTACGT\n')
    params = {'preset': 'asm20', 'k': 0, 'w': 0}
    aln = alignment_from_tool_output(
        'minimap2', 'q1\t8\t0\t8\t+\tt1\t12\t0\t8\t8\t8\t60\n'
    )

    assert cache.get_paf('minimap2', params, query.digest, target.digest) is None
    cache.put_paf('minimap2', params, aln, query.digest, target.digest)
    hit = cache.get_paf('minimap2', params, query.digest, target.digest)
    assert hit is aln
    # Param order must not matter.
    reordered = {'w': 0, 'k': 0, 'preset': 'asm20'}
    assert cache.get_paf('minimap2', reordered, query.digest, target.digest) is aln
    # Different params, method, or inputs miss.
    assert (
        cache.get_paf('minimap2', {**params, 'k': 19}, query.digest, target.digest)
        is None
    )
    assert cache.get_paf('lastz', params, query.digest, target.digest) is None
    assert cache.get_paf('minimap2', params, target.digest, query.digest) is None


def test_cache_roundtrip_via_build_args_params():
    # The exact params dicts used by the app (input to build_tool_args)
    # must be usable as cache keys for every biowasm tool.
    cache = SessionCache()
    q = parse_fasta_bytes(b'>q\nAAAA\n')
    t = parse_fasta_bytes(b'>t\nTTTT\n')
    all_params = {
        'minimap2': {'preset': 'asm5', 'k': 19, 'w': 10},
        'lastz': {'step': 5, 'gapped': True, 'notransition': False},
        'nucmer': {'l': 20, 'c': 65, 'maxmatch': True},
    }
    aln = alignment_from_tool_output('nucmer', DELTA_FWD)
    for tool, params in all_params.items():
        build_tool_args(tool, params)  # must not raise
        cache.put_paf(tool, params, aln, q.digest, t.digest)
        assert cache.get_paf(tool, params, q.digest, t.digest) is aln
