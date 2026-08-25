"""Tests for the browser app's GenBank input (app/core/genbank.py)."""

import gzip
from pathlib import Path
import sys

import pytest

APP_DIR = Path(__file__).resolve().parent.parent / 'app'
sys.path.insert(0, str(APP_DIR))

from core.genbank import (  # noqa: E402
    parse_genbank_bytes,
    parse_location,
)

SAMPLE = (Path(__file__).resolve().parent / 'data' / 'sample.gbk').read_bytes()


@pytest.fixture(scope='module')
def parsed():
    return parse_genbank_bytes(SAMPLE)


# ------------------------------------------------------------- locations


@pytest.mark.parametrize(
    ('text', 'parts', 'strand'),
    [
        ('101..400', ((100, 400),), '+'),
        ('complement(101..400)', ((100, 400),), '-'),
        ('join(101..200,301..400)', ((100, 200), (300, 400)), '+'),
        ('complement(join(801..830,861..900))', ((800, 830), (860, 900)), '-'),
        ('order(1..10,20..30)', ((0, 10), (19, 30)), '+'),
        ('950', ((949, 950),), '+'),
        ('<501..>700', ((500, 700),), '-'),  # strand from the caller's default
        ('750^751', ((749, 750),), '+'),
    ],
)
def test_parse_location_forms(text, parts, strand):
    loc = parse_location(text)
    assert loc.parts == parts
    if 'complement' in text:
        assert loc.strand == '-'
    else:
        assert loc.strand == '+'


def test_double_complement_returns_to_plus_strand():
    assert parse_location('complement(complement(1001..1100))').strand == '+'


def test_fuzzy_and_ordered_flags():
    assert parse_location('<501..>700').fuzzy
    assert not parse_location('501..700').fuzzy
    assert parse_location('order(1..10,20..30)').ordered
    assert not parse_location('join(1..10,20..30)').ordered


def test_remote_reference_selects_nothing():
    """A location in another record yields no parts -- but is not an error.

    The caller skips such features; aborting the whole upload over one
    cross-record reference would be hostile.
    """
    assert parse_location('NC_999999.1:10..20').parts == ()


@pytest.mark.parametrize('bad', ['', 'abc..def', '100..50', 'join()'])
def test_bad_locations_raise(bad):
    with pytest.raises(ValueError):
        parse_location(bad)


# --------------------------------------------------------------- records


def test_record_names_prefer_version_then_accession(parsed):
    # SCAF01 has VERSION SCAF01.2; SCAF02 has only ACCESSION.
    assert parsed.fasta.names == ['SCAF01.2', 'SCAF02']


def test_origin_sequence_stripped_of_offsets_and_spaces(parsed):
    seqs = dict(parsed.fasta.records)
    assert seqs['SCAF01.2'] == 'ACGT' * 300  # matches the declared 1200 bp
    assert seqs['SCAF02'] == ('TTTTGGGGCCCCAAAA' * 4)[:60]


def test_length_mismatch_warns_but_keeps_the_origin_sequence(caplog):
    """A hand-edited or truncated file is still worth plotting."""
    truncated = SAMPLE.replace(b'1200 bp', b'9999 bp')
    with caplog.at_level('WARNING'):
        out = parse_genbank_bytes(truncated)
    assert len(dict(out.fasta.records)['SCAF01.2']) == 1200
    assert any('declares 9999 bp' in r.getMessage() for r in caplog.records)


def test_digest_is_of_the_raw_upload(parsed):
    from core.fasta import content_digest

    assert parsed.fasta.digest == content_digest(SAMPLE)


def test_gzip_roundtrip():
    gz = parse_genbank_bytes(gzip.compress(SAMPLE))
    assert gz.fasta.records == parse_genbank_bytes(SAMPLE).fasta.records
    assert gz.fasta.digest != parse_genbank_bytes(SAMPLE).fasta.digest


def test_crlf_line_endings():
    crlf = SAMPLE.replace(b'\n', b'\r\n')
    assert parse_genbank_bytes(crlf).fasta.names == ['SCAF01.2', 'SCAF02']


@pytest.mark.parametrize(
    ('data', 'match'),
    [
        (b'', 'No GenBank records'),
        (b'DEFINITION nothing\n//\n', 'No GenBank records'),
        (b'\x1f\x8btruncated', 'gzip'),
        (b'\xff\xfe\x00bad', 'not UTF-8'),
    ],
)
def test_parse_errors(data, match):
    with pytest.raises(ValueError, match=match):
        parse_genbank_bytes(data)


def test_duplicate_record_names_rejected():
    doubled = SAMPLE + SAMPLE
    with pytest.raises(ValueError, match='Duplicate record name'):
        parse_genbank_bytes(doubled)


def test_record_without_origin_rejected():
    trimmed = SAMPLE.split(b'ORIGIN')[0] + b'//\n'
    with pytest.raises(ValueError, match='no ORIGIN sequence'):
        parse_genbank_bytes(trimmed)


# -------------------------------------------------------------- features


def _annotation(parsed):
    from rusty_dot.annotation import GffAnnotation

    return GffAnnotation.from_text(parsed.gff_text)


def test_emitted_gff_roundtrips_through_the_library(parsed):
    ann = _annotation(parsed)
    types = ann.feature_types()
    assert 'CDS' in types and 'gene' in types and 'tRNA' in types
    # source is dropped by default: it spans the record and would shade
    # the whole diagonal.
    assert 'source' not in types


def test_include_source_keeps_the_whole_record_feature():
    with_source = parse_genbank_bytes(SAMPLE, include_source=True)
    assert 'source' in _annotation(with_source).feature_types()


def test_coordinates_survive_the_gff_round_trip(parsed):
    ann = _annotation(parsed)
    gene = next(
        r
        for r in ann.records
        if r.feature_type == 'gene' and r.attr_dict().get('locus_tag') == 'TST_0001'
    )
    assert (gene.start, gene.end) == (100, 400)  # 0-based half-open
    assert gene.strand == '-'


def test_join_parts_group_into_one_feature(parsed):
    """The synthetic Parent matters: iter_groups needs both ID and Parent.

    Without it a join(...) CDS draws as unrelated fragments with no
    connector line.
    """
    ann = _annotation(parsed)
    groups = ann.iter_groups('SCAF01.2')
    cds_groups = [parts for _key, parts in groups if parts[0].feature_type == 'CDS']
    assert len(cds_groups) == 1
    assert [(p.start, p.end) for p in cds_groups[0]] == [(100, 200), (300, 400)]

    trna = [parts for _key, parts in groups if parts[0].feature_type == 'tRNA'][0]
    assert len(trna) == 2
    assert all(p.strand == '-' for p in trna)


def test_multiline_qualifiers_are_joined_with_spaces(parsed):
    ann = _annotation(parsed)
    cds = next(r for r in ann.records if r.feature_type == 'CDS')
    attrs = cds.attr_dict()
    assert attrs['product'] == 'ABC transporter, ATP-binding protein subunit A'
    assert attrs['note'] == 'first line of note second line of note'


def test_translation_qualifier_is_dropped(parsed):
    assert 'translation' not in parsed.gff_text
    assert 'MKKLLTAAAAV' not in parsed.gff_text


def test_reserved_characters_are_percent_encoded(parsed):
    """A note containing ';' and '=' must not corrupt column 9."""
    ann = _annotation(parsed)
    rep = next(r for r in ann.records if r.feature_type == 'repeat_region')
    assert rep.attr_dict()['note'] == 'partial repeat; with=reserved,chars'


def test_fuzzy_bounds_are_flagged_and_stripped(parsed):
    ann = _annotation(parsed)
    rep = next(r for r in ann.records if r.feature_type == 'repeat_region')
    assert (rep.start, rep.end) == (500, 700)
    assert rep.attr_dict().get('partial') == 'true'


def test_zero_width_and_single_base_features_stay_drawable(parsed):
    ann = _annotation(parsed)
    misc = next(r for r in ann.records if r.feature_type == 'misc_feature')
    assert misc.end > misc.start  # 750^751 must not collapse to zero width
    mrna = next(r for r in ann.records if r.feature_type == 'mRNA')
    assert (mrna.start, mrna.end) == (949, 950)


def test_bare_qualifier_without_value_is_tolerated(parsed):
    ann = _annotation(parsed)
    dbl = next(r for r in ann.records if r.attr_dict().get('gene') == 'dblComp')
    assert dbl.strand == '+'  # complement(complement(...))


def test_remote_reference_feature_is_skipped_not_fatal(parsed):
    ann = _annotation(parsed)
    assert 'misc_binding' not in ann.feature_types()


def test_seqnames_match_the_fasta_record_names(parsed):
    """Otherwise every feature silently vanishes from the plot."""
    ann = _annotation(parsed)
    assert {r.seqname for r in ann.records} <= set(parsed.fasta.names)


# ----------------------------------------------------------- end-to-end


def test_genbank_input_plots_like_a_fasta_upload(parsed):
    """The exact call path the app uses: GenBank -> index + annotation."""
    from rusty_dot import DotPlotter, SequenceIndex

    idx = SequenceIndex(k=11)
    for name, seq in parsed.fasta.records:
        idx.add_sequence(name, seq)
    ann = _annotation(parsed)
    name = parsed.fasta.names[0]
    fig = DotPlotter(idx).plot(
        query_names=[name],
        target_names=[name],
        annotation=ann,
        annotation_query=ann,
        annotation_target=ann,
        annotation_tracks=True,
    )
    assert fig.axes
