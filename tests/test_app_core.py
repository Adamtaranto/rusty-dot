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
