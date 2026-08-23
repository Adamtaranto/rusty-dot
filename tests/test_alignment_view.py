"""Tests for the CIGAR-based aligned-sequence rendering."""

import pytest

from rusty_dot.alignment_view import aligned_text
from rusty_dot.paf_io import PafRecord


def _rec(cigar, *, q_start=0, q_end=10, t_start=0, t_end=10, strand='+', q_len=10):
    line = (
        f'q\t{q_len}\t{q_start}\t{q_end}\t{strand}'
        f'\tt\t100\t{t_start}\t{t_end}\t10\t10\t60\tcg:Z:{cigar}'
    )
    return PafRecord.from_line(line)


class TestAlignedText:
    def test_perfect_match(self):
        out = aligned_text(_rec('10M'), 'ACGTACGTAC', 'ACGTACGTACTTTT')
        assert out['text'] == 'ACGTACGTAC\n||||||||||\nACGTACGTAC'
        assert out['cols'] == 10
        assert not out['truncated']

    def test_m_op_marks_mismatches(self):
        out = aligned_text(_rec('10M'), 'ACGTACGTAC', 'ACGAACGTAC')
        q_line, m_line, t_line = out['text'].split('\n')
        assert m_line == '||| ||||||'
        assert q_line[3] == 'T' and t_line[3] == 'A'

    def test_eqx_ops(self):
        out = aligned_text(_rec('4=2X4='), 'ACGTTTACGT', 'ACGTAAACGT')
        assert out['text'].split('\n')[1] == '||||  ||||'

    def test_insertion_gaps_target(self):
        rec = _rec('4M2I4M', q_end=10, t_end=8)
        out = aligned_text(rec, 'ACGTTTACGT', 'ACGTACGT')
        q_line, m_line, t_line = out['text'].split('\n')
        assert q_line == 'ACGTTTACGT'
        assert t_line == 'ACGT--ACGT'
        assert m_line == '||||  ||||'

    def test_deletion_gaps_query(self):
        rec = _rec('4M2D4M', q_end=8, t_end=10)
        out = aligned_text(rec, 'ACGTACGT', 'ACGTTTACGT')
        q_line, m_line, t_line = out['text'].split('\n')
        assert q_line == 'ACGT--ACGT'
        assert t_line == 'ACGTTTACGT'

    def test_minus_strand_revcomps_query(self):
        # Query slice AAAACCCC revcomped -> GGGGTTTT aligns to the target.
        rec = _rec('8M', q_end=8, t_end=8, q_len=8, strand='-')
        out = aligned_text(rec, 'AAAACCCC', 'GGGGTTTT')
        q_line, m_line, t_line = out['text'].split('\n')
        assert q_line == 'GGGGTTTT'
        assert m_line == '||||||||'

    def test_wrapping(self):
        rec = _rec('100M', q_end=100, t_end=100, q_len=100)
        out = aligned_text(rec, 'A' * 100, 'A' * 100, width=40)
        blocks = out['text'].split('\n\n')
        assert len(blocks) == 3  # 40 + 40 + 20 columns
        assert blocks[0] == 'A' * 40 + '\n' + '|' * 40 + '\n' + 'A' * 40
        assert blocks[2].split('\n')[0] == 'A' * 20

    def test_truncation(self):
        rec = _rec('100M', q_end=100, t_end=100, q_len=100)
        out = aligned_text(rec, 'A' * 100, 'C' * 100, max_cols=30)
        assert out['truncated']
        assert out['cols'] == 30
        assert '[truncated at 30 alignment columns]' in out['text']

    def test_no_cigar_raises(self):
        rec = PafRecord.from_line('q\t10\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60')
        with pytest.raises(ValueError, match='no CIGAR'):
            aligned_text(rec, 'A' * 10, 'A' * 10)
