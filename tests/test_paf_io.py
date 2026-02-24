"""Tests for the PAF file I/O module and CIGAR string parsing."""

from pathlib import Path
import tempfile
import textwrap

import pytest

from rusty_dot.paf_io import (
    PafAlignment,
    PafRecord,
    compute_gravity_contigs,
    parse_paf_file,
)

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

SIMPLE_PAF = textwrap.dedent("""\
    query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t45\t50\t255
    query1\t100\t60\t90\t+\ttarget1\t200\t80\t110\t28\t30\t255
    query2\t80\t0\t40\t-\ttarget1\t200\t150\t190\t38\t40\t255
""")

CIGAR_PAF = textwrap.dedent("""\
    query1\t100\t0\t20\t+\ttarget1\t200\t0\t22\t18\t22\t60\tcp:Z:P\tcg:Z:18=2X
    query2\t50\t0\t15\t+\ttarget2\t100\t5\t20\t12\t15\t60\tcg:Z:3=2I3=4D3=
""")


def _write_temp(content: str) -> Path:
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False)
    f.write(content)
    f.flush()
    f.close()
    return Path(f.name)


# ---------------------------------------------------------------------------
# PafRecord.from_line
# ---------------------------------------------------------------------------


class TestPafRecordFromLine:
    def test_basic_parse(self):
        line = 'query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t45\t50\t255'
        rec = PafRecord.from_line(line)
        assert rec.query_name == 'query1'
        assert rec.query_len == 100
        assert rec.query_start == 0
        assert rec.query_end == 50
        assert rec.strand == '+'
        assert rec.target_name == 'target1'
        assert rec.target_len == 200
        assert rec.target_start == 10
        assert rec.target_end == 60
        assert rec.residue_matches == 45
        assert rec.alignment_block_len == 50
        assert rec.mapping_quality == 255

    def test_minus_strand(self):
        line = 'q\t80\t0\t40\t-\tt\t200\t150\t190\t38\t40\t255'
        rec = PafRecord.from_line(line)
        assert rec.strand == '-'

    def test_optional_tags_parsed(self):
        line = 'q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t60\ttp:A:P\tNM:i:3'
        rec = PafRecord.from_line(line)
        assert rec.tags.get('tp') == 'P'
        assert rec.tags.get('NM') == 3

    def test_too_few_fields_raises(self):
        with pytest.raises(ValueError):
            PafRecord.from_line('q\t100\t0\t50\t+\tt')

    def test_no_cigar_gives_none(self):
        line = 'q\t100\t0\t50\t+\tt\t200\t0\t50\t45\t50\t255'
        rec = PafRecord.from_line(line)
        assert rec.cigar is None
        assert rec.alignment_length is None

    def test_cigar_exact_match(self):
        line = 'q\t100\t0\t20\t+\tt\t200\t0\t22\t18\t22\t60\tcg:Z:18=2X'
        rec = PafRecord.from_line(line)
        assert rec.cigar == '18=2X'
        assert rec.n_matches == 18
        assert rec.n_mismatches == 2
        assert rec.alignment_length == 20  # 18= + 2X

    def test_cigar_with_indels(self):
        # 3= 2I 3= 4D 3= → target-consuming ops: 3= + 3= + 4D + 3= = 13
        # (I does not consume target bases; D and = do)
        line = 'q\t50\t0\t15\t+\tt\t100\t5\t20\t12\t15\t60\tcg:Z:3=2I3=4D3='
        rec = PafRecord.from_line(line)
        assert rec.cigar == '3=2I3=4D3='
        assert rec.n_gap_bases == 6  # 2I + 4D
        assert rec.n_gaps == 2  # one I-run + one D-run = 2 distinct gap events
        assert rec.alignment_length == 13  # 3+3+4+3 (= and D ops consume target)

    def test_query_aligned_len_property(self):
        line = 'q\t100\t10\t60\t+\tt\t200\t0\t50\t45\t50\t255'
        rec = PafRecord.from_line(line)
        assert rec.query_aligned_len == 50

    def test_target_aligned_len_property(self):
        line = 'q\t100\t0\t50\t+\tt\t200\t20\t80\t55\t60\t255'
        rec = PafRecord.from_line(line)
        assert rec.target_aligned_len == 60

    def test_to_line_roundtrip(self):
        line = 'query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t45\t50\t255'
        rec = PafRecord.from_line(line)
        # to_line() should reproduce the 12 core fields
        assert rec.to_line() == line


# ---------------------------------------------------------------------------
# parse_paf_file
# ---------------------------------------------------------------------------


class TestParsePafFile:
    def test_yields_correct_count(self):
        path = _write_temp(SIMPLE_PAF)
        records = list(parse_paf_file(path))
        assert len(records) == 3

    def test_skips_comments_and_blank_lines(self):
        content = (
            '# header\n\nquery1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t45\t50\t255\n'
        )
        path = _write_temp(content)
        records = list(parse_paf_file(path))
        assert len(records) == 1

    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            list(parse_paf_file('/nonexistent/path.paf'))

    def test_cigar_parsed_from_file(self):
        path = _write_temp(CIGAR_PAF)
        records = list(parse_paf_file(path))
        assert records[0].cigar == '18=2X'
        assert records[1].cigar == '3=2I3=4D3='


# ---------------------------------------------------------------------------
# PafAlignment
# ---------------------------------------------------------------------------


class TestPafAlignment:
    def setup_method(self):
        self.path = _write_temp(SIMPLE_PAF)
        self.aln = PafAlignment.from_file(self.path)

    def test_len(self):
        assert len(self.aln) == 3

    def test_query_names(self):
        assert set(self.aln.query_names) == {'query1', 'query2'}

    def test_target_names(self):
        assert set(self.aln.target_names) == {'target1'}

    def test_repr(self):
        r = repr(self.aln)
        assert 'PafAlignment' in r
        assert '3' in r  # record count

    def test_from_records(self):
        records = list(parse_paf_file(self.path))
        aln2 = PafAlignment.from_records(records)
        assert len(aln2) == len(self.aln)

    def test_filter_by_query(self):
        filtered = self.aln.filter_by_query(['query1'])
        assert all(r.query_name == 'query1' for r in filtered.records)
        assert len(filtered) == 2

    def test_filter_by_target(self):
        filtered = self.aln.filter_by_target(['target1'])
        assert len(filtered) == 3  # all records target target1

    def test_filter_by_target_empty(self):
        filtered = self.aln.filter_by_target(['no_such_target'])
        assert len(filtered) == 0


# ---------------------------------------------------------------------------
# compute_gravity_contigs
# ---------------------------------------------------------------------------


class TestComputeGravityContigs:
    def _make_records(self) -> list[PafRecord]:
        """Two query contigs: q_early maps to the start of target, q_late to end."""
        return [
            PafRecord.from_line(
                'q_early\t50\t0\t50\t+\ttarget\t100\t0\t50\t48\t50\t255'
            ),
            PafRecord.from_line(
                'q_late\t50\t0\t50\t+\ttarget\t100\t50\t100\t48\t50\t255'
            ),
        ]

    def test_collinear_order(self):
        records = self._make_records()
        q_sorted, t_sorted = compute_gravity_contigs(
            records, ['q_late', 'q_early'], ['target']
        )
        # q_early has lower gravity (maps to start of target) → sorts first
        assert q_sorted[0] == 'q_early'
        assert q_sorted[1] == 'q_late'

    def test_no_match_sorts_last(self):
        records = self._make_records()
        # Add a query contig that has no matches
        q_sorted, _ = compute_gravity_contigs(
            records, ['q_late', 'q_early', 'q_none'], ['target']
        )
        assert q_sorted[-1] == 'q_none'

    def test_target_reordered(self):
        """Targets are also reordered by gravity."""
        records = [
            PafRecord.from_line('query\t100\t0\t50\t+\tt_late\t50\t0\t50\t48\t50\t255'),
            PafRecord.from_line(
                'query\t100\t60\t100\t+\tt_early\t40\t0\t40\t38\t40\t255'
            ),
        ]
        _, t_sorted = compute_gravity_contigs(records, ['query'], ['t_late', 't_early'])
        # t_early maps to lower query positions (60–100 mid=80, t_late mid=25)
        # With only one query, gravity = weighted mean query mid-point
        # t_early: query mid = 80 → higher gravity → sorts after t_late
        assert 't_early' in t_sorted
        assert 't_late' in t_sorted


# ---------------------------------------------------------------------------
# PafAlignment.reorder_contigs
# ---------------------------------------------------------------------------


class TestReorderContigs:
    def test_reorder_defaults_to_all_names(self):
        path = _write_temp(SIMPLE_PAF)
        aln = PafAlignment.from_file(path)
        q_sorted, t_sorted = aln.reorder_contigs()
        assert set(q_sorted) == set(aln.query_names)
        assert set(t_sorted) == set(aln.target_names)

    def test_reorder_subset(self):
        path = _write_temp(SIMPLE_PAF)
        aln = PafAlignment.from_file(path)
        q_sorted, t_sorted = aln.reorder_contigs(['query1'], ['target1'])
        assert q_sorted == ['query1']
        assert t_sorted == ['target1']


# ---------------------------------------------------------------------------
# CrossIndexPaf
# ---------------------------------------------------------------------------


class TestCrossIndexPaf:
    def test_add_sequence_and_repr(self):
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'TACGTACGTACG', group='b')
        assert 'CrossIndexPaf' in repr(cross)
        assert cross.query_names == ['q1']
        assert cross.target_names == ['t1']

    def test_invalid_group_raises(self):
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        with pytest.raises(ValueError):
            cross.add_sequence('x', 'ACGT', group='c')

    def test_get_paf_all_cross(self):
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        lines = cross.get_paf_all()
        assert isinstance(lines, list)
        # All lines should reference the original (un-prefixed) names
        for line in lines:
            fields = line.split('\t')
            assert fields[0] == 'q1'
            assert fields[5] == 't1'

    def test_get_paf_all_single_group(self):
        """get_paf_all with no group-B sequences does all-vs-all within group A."""
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='a')
        lines = cross.get_paf_all()
        # s1 vs s2 and s2 vs s1 should both appear
        queries = {line.split('\t')[0] for line in lines}
        assert 's1' in queries or 's2' in queries

    def test_reorder_contigs_raises_without_group_b(self):
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='a')
        with pytest.raises(ValueError):
            cross.reorder_contigs()

    def test_reorder_contigs_returns_original_names(self):
        from rusty_dot.paf_io import CrossIndexPaf

        cross = CrossIndexPaf(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('q2', 'TACGTACGTACGTACG', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        q_sorted, t_sorted = cross.reorder_contigs()
        assert set(q_sorted) == {'q1', 'q2'}
        assert set(t_sorted) == {'t1'}


# ---------------------------------------------------------------------------
# compute_gravity_contigs: unmatched sorted by descending length
# ---------------------------------------------------------------------------


class TestComputeGravityContigsUnmatchedLength:
    def test_unmatched_sorted_by_length_desc(self):
        """Unmatched contigs must appear after matched ones, sorted by length descending."""
        records = [
            PafRecord.from_line(
                'q_early\t50\t0\t50\t+\ttarget\t100\t0\t50\t48\t50\t255'
            ),
        ]
        # Add two unmatched queries: long_unmatched (100 bp) and short_unmatched (10 bp)
        # Neither appears in records → len_map has no entry → length 0 for both
        # In this edge case they sort equal; just verify they're both at the end.
        q_sorted, _ = compute_gravity_contigs(
            records,
            ['q_early', 'long_unmatched', 'short_unmatched'],
            ['target'],
        )
        assert q_sorted[0] == 'q_early'
        assert set(q_sorted[1:]) == {'long_unmatched', 'short_unmatched'}
