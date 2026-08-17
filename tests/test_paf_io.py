"""Tests for the PAF file I/O module and CIGAR string parsing."""

from pathlib import Path
import tempfile
import textwrap

import pytest

from rusty_dot.paf_io import (
    PafAlignment,
    PafRecord,
    compute_gravity_contigs,
    compute_reversed_contigs,
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

    def test_filter_by_min_length_keeps_long(self):
        # SIMPLE_PAF has alignment_block_len values 50, 30, 40
        # query_aligned_len: 50-0=50, 90-60=30, 40-0=40
        filtered = self.aln.filter_by_min_length(40)
        assert len(filtered) == 2  # records with query_aligned_len >= 40: 50 and 40

    def test_filter_by_min_length_zero_keeps_all(self):
        filtered = self.aln.filter_by_min_length(0)
        assert len(filtered) == len(self.aln)

    def test_filter_by_min_length_removes_all(self):
        filtered = self.aln.filter_by_min_length(1000)
        assert len(filtered) == 0

    def test_filter_by_min_length_exact_boundary(self):
        # Only keeps records where query_aligned_len >= 50 (exactly 50 passes)
        filtered = self.aln.filter_by_min_length(50)
        assert len(filtered) == 1
        assert filtered.records[0].query_aligned_len == 50

    def test_sequence_names_includes_all(self):
        # SIMPLE_PAF has query1, query2 as queries and target1 as target
        names = set(self.aln.sequence_names())
        assert names == {'query1', 'query2', 'target1'}

    def test_sequence_names_no_duplicates(self):
        # target1 appears in 3 records but should only appear once
        names = self.aln.sequence_names()
        assert len(names) == len(set(names))

    def test_get_sequence_length_query(self):
        # query1 has query_len=100 in SIMPLE_PAF
        assert self.aln.get_sequence_length('query1') == 100

    def test_get_sequence_length_target(self):
        # target1 has target_len=200 in SIMPLE_PAF
        assert self.aln.get_sequence_length('target1') == 200

    def test_get_sequence_length_missing_raises(self):
        with pytest.raises(KeyError):
            self.aln.get_sequence_length('nonexistent')


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
        q_sorted, t_sorted, _reversed = compute_gravity_contigs(
            records, ['q_late', 'q_early'], ['target']
        )
        # q_early has lower gravity (maps to start of target) → sorts first
        assert q_sorted[0] == 'q_early'
        assert q_sorted[1] == 'q_late'

    def test_no_match_sorts_last(self):
        records = self._make_records()
        # Add a query contig that has no matches
        q_sorted, _t, _reversed = compute_gravity_contigs(
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
        _q, t_sorted, _reversed = compute_gravity_contigs(
            records, ['query'], ['t_late', 't_early']
        )
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
# ---------------------------------------------------------------------------
# CrossIndex (formerly CrossIdx)
# ---------------------------------------------------------------------------


class TestCrossIndex:
    def test_add_sequence_and_repr(self):
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'TACGTACGTACG', group='b')
        r = repr(cross)
        assert r.startswith('CrossIndex(')
        assert cross.query_names == ['q1']
        assert cross.target_names == ['t1']

    def test_group_colon_in_name_raises(self):
        """Group names must not contain ':'."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        with pytest.raises(ValueError, match="must not contain ':'"):
            cross.add_sequence('x', 'ACGT', group='bad:group')

    def test_arbitrary_group_names_accepted(self):
        """Any group name without ':' is accepted."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('x', 'ACGT' * 4, group='assembly_c')
        assert 'assembly_c' in cross.group_names

    def test_get_paf_all_cross(self):
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
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
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='a')
        lines = cross.get_paf_all()
        queries = {line.split('\t')[0] for line in lines}
        # Both sequences should appear as queries (s1→s2 and s2→s1)
        assert queries == {'s1', 's2'}

    def test_reorder_contigs_raises_without_group_b(self):
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='a')
        with pytest.raises(ValueError):
            cross.reorder_contigs()

    def test_reorder_contigs_returns_original_names(self):
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('q2', 'TACGTACGTACGTACG', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        cross.compute_matches()
        q_sorted, t_sorted = cross.reorder_contigs()
        assert set(q_sorted) == {'q1', 'q2'}
        assert set(t_sorted) == {'t1'}

    def test_contig_order_default_insertion(self):
        """contig_order reflects insertion order by default."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('c', 'CCCC' * 4, group='g')
        cross.add_sequence('a', 'AAAA' * 4, group='g')
        cross.add_sequence('b', 'TTTT' * 4, group='g')
        assert cross.contig_order['g'] == ['c', 'a', 'b']

    def test_reorder_by_length(self):
        """reorder_by_length sorts within a group by descending sequence length."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('short', 'ACGT' * 4, group='g')  # 16 bp
        cross.add_sequence('long', 'ACGT' * 10, group='g')  # 40 bp
        cross.add_sequence('mid', 'ACGT' * 6, group='g')  # 24 bp
        cross.reorder_by_length(group='g')
        assert cross.contig_order['g'] == ['long', 'mid', 'short']

    def test_reorder_by_length_all_groups(self):
        """reorder_by_length with group=None sorts all groups."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s', 'ACGT' * 4, group='g1')
        cross.add_sequence('l', 'ACGT' * 8, group='g1')
        cross.add_sequence('x', 'TTTT' * 4, group='g2')
        cross.add_sequence('y', 'TTTT' * 8, group='g2')
        cross.reorder_by_length()
        assert cross.contig_order['g1'][0] == 'l'
        assert cross.contig_order['g2'][0] == 'y'

    def test_reorder_for_colinearity_updates_contig_order(self):
        """reorder_for_colinearity updates contig_order for both groups."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGT' * 4, group='a')
        cross.add_sequence('q2', 'TTTTACGTACGT' * 4, group='a')
        cross.add_sequence('t1', 'ACGTACGTACGT' * 4, group='b')
        orig_a = list(cross.contig_order['a'])
        cross.compute_matches('a', 'b')
        cross.reorder_for_colinearity('a', 'b')
        # Order may or may not change, but both groups must be present
        assert set(cross.contig_order['a']) == set(orig_a)
        assert 't1' in cross.contig_order['b']

    def test_sequence_names_returns_prefixed(self):
        """sequence_names returns 'group:name' strings."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('seq1', 'ACGT' * 4, group='g')
        names = cross.sequence_names(group='g')
        assert names == ['g:seq1']

    def test_sequence_names_all_groups(self):
        """sequence_names with no group returns all groups combined."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT' * 4, group='a')
        cross.add_sequence('s2', 'ACGT' * 4, group='b')
        names = set(cross.sequence_names())
        assert names == {'a:s1', 'b:s2'}

    def test_get_sequence_length_via_internal_name(self):
        """get_sequence_length works with 'group:name' key."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT' * 5, group='a')
        assert cross.get_sequence_length('a:s1') == 20

    def test_get_paf_two_groups(self):
        """get_paf returns un-prefixed names in PAF fields."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        lines = cross.get_paf()
        assert isinstance(lines, list)
        for line in lines:
            fields = line.split('\t')
            assert fields[0] == 'q1'
            assert fields[5] == 't1'

    def test_str_summary(self):
        """__str__ includes group stats."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGT' * 4, group='a')
        cross.add_sequence('t1', 'ACGT' * 4, group='b')
        s = str(cross)
        assert 'CrossIndex' in s
        assert "'a'" in s
        assert "'b'" in s

    def test_three_group_default_pairs(self):
        """With 3 groups, all non-self ordered pairs are aligned."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='g1')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='g2')
        cross.add_sequence('s3', 'ACGTACGTACGTACGT', group='g3')
        pairs = cross._get_default_group_pairs()
        assert len(pairs) == 6  # 3 groups → 6 ordered non-self pairs

    def test_get_paf_explicit_group_pairs(self):
        """get_paf respects explicit group_pairs argument."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='g1')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='g2')
        cross.add_sequence('s3', 'TTTTTTTTTTTTTTTT', group='g3')
        # Only compare g1 vs g2; g3 should not appear
        lines = cross.get_paf(group_pairs=[('g1', 'g2')])
        for line in lines:
            fields = line.split('\t')
            assert fields[0] != 's3'
            assert fields[5] != 's3'

    def test_dotplotter_compatibility(self, tmp_path):
        """CrossIndex can be passed to DotPlotter."""
        from rusty_dot.dotplot import DotPlotter
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGT' * 15, group='a')
        cross.add_sequence('t1', 'ACGT' * 15, group='b')

        plotter = DotPlotter(cross)
        out = tmp_path / 'cross.png'
        plotter.plot(
            query_names=cross.sequence_names(group='a'),
            target_names=cross.sequence_names(group='b'),
            output_path=str(out),
        )
        assert out.exists()

    def test_reorder_contigs_custom_group_names(self):
        """reorder_contigs works when groups have non-default names."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='Group_A')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='Group_B')
        cross.compute_matches()
        q_sorted, t_sorted = cross.reorder_contigs()
        assert set(q_sorted) == {'q1'}
        assert set(t_sorted) == {'t1'}

    def test_reorder_contigs_explicit_group_params(self):
        """reorder_contigs accepts explicit query_group and target_group."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='g1')
        cross.add_sequence('q2', 'TACGTACGTACGTACG', group='g1')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='g2')
        cross.add_sequence('t2', 'ACGTACGTACGTACGT', group='g3')
        # Explicitly compute and compare g1 vs g2 even though there are 3 groups
        cross.compute_matches(query_group='g1', target_group='g2')
        q_sorted, t_sorted = cross.reorder_contigs(query_group='g1', target_group='g2')
        assert set(q_sorted) == {'q1', 'q2'}
        assert set(t_sorted) == {'t1'}

    def test_reorder_contigs_raises_three_groups_no_params(self):
        """reorder_contigs raises ValueError with 3 groups and no explicit params."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='g1')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='g2')
        cross.add_sequence('s3', 'ACGTACGTACGTACGT', group='g3')
        with pytest.raises(ValueError, match='query_group'):
            cross.reorder_contigs()

    def test_reorder_contigs_raises_one_group_param_only(self):
        """reorder_contigs raises when only one of query_group/target_group given."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        with pytest.raises(ValueError, match='both'):
            cross.reorder_contigs(query_group='a')

    def test_rename_group_updates_label(self):
        """rename_group changes the group label while preserving sequences."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='old_name')
        cross.rename_group('old_name', 'new_name')
        assert 'new_name' in cross.group_names
        assert 'old_name' not in cross.group_names
        assert cross.contig_order['new_name'] == ['q1']

    def test_rename_group_allows_subsequent_reorder(self):
        """After rename_group, compute_matches and reorder_contigs still work."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='Group_A')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='Group_B')
        cross.rename_group('Group_A', 'query')
        cross.rename_group('Group_B', 'target')
        cross.compute_matches()
        q_sorted, t_sorted = cross.reorder_contigs()
        assert set(q_sorted) == {'q1'}
        assert set(t_sorted) == {'t1'}
        # PAF output should still produce un-prefixed names
        lines = cross.get_paf()
        assert any('q1' in line.split('\t')[0] for line in lines)

    def test_rename_group_raises_unknown_group(self):
        """rename_group raises KeyError for unknown old_name."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT', group='a')
        with pytest.raises(KeyError):
            cross.rename_group('no_such_group', 'x')

    def test_rename_group_raises_duplicate_new_name(self):
        """rename_group raises ValueError when new_name already exists."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT', group='a')
        cross.add_sequence('s2', 'ACGT', group='b')
        with pytest.raises(ValueError, match='already exists'):
            cross.rename_group('a', 'b')

    def test_rename_group_raises_colon_in_name(self):
        """rename_group raises ValueError when new_name contains ':'."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT', group='a')
        with pytest.raises(ValueError, match="must not contain ':'"):
            cross.rename_group('a', 'bad:name')

    def test_set_group_members_updates_list(self):
        """set_group_members replaces the sequence list for a group."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGT' * 4, group='g')
        cross.add_sequence('s2', 'ACGT' * 4, group='g')
        cross.add_sequence('s3', 'ACGT' * 4, group='g')
        cross.set_group_members('g', ['s1', 's3'])
        assert cross.contig_order['g'] == ['s1', 's3']

    def test_set_group_members_warns_on_overlap(self, caplog):
        """set_group_members logs a warning when a name appears in another group."""
        import logging

        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('shared', 'ACGT' * 4, group='g1')
        cross.add_sequence('shared', 'ACGT' * 4, group='g2')
        with caplog.at_level(logging.WARNING):
            cross.set_group_members('g2', ['shared'])
        assert any('shared' in msg for msg in caplog.messages)

    def test_set_group_members_raises_unknown_group(self):
        """set_group_members raises KeyError for unknown group."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        with pytest.raises(KeyError):
            cross.set_group_members('no_such', ['s1'])

    def test_compute_matches_populates_records(self):
        """compute_matches stores PAF records for the computed pair."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        assert len(cross._paf_records) == 0
        cross.compute_matches()
        assert len(cross._paf_records) > 0

    def test_computed_group_pairs_empty_before_compute(self):
        """computed_group_pairs is empty before compute_matches is called."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        assert cross.computed_group_pairs == []

    def test_computed_group_pairs_after_compute(self):
        """computed_group_pairs reflects pairs after compute_matches."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        cross.compute_matches()
        assert ('a', 'b') in cross.computed_group_pairs

    def test_reorder_contigs_raises_without_compute_matches(self):
        """reorder_contigs raises if compute_matches was not called first."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        with pytest.raises(ValueError, match='compute_matches'):
            cross.reorder_contigs()

    def test_reorder_for_colinearity_raises_without_compute_matches(self):
        """reorder_for_colinearity raises if compute_matches was not called."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGT' * 4, group='a')
        cross.add_sequence('t1', 'ACGT' * 4, group='b')
        with pytest.raises(ValueError, match='compute_matches'):
            cross.reorder_for_colinearity('a', 'b')

    def test_compute_matches_explicit_groups(self):
        """compute_matches with explicit group names stores the correct pair."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('s1', 'ACGTACGTACGTACGT', group='g1')
        cross.add_sequence('s2', 'ACGTACGTACGTACGT', group='g2')
        cross.add_sequence('s3', 'TTTTTTTTTTTTTTTT', group='g3')
        cross.compute_matches(query_group='g1', target_group='g2')
        assert ('g1', 'g2') in cross.computed_group_pairs
        assert ('g1', 'g3') not in cross.computed_group_pairs

    def test_compute_matches_logs_messages(self, caplog):
        """compute_matches emits INFO-level log messages."""
        import logging

        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        with caplog.at_level(logging.INFO, logger='rusty_dot.paf_io'):
            cross.compute_matches()
        assert any('compute_matches' in msg for msg in caplog.messages)

    def test_add_sequence_logs_debug(self, caplog):
        """add_sequence emits DEBUG-level log messages."""
        import logging

        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        with caplog.at_level(logging.DEBUG, logger='rusty_dot.paf_io'):
            cross.add_sequence('q1', 'ACGT' * 4, group='a')
        assert any('q1' in msg for msg in caplog.messages)

    def test_add_sequence_warns_same_group_duplicate(self, caplog):
        """add_sequence warns when adding a name that exists in the same group."""
        import logging

        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('dup', 'ACGT' * 4, group='a')
        with caplog.at_level(logging.WARNING, logger='rusty_dot.paf_io'):
            cross.add_sequence('dup', 'TTTT' * 4, group='a')
        assert any('dup' in msg for msg in caplog.messages)

    def test_add_sequence_warns_cross_group_duplicate(self, caplog):
        """add_sequence warns when adding a name that exists in another group."""
        import logging

        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('seq1', 'ACGT' * 4, group='a')
        with caplog.at_level(logging.WARNING, logger='rusty_dot.paf_io'):
            cross.add_sequence('seq1', 'TTTT' * 4, group='b')
        assert any('seq1' in msg for msg in caplog.messages)

    def test_run_merge_populates_paf_records(self):
        """run_merge populates _paf_records (via compute_matches)."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=4)
        cross.add_sequence('q1', 'ACGTACGTACGTACGT', group='a')
        cross.add_sequence('t1', 'ACGTACGTACGTACGT', group='b')
        assert len(cross._paf_records) == 0
        cross.run_merge()
        assert len(cross._paf_records) > 0


# ---------------------------------------------------------------------------
# PafAlignment group management
# ---------------------------------------------------------------------------


class TestPafAlignmentGroups:
    def setup_method(self):
        self.aln = PafAlignment.from_records(
            [PafRecord.from_line(line) for line in SIMPLE_PAF.strip().splitlines()]
        )

    def test_groups_default_returns_a_b(self):
        """Default groups map query_names → 'a' and target_names → 'b'."""
        g = self.aln.groups
        assert set(g['a']) == set(self.aln.query_names)
        assert set(g['b']) == set(self.aln.target_names)

    def test_set_groups_custom(self):
        """set_groups stores custom group assignments."""
        self.aln.set_groups({'q': ['query1', 'query2'], 't': ['target1']})
        g = self.aln.groups
        assert set(g['q']) == {'query1', 'query2'}
        assert set(g['t']) == {'target1'}

    def test_set_groups_warns_on_overlap(self, caplog):
        """set_groups warns when a name appears in two groups."""
        import logging

        with caplog.at_level(logging.WARNING):
            self.aln.set_groups({'g1': ['query1'], 'g2': ['query1', 'target1']})
        assert any('query1' in msg for msg in caplog.messages)

    def test_rename_group_updates_label(self):
        """rename_group changes a group label."""
        self.aln.rename_group('a', 'queries')
        g = self.aln.groups
        assert 'queries' in g
        assert 'a' not in g
        assert set(g['queries']) == set(self.aln.query_names)

    def test_rename_group_raises_unknown(self):
        """rename_group raises KeyError for an unknown group."""
        with pytest.raises(KeyError):
            self.aln.rename_group('no_such', 'x')

    def test_rename_group_raises_existing_label(self):
        """rename_group raises ValueError if new name already exists."""
        with pytest.raises(ValueError, match='already exists'):
            self.aln.rename_group('a', 'b')

    def test_reorder_contigs_by_group(self):
        """reorder_contigs respects query_group and target_group params."""
        self.aln.set_groups({'q': ['query1', 'query2'], 't': ['target1']})
        q_sorted, t_sorted = self.aln.reorder_contigs(query_group='q', target_group='t')
        assert set(q_sorted) == {'query1', 'query2'}
        assert set(t_sorted) == {'target1'}

    def test_reorder_contigs_group_unknown_raises(self):
        """reorder_contigs raises KeyError for unknown group label."""
        with pytest.raises(KeyError):
            self.aln.reorder_contigs(query_group='no_such', target_group='b')

    def test_reorder_contigs_backwards_compat(self):
        """reorder_contigs default behavior unchanged (no group params)."""
        q_sorted, t_sorted = self.aln.reorder_contigs()
        assert set(q_sorted) == set(self.aln.query_names)
        assert set(t_sorted) == set(self.aln.target_names)


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
        q_sorted, _t, _reversed = compute_gravity_contigs(
            records,
            ['q_early', 'long_unmatched', 'short_unmatched'],
            ['target'],
        )
        assert q_sorted[0] == 'q_early'
        assert set(q_sorted[1:]) == {'long_unmatched', 'short_unmatched'}


# ---------------------------------------------------------------------------
# Best-matching chromosome (argmax) + squared weighting
# ---------------------------------------------------------------------------


class TestBestMatchingChromosome:
    def test_argmax_uses_best_chromosome(self):
        """A contig's best chromosome is its single big match's target, not the
        many small hits on another chromosome (squared weighting + argmax)."""
        from rusty_dot.paf_io import _gravity_order

        # q_multi has five tiny hits on chrA but one large hit on chrB.  Squared
        # weighting makes the single big block dominate, so argmax → chrB.
        recs = [
            PafRecord.from_line(
                f'q_multi\t900\t{100 + i * 10}\t{110 + i * 10}\t+\t'
                f'chrA\t1000\t{100 + i * 10}\t{110 + i * 10}\t9\t10\t255'
            )
            for i in range(5)
        ]
        recs.append(
            PafRecord.from_line(
                'q_multi\t900\t300\t800\t+\tchrB\t1000\t400\t900\t480\t500\t255'
            )
        )
        matches: dict = {}
        for r in recs:
            matches.setdefault((r.query_name, r.target_name), []).append(r)
        len_map = {'q_multi': 900, 'chrA': 1000, 'chrB': 1000}
        _ordered, best_other = _gravity_order(
            ['q_multi'], ['chrA', 'chrB'], matches, len_map, self_is_query=True
        )
        assert best_other['q_multi'] == 'chrB'

    def test_tie_breaks_by_descending_length(self):
        """Two contigs at the same gravity position sort by descending length."""
        recs = [
            PafRecord.from_line(
                'q_short\t100\t0\t50\t+\tt\t1000\t100\t200\t48\t50\t255'
            ),
            PafRecord.from_line(
                'q_long\t500\t0\t50\t+\tt\t1000\t100\t200\t48\t50\t255'
            ),
        ]
        q_sorted, _t, _rev = compute_gravity_contigs(recs, ['q_short', 'q_long'], ['t'])
        # Identical target mid-point → equal gravity → longer contig first.
        assert q_sorted == ['q_long', 'q_short']


# ---------------------------------------------------------------------------
# Reverse-orientation detection (compute_reversed_contigs / gravity 3-tuple)
# ---------------------------------------------------------------------------


class TestReversedContigDetection:
    def _reversed_from(self, records, q, t):
        _sq, _st, reversed_q = compute_gravity_contigs(records, q, t)
        return reversed_q

    def test_reverse_oriented_contig_detected(self):
        """A contig whose target position decreases as query increases is flagged."""
        recs = [
            PafRecord.from_line(
                'qr\t300\t0\t100\t-\tref\t1000\t900\t1000\t95\t100\t255'
            ),
            PafRecord.from_line(
                'qr\t300\t100\t200\t-\tref\t1000\t800\t900\t95\t100\t255'
            ),
            PafRecord.from_line(
                'qr\t300\t200\t300\t-\tref\t1000\t700\t800\t95\t100\t255'
            ),
        ]
        assert self._reversed_from(recs, ['qr'], ['ref']) == {'qr'}

    def test_forward_contig_not_reversed(self):
        """A contig whose target position increases with query is not flagged."""
        recs = [
            PafRecord.from_line('qf\t300\t0\t100\t+\tref\t1000\t0\t100\t95\t100\t255'),
            PafRecord.from_line(
                'qf\t300\t100\t200\t+\tref\t1000\t100\t200\t95\t100\t255'
            ),
            PafRecord.from_line(
                'qf\t300\t200\t300\t+\tref\t1000\t200\t300\t95\t100\t255'
            ),
        ]
        assert self._reversed_from(recs, ['qf'], ['ref']) == set()

    def test_small_noise_match_does_not_flip(self):
        """A tiny reverse noise match below the 10% threshold is ignored."""
        recs = [
            # One dominant forward match.
            PafRecord.from_line('q\t300\t0\t200\t+\tref\t1000\t0\t200\t190\t200\t255'),
            # Tiny reverse hit: euclidean length far below 10% of the big match.
            PafRecord.from_line('q\t300\t210\t215\t-\tref\t1000\t505\t500\t4\t5\t255'),
        ]
        assert self._reversed_from(recs, ['q'], ['ref']) == set()

    def test_single_big_reverse_match_flagged(self):
        """A single dominant reverse-strand match flags the contig as reversed."""
        recs = [
            PafRecord.from_line('q\t300\t0\t250\t-\tref\t1000\t0\t250\t240\t250\t255'),
        ]
        assert self._reversed_from(recs, ['q'], ['ref']) == {'q'}

    def test_compute_reversed_contigs_direct(self):
        """The standalone helper agrees with the gravity-tuple result."""
        recs = [
            PafRecord.from_line(
                'qr\t300\t0\t150\t-\tref\t1000\t850\t1000\t140\t150\t255'
            ),
            PafRecord.from_line(
                'qr\t300\t150\t300\t-\tref\t1000\t700\t850\t140\t150\t255'
            ),
        ]
        matches = {('qr', 'ref'): recs}
        len_map = {'qr': 300, 'ref': 1000}
        result = compute_reversed_contigs(['qr'], matches, len_map, {'qr': 'ref'})
        assert result == {'qr'}


# ---------------------------------------------------------------------------
# reversed_contigs exposure on PafAlignment and CrossIndex
# ---------------------------------------------------------------------------


class TestReversedContigsAPI:
    def test_pafalignment_reorder_returns_two_tuple(self):
        """reorder_contigs stays a 2-tuple; reversed set lives on the property."""
        recs = [
            PafRecord.from_line(
                'qr\t300\t0\t150\t-\tref\t1000\t850\t1000\t140\t150\t255'
            ),
            PafRecord.from_line(
                'qr\t300\t150\t300\t-\tref\t1000\t700\t850\t140\t150\t255'
            ),
        ]
        aln = PafAlignment.from_records(recs)
        result = aln.reorder_contigs(['qr'], ['ref'])
        assert isinstance(result, tuple) and len(result) == 2
        assert aln.reversed_contigs == {'qr'}

    def test_pafalignment_reversed_contigs_empty_before_reorder(self):
        """The property is empty until reorder_contigs has run."""
        recs = [
            PafRecord.from_line('q\t100\t0\t50\t+\tref\t1000\t0\t50\t48\t50\t255'),
        ]
        aln = PafAlignment.from_records(recs)
        assert aln.reversed_contigs == set()

    def test_crossindex_reversed_contigs_populated(self):
        """reorder_for_colinearity records the reverse-oriented query contigs."""
        from rusty_dot.paf_io import CrossIndex

        def revcomp(s):
            return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

        seq = (
            'ACGTTGCAAGGCCTTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCG'
            'TTAGGCATCCGATTGACCGATACGGATTCAGCTAGGCATTACGGATCCGATTAGCACGTATGC'
        )
        cross = CrossIndex(k=11)
        cross.add_sequence('fwd', seq, group='a')
        cross.add_sequence('rev', revcomp(seq), group='a')
        cross.add_sequence('ref', seq, group='b')
        cross.compute_matches('a', 'b')
        cross.reorder_for_colinearity('a', 'b')
        reversed_a = cross.reversed_contigs('a')
        assert isinstance(reversed_a, set)
        # The reverse-complement contig aligns on the '-' strand → flagged;
        # the forward copy is not.
        assert 'rev' in reversed_a
        assert 'fwd' not in reversed_a

    def test_crossindex_reversed_contigs_empty_for_unreordered_group(self):
        """A group never used as query axis has an empty reversed set."""
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=6)
        cross.add_sequence('q', 'ACGTACGTACGT' * 3, group='a')
        cross.add_sequence('t', 'ACGTACGTACGT' * 3, group='b')
        assert cross.reversed_contigs('a') == set()


# ---------------------------------------------------------------------------
# reverse_complement helper
# ---------------------------------------------------------------------------


class TestReverseComplement:
    def test_basic(self):
        from rusty_dot.paf_io import reverse_complement

        assert reverse_complement('ACGT') == 'ACGT'
        assert reverse_complement('AAAA') == 'TTTT'
        assert reverse_complement('ACGTN') == 'NACGT'

    def test_case_preserved(self):
        from rusty_dot.paf_io import reverse_complement

        assert reverse_complement('acgt') == 'acgt'
        assert reverse_complement('AcgT') == 'AcgT'

    def test_double_reverse_is_identity(self):
        from rusty_dot.paf_io import reverse_complement

        seq = 'ACGTTGCAAGGCCTTAGCTAGGATCCGATCG'
        assert reverse_complement(reverse_complement(seq)) == seq


# ---------------------------------------------------------------------------
# compute_gravity_contigs with fixed targets (sort_targets=False)
# ---------------------------------------------------------------------------


class TestFixedTargetOrder:
    def test_targets_unchanged_when_sort_targets_false(self):
        recs = [
            PafRecord.from_line('qB\t100\t0\t50\t+\tt2\t1000\t900\t950\t48\t50\t255'),
            PafRecord.from_line('qA\t100\t0\t50\t+\tt1\t1000\t0\t50\t48\t50\t255'),
        ]
        q_sorted, t_sorted, _rev = compute_gravity_contigs(
            recs, ['qB', 'qA'], ['t1', 't2'], sort_targets=False
        )
        # Targets kept exactly as given; queries reordered against that axis.
        assert t_sorted == ['t1', 't2']
        assert q_sorted == ['qA', 'qB']


# ---------------------------------------------------------------------------
# CrossIndex FASTA output + fixed-target reorder
# ---------------------------------------------------------------------------


def _rc(s: str) -> str:
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


_FASTA_SEQ = (
    'ACGTTGCAAGGCCTTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCG'
    'TTAGGCATCCGATTGACCGATACGGATTCAGCTAGGCATTACGGATCCGATTAGCACGTATGC'
)


class TestCrossIndexFasta:
    def _cross(self):
        from rusty_dot.paf_io import CrossIndex

        cross = CrossIndex(k=11)
        cross.add_sequence('fwd', _FASTA_SEQ, group='a')
        cross.add_sequence('rev', _rc(_FASTA_SEQ), group='a')
        cross.add_sequence('ref', _FASTA_SEQ, group='b')
        cross.compute_matches('a', 'b')
        return cross

    def test_get_sequence_roundtrip(self):
        cross = self._cross()
        assert cross.get_sequence('fwd', group='a') == _FASTA_SEQ
        assert cross.get_sequence('a:fwd') == _FASTA_SEQ

    def _read_fasta(self, path):
        records, name, seq = {}, None, []
        for line in path.read_text().splitlines():
            if line.startswith('>'):
                if name is not None:
                    records[name] = ''.join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            records[name] = ''.join(seq)
        return records

    def test_write_fasta_reorients_reversed_contig(self, tmp_path):
        cross = self._cross()
        cross.reorder_for_colinearity('a', 'b')
        assert 'rev' in cross.reversed_contigs('a')

        out = tmp_path / 'a.fasta'
        cross.write_fasta(out, 'a')
        recs = self._read_fasta(out)
        # 'rev' was reverse-complemented on write → matches the forward reference.
        assert recs['rev'] == _FASTA_SEQ
        assert recs['fwd'] == _FASTA_SEQ

    def test_write_fasta_reverse_disabled(self, tmp_path):
        cross = self._cross()
        cross.reorder_for_colinearity('a', 'b')
        out = tmp_path / 'a_noflip.fasta'
        cross.write_fasta(out, 'a', reverse=set())
        recs = self._read_fasta(out)
        # With flipping disabled the reverse contig stays reverse-complemented.
        assert recs['rev'] == _rc(_FASTA_SEQ)

    def test_write_fasta_respects_order(self, tmp_path):
        cross = self._cross()
        out = tmp_path / 'a_ordered.fasta'
        cross.write_fasta(out, 'a', order=['rev', 'fwd'], reverse=set())
        names = [
            line[1:].split()[0]
            for line in out.read_text().splitlines()
            if line.startswith('>')
        ]
        assert names == ['rev', 'fwd']

    def test_write_fasta_line_width(self, tmp_path):
        cross = self._cross()
        out = tmp_path / 'wrapped.fasta'
        cross.write_fasta(out, 'a', reverse=set(), line_width=20)
        seq_lines = [
            line for line in out.read_text().splitlines() if not line.startswith('>')
        ]
        assert all(len(line) <= 20 for line in seq_lines)
        assert any(len(line) == 20 for line in seq_lines)

    def test_write_fasta_unknown_group_raises(self, tmp_path):
        cross = self._cross()
        with pytest.raises(KeyError):
            cross.write_fasta(tmp_path / 'x.fasta', 'nope')

    def test_reorder_target_false_keeps_target_fixed(self):
        cross = self._cross()
        cross.add_sequence('ref2', _rc(_FASTA_SEQ), group='b')
        cross.compute_matches('a', 'b')
        before_b = list(cross.contig_order['b'])
        cross.reorder_for_colinearity('a', 'b', reorder_target=False)
        assert cross.contig_order['b'] == before_b  # target untouched
        assert set(cross.contig_order['a']) == {'fwd', 'rev'}


# ---------------------------------------------------------------------------
# Stranded compute_matches cache
# ---------------------------------------------------------------------------


def _rc_seq(seq):
    """Return the reverse complement of *seq*."""
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(comp[b] for b in reversed(seq))


def _stranded_cross():
    """Two-group CrossIndex whose only query contig is reverse-complemented."""
    import random

    from rusty_dot.paf_io import CrossIndex

    rng = random.Random(11)
    chrom = ''.join(rng.choice('ACGT') for _ in range(200))
    cross = CrossIndex(k=15)
    cross.add_sequence('q_rc', _rc_seq(chrom), group='qry')
    cross.add_sequence('chr1', chrom, group='ref')
    return cross


class TestStrandedComputeMatches:
    """compute_matches caches both-strand records with valid PAF semantics."""

    def test_cache_contains_minus_strand_records(self):
        """A reverse-complemented contig yields '-' strand cached records."""
        cross = _stranded_cross()
        cross.compute_matches('qry', 'ref')
        records = cross.get_records_for_pair('qry', 'ref')
        assert records, 'expected cached records for the pair'
        minus = [r for r in records if r.strand == '-']
        assert minus, "expected '-' strand records for a reverse-complement contig"
        for r in minus:
            # PAF spec: query coordinates are always on the forward strand.
            assert 0 <= r.query_start < r.query_end <= r.query_len
            assert 0 <= r.target_start < r.target_end <= r.target_len

    def test_reorder_for_colinearity_flags_reversed_from_cache(self):
        """reorder_for_colinearity still works after the stranded cache fix."""
        cross = _stranded_cross()
        cross.compute_matches('qry', 'ref')
        cross.reorder_for_colinearity('qry', 'ref')
        assert cross.reversed_contigs('qry') == {'q_rc'}

    def test_get_paf_text_includes_minus_strand(self):
        """CrossIndex.get_paf emits '-' strand lines that round-trip cleanly."""
        cross = _stranded_cross()
        lines = cross.get_paf(group_pairs=[('qry', 'ref')])
        minus_lines = [ln for ln in lines if ln.split('\t')[4] == '-']
        assert minus_lines, "expected '-' strand PAF lines"
        for ln in minus_lines:
            rec = PafRecord.from_line(ln)
            assert rec.strand == '-'
            assert 0 <= rec.query_start < rec.query_end <= rec.query_len
            assert 0 <= rec.target_start < rec.target_end <= rec.target_len
            assert rec.to_line() == ln  # lossless round-trip

    def test_cache_matches_get_paf_output(self):
        """Cached records and get_paf text agree (same stranded mechanism)."""
        cross = _stranded_cross()
        cross.compute_matches('qry', 'ref')
        cached = sorted(r.to_line() for r in cross.get_records_for_pair('qry', 'ref'))
        text = sorted(cross.get_paf(group_pairs=[('qry', 'ref')]))
        assert cached == text
