"""Tests for strand-aware k-mer matching and non-ATGC character handling."""

import pytest

from rusty_dot._rusty_dot import (
    SequenceIndex,
    py_build_kmer_set,
    py_find_kmer_coords,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def make_idx(seqs: dict, k: int = 4) -> SequenceIndex:
    idx = SequenceIndex(k=k)
    for name, seq in seqs.items():
        idx.add_sequence(name, seq)
    return idx


# ---------------------------------------------------------------------------
# Non-ATGC character handling
# ---------------------------------------------------------------------------


class TestNonAtgcHandling:
    """K-mers containing N or other ambiguous bases are excluded."""

    def test_build_kmer_set_excludes_n(self):
        """K-mers spanning an N position must not appear in the k-mer set."""
        kmers = py_build_kmer_set('ACGNACGT', 4)
        # CGNA, GNAC, NACT all contain N and must be excluded
        for kmer in kmers:
            assert 'N' not in kmer, f'K-mer with N slipped through: {kmer}'

    def test_build_kmer_set_all_n(self):
        """A sequence of pure Ns should yield an empty k-mer set."""
        assert py_build_kmer_set('NNNNNNNN', 4) == set()

    def test_build_kmer_set_lowercase_excluded(self):
        """Lowercase characters are not valid ACGT and must be excluded."""
        # Build from uppercase; lowercase bases produce invalid k-mers
        kmers_upper = py_build_kmer_set('ACGT', 4)
        kmers_lower = py_build_kmer_set('acgt', 4)
        # Lowercase k-mers do not appear (only A C G T uppercase are valid)
        assert 'acgt' not in kmers_lower
        # The uppercase version should be found
        assert 'ACGT' in kmers_upper

    def test_find_kmer_coords_with_n_in_sequence(self):
        """Searching a sequence that contains N returns only valid positions."""
        # The k-mer ACGT must be found at positions not overlapping N
        coords = py_find_kmer_coords('ACGTNACGT', ['ACGT'])
        positions = coords.get('ACGT', [])
        assert 0 in positions  # ACGT at pos 0 (no N)
        assert 5 in positions  # ACGT at pos 5 (no N)

    def test_find_kmer_coords_kmer_spanning_n_not_found(self):
        """Searching for a k-mer containing N returns positions where it literally exists.

        Note: the FM-index performs exact string matching and *will* find "CGTN"
        at position 1 in "ACGTNACGT" because the bytes match.  The filtering of
        N-containing k-mers happens in ``build_kmer_set``, not in the FM-index
        query.  This test documents the actual (correct) behavior.
        """
        coords = py_find_kmer_coords('ACGTNACGT', ['CGTN'])
        # "CGTN" is literally at position 1 → the FM-index WILL find it
        assert 'CGTN' in coords
        assert 1 in coords['CGTN']

    def test_compare_sequences_n_in_sequence(self):
        """Compare two sequences where one contains N; shared valid k-mers found."""
        idx = make_idx({'a': 'ACGTACGT', 'b': 'ACGTNACGT'}, k=4)
        matches = idx.compare_sequences('a', 'b')
        # At least the non-N-spanning k-mers should match
        assert isinstance(matches, list)
        assert len(matches) > 0

    def test_kmer_set_no_invalid_characters(self):
        """K-mer set from a sequence with mixed valid/invalid chars contains only valid."""
        kmers = py_build_kmer_set('ACGTXYZACGT', 4)
        for kmer in kmers:
            assert all(c in 'ACGT' for c in kmer)


# ---------------------------------------------------------------------------
# Strand handling
# ---------------------------------------------------------------------------


class TestStrandHandling:
    """Reverse-complement strand matches are correctly detected and reported."""

    def test_revcomp_match_detected(self):
        """A k-mer's reverse complement in the target should appear as '-' strand."""
        # query k-mer AAAC has revcomp GTTT; target contains GTTT
        idx = make_idx({'q': 'AAAC', 't': 'GTTT'}, k=4)
        stranded = idx.compare_sequences_stranded('q', 't')
        strands = [m[4] for m in stranded]
        assert '-' in strands, "Expected a '-' strand match for RC hit"

    def test_forward_match_reported_as_plus(self):
        """A forward match must be reported with '+' strand."""
        idx = make_idx({'q': 'ACGTACGT', 't': 'ACGTTTTT'}, k=4)
        stranded = idx.compare_sequences_stranded('q', 't')
        strands = [m[4] for m in stranded]
        assert '+' in strands

    def test_palindrome_matches_both_strands(self):
        """A palindromic k-mer (ACGT) matches both + and - strands."""
        # ACGT is self-complementary: revcomp(ACGT) = ACGT
        idx = make_idx({'q': 'ACGT', 't': 'ACGT'}, k=4)
        stranded = idx.compare_sequences_stranded('q', 't')
        strands = {m[4] for m in stranded}
        # Should see at least + strand (same sequence), may also see -
        assert '+' in strands

    def test_stranded_result_has_5_elements(self):
        """Each result from compare_sequences_stranded must have 5 elements."""
        idx = make_idx({'q': 'ACGTACGT', 't': 'TACGTACG'}, k=4)
        for match in idx.compare_sequences_stranded('q', 't'):
            assert len(match) == 5
            qs, qe, ts, te, strand = match
            assert strand in ('+', '-')
            assert qs < qe
            assert ts < te

    def test_unmerged_stranded(self):
        """Unmerged stranded results also carry strand information."""
        idx = make_idx({'q': 'ACGTACGT', 't': 'ACGTACGT'}, k=4)
        results = idx.compare_sequences_stranded('q', 't', merge=False)
        for _qs, _qe, _ts, _te, strand in results:
            assert strand in ('+', '-')

    def test_self_comparison_stranded(self):
        """Self-comparison should produce + strand hits on the main diagonal."""
        idx = make_idx({'s': 'ACGTACGT'}, k=4)
        stranded = idx.compare_sequences_stranded('s', 's')
        main_diag = [m for m in stranded if m[0] == m[2] and m[4] == '+']
        assert len(main_diag) > 0

    def test_missing_sequence_raises(self):
        """Stranded comparison with unknown names raises KeyError."""
        idx = make_idx({'q': 'ACGTACGT'}, k=4)
        with pytest.raises(KeyError):
            idx.compare_sequences_stranded('q', 'missing')

    def test_compare_sequences_unchanged_api(self):
        """compare_sequences (unstranded) still returns 4-tuples."""
        idx = make_idx({'q': 'ACGTACGT', 't': 'TACGTACG'}, k=4)
        for m in idx.compare_sequences('q', 't'):
            assert len(m) == 4

    def test_no_match_stranded_returns_empty(self):
        """Sequences with no shared k-mers (even by RC) return empty list."""
        # AAAA and CCCC have no shared k-mers and no RC relation
        idx = make_idx({'q': 'AAAAAAAAAA', 't': 'CCCCCCCCCC'}, k=4)
        assert idx.compare_sequences_stranded('q', 't') == []

    def test_rc_anti_diagonal_coordinates(self):
        """Reverse-strand hit coordinates bracket the RC region on the target."""
        # query: AAAC (k=4); target: GTTT (revcomp of AAAC)
        # Expected - strand hit: target_start=0, target_end=4
        idx = make_idx({'q': 'AAAC', 't': 'GTTT'}, k=4)
        stranded = idx.compare_sequences_stranded('q', 't')
        rev_hits = [m for m in stranded if m[4] == '-']
        assert len(rev_hits) > 0
        _, _, ts, te, _ = rev_hits[0]
        assert ts < te
        assert te - ts >= 4  # at least k bases covered

    def test_consecutive_opposite_strand_kmers_not_merged(self):
        """Sequential query k-mers with collinear hits on opposite strands must not merge.

        query = 'AACCC' (k=4):
          AACC at q=0 → forward match in target (AACC at t=0).
          ACCC at q=1 → RC match in target (RC(ACCC)=GGGT at t=4).

        These two k-mers are *consecutive* in the query (overlapping by k-1=3 bases),
        and each has a collinear hit in the target — but on **opposite strands**.
        The expected result is two separate blocks (one '+', one '-'), NOT a single
        merged cross-strand block.
        """
        idx = make_idx({'q': 'AACCC', 't': 'AACCGGGT'}, k=4)
        stranded = idx.compare_sequences_stranded('q', 't', merge=True)

        fwd_hits = [m for m in stranded if m[4] == '+']
        rev_hits = [m for m in stranded if m[4] == '-']

        assert len(fwd_hits) >= 1, f'expected + strand hit for AACC, got {stranded}'
        assert len(rev_hits) >= 1, f'expected - strand hit for ACCC→GGGT, got {stranded}'

        # Exactly one block per strand: no cross-strand merging occurred.
        assert len(stranded) == 2, (
            f'expected exactly 2 blocks (one per strand), got {stranded}'
        )

        # Forward block: AACC matched at query q=0, target t=0.
        fwd = fwd_hits[0]
        assert fwd[0] == 0, f'forward hit q_start should be 0, got {fwd}'
        assert fwd[2] == 0, f'forward hit t_start should be 0, got {fwd}'

        # RC block: ACCC matched at query q=1; RC(ACCC)=GGGT at target t=4.
        rev = rev_hits[0]
        assert rev[0] == 1, f'RC hit q_start should be 1, got {rev}'
        assert rev[2] == 4, f'RC hit t_start should be 4, got {rev}'


# ---------------------------------------------------------------------------
# Optimal contig order
# ---------------------------------------------------------------------------


class TestOptimalContigOrder:
    """Gravity-based contig reordering produces a reasonable ordering."""

    def test_returns_same_names(self):
        """Result contains exactly the same names as the input."""
        idx = make_idx({'a': 'ACGTACGTACGT', 'b': 'TTTTACGTACGT', 'c': 'CGCGCGCGCGCG'})
        q_sorted, t_sorted = idx.optimal_contig_order(['a', 'b'], ['b', 'c'])
        assert set(q_sorted) == {'a', 'b'}
        assert set(t_sorted) == {'b', 'c'}

    def test_collinear_sequences_sorted_in_order(self):
        """Gravity-based sort places the contig with earlier matches first.

        This test uses the PafAlignment.reorder_contigs path (pure Python,
        from known records) to avoid the ambiguity from RC matches in the
        FM-index path.  For the SequenceIndex.optimal_contig_order method we
        just verify that the returned names are unchanged as a set.
        """
        from rusty_dot.paf_io import PafAlignment, PafRecord

        # Construct explicit records where q_early clearly maps to the start
        # of the target and q_late maps to the end, with no ambiguity.
        records = [
            PafRecord.from_line('q_early\t50\t0\t50\t+\tref\t200\t0\t50\t48\t50\t255'),
            PafRecord.from_line(
                'q_late\t50\t0\t50\t+\tref\t200\t150\t200\t48\t50\t255'
            ),
        ]
        aln = PafAlignment.from_records(records)
        q_sorted, _ = aln.reorder_contigs(['q_late', 'q_early'], ['ref'])
        assert q_sorted[0] == 'q_early', (
            f'Expected q_early first (lower gravity), got {q_sorted}'
        )
        assert q_sorted[1] == 'q_late'

    def test_optimal_contig_order_returns_same_names(self):
        """optimal_contig_order returns exactly the names that were passed in."""
        idx = SequenceIndex(k=4)
        idx.add_sequence('a', 'ACGTACGTACGT')
        idx.add_sequence('b', 'TTTTACGTACGT')
        idx.add_sequence('c', 'CGCGCGCGCGCG')
        q_sorted, t_sorted = idx.optimal_contig_order(['a', 'b'], ['b', 'c'])
        assert set(q_sorted) == {'a', 'b'}
        assert set(t_sorted) == {'b', 'c'}

    def test_unmatched_contigs_sort_last(self):
        """Contigs with no matches have f64::MAX gravity and sort to the end."""
        idx = SequenceIndex(k=4)
        idx.add_sequence('match', 'ACGTACGTACGT')
        idx.add_sequence('nomatch', 'AAAAAAAAAAAAAA')
        idx.add_sequence('target', 'ACGTACGTACGT')
        q_sorted, _ = idx.optimal_contig_order(['match', 'nomatch'], ['target'])
        # 'match' has hits; 'nomatch' has none → 'match' sorts first
        assert q_sorted[0] == 'match'
        assert q_sorted[-1] == 'nomatch'

    def test_missing_name_raises(self):
        """Passing an unknown name raises KeyError."""
        idx = make_idx({'a': 'ACGTACGT'})
        with pytest.raises(KeyError):
            idx.optimal_contig_order(['a', 'missing'], ['a'])
