"""Tests for the SequenceIndex class."""

import os

import pytest

from rusty_dot._rusty_dot import SequenceIndex


@pytest.fixture
def simple_index():
    """Create a simple index with two sequences."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGTACGTACGTACGT')
    idx.add_sequence('seq2', 'TACGTACGTACGTACG')
    return idx


def test_index_creation():
    """Test basic index creation."""
    idx = SequenceIndex(k=5)
    assert idx.k == 5
    assert len(idx) == 0


def test_index_invalid_k():
    """Test that k=0 raises ValueError."""
    with pytest.raises(ValueError):
        SequenceIndex(k=0)


def test_add_sequence(simple_index):
    """Test that sequences are added correctly."""
    assert len(simple_index) == 2
    names = simple_index.sequence_names()
    assert 'seq1' in names
    assert 'seq2' in names


def test_get_sequence_length(simple_index):
    """Test sequence length retrieval."""
    assert simple_index.get_sequence_length('seq1') == 16


def test_get_sequence_length_missing(simple_index):
    """Test that missing sequence raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.get_sequence_length('nonexistent')


def test_get_kmer_set(simple_index):
    """Test k-mer set retrieval."""
    kmers = simple_index.get_kmer_set('seq1')
    assert isinstance(kmers, set)
    assert 'ACGT' in kmers
    assert len(kmers) > 0


def test_get_kmer_set_missing(simple_index):
    """Test that missing sequence raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.get_kmer_set('nonexistent')


def test_compare_sequences_returns_coords(simple_index):
    """Test that compare_sequences returns coordinate tuples."""
    matches = simple_index.compare_sequences('seq1', 'seq2')
    assert isinstance(matches, list)
    for match in matches:
        assert len(match) == 4
        qs, qe, ts, te = match
        assert qs < qe
        assert ts < te


def test_compare_sequences_no_overlap():
    """Test comparison with no shared k-mers."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('all_a', 'AAAAAAAAAAAAAAAA')
    idx.add_sequence('all_c', 'CCCCCCCCCCCCCCCC')
    matches = idx.compare_sequences('all_a', 'all_c')
    assert matches == []


def test_compare_sequences_caching(simple_index):
    """Test that results are cached (same result on second call)."""
    matches1 = simple_index.compare_sequences('seq1', 'seq2')
    matches2 = simple_index.compare_sequences('seq1', 'seq2')
    assert matches1 == matches2


def test_compare_sequences_unmerged(simple_index):
    """Test unmerged comparison returns individual k-mer hits."""
    merged = simple_index.compare_sequences('seq1', 'seq2', merge=True)
    # Unmerged should have >= merged count
    unmerged_results = simple_index.compare_sequences('seq1', 'seq2', merge=False)
    assert len(unmerged_results) >= len(merged)


def test_get_paf(simple_index):
    """Test PAF output format."""
    paf_lines = simple_index.get_paf('seq1', 'seq2')
    for line in paf_lines:
        fields = line.split('\t')
        assert len(fields) == 12
        assert fields[0] == 'seq1'
        assert fields[5] == 'seq2'


def test_load_fasta(fasta_file):
    """Test loading sequences from FASTA file."""
    idx = SequenceIndex(k=4)
    names = idx.load_fasta(fasta_file)
    assert len(names) == 3
    assert 'seq1' in names


def test_load_fasta_gzip(gzip_fasta_file):
    """Test loading sequences from gzipped FASTA file."""
    idx = SequenceIndex(k=4)
    names = idx.load_fasta(gzip_fasta_file)
    assert len(names) == 3


def test_save_and_load(simple_index, tmp_path):
    """Test saving and loading an index."""
    path = str(tmp_path / 'index.bin')
    simple_index.save(path)
    assert os.path.exists(path)

    # Load into new index with same k
    new_idx = SequenceIndex(k=4)
    new_idx.load(path)
    assert len(new_idx) == 2
    assert set(new_idx.sequence_names()) == set(simple_index.sequence_names())


def test_load_wrong_k(simple_index, tmp_path):
    """Test that loading with wrong k raises ValueError."""
    path = str(tmp_path / 'index.bin')
    simple_index.save(path)

    new_idx = SequenceIndex(k=5)
    with pytest.raises(ValueError):
        new_idx.load(path)


def test_precompute_all_pairs():
    """Test precomputing all pairwise comparisons."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('a', 'ACGTACGTACGT')
    idx.add_sequence('b', 'TACGTACGTACG')
    idx.add_sequence('c', 'CGCGCGCGCGCG')
    pairs = idx.precompute_all_pairs()
    # Should compute n*(n-1) = 3*2 = 6 pairs
    assert len(pairs) == 6


def test_repr(simple_index):
    """Test string representation."""
    r = repr(simple_index)
    assert 'SequenceIndex' in r
    assert 'k=4' in r


def test_compare_sequences_cache_isolated_by_merge_flag(simple_index):
    """Test that merge=True and merge=False cache entries are independent."""
    merged = simple_index.compare_sequences('seq1', 'seq2', merge=True)
    unmerged = simple_index.compare_sequences('seq1', 'seq2', merge=False)
    # Calling again should return the same cached values for each flag
    assert simple_index.compare_sequences('seq1', 'seq2', merge=True) == merged
    assert simple_index.compare_sequences('seq1', 'seq2', merge=False) == unmerged
    # Merging reduces or preserves the count
    assert len(unmerged) >= len(merged)


def test_compare_sequence_self():
    """Test self-comparison produces matches on the main diagonal."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('seq1', 'ACGTACGT')
    matches = idx.compare_sequences('seq1', 'seq1')
    assert len(matches) > 0
    # The main diagonal (q_start == t_start) must be present for an exact self-match
    main_diag = [m for m in matches if m[0] == m[2]]
    assert len(main_diag) > 0


def test_compare_sequences_missing_query(simple_index):
    """Test that a missing query name raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.compare_sequences('no_such_seq', 'seq2')


def test_compare_sequences_missing_target(simple_index):
    """Test that a missing target name raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.compare_sequences('seq1', 'no_such_seq')


def test_get_paf_all_returns_paf_lines():
    """Test that get_paf_all returns PAF lines for all pairs."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('a', 'ACGTACGTACGTACGT')
    idx.add_sequence('b', 'TACGTACGTACGTACG')
    paf_lines = idx.get_paf_all()
    # Each line must have 12 tab-separated fields
    assert isinstance(paf_lines, list)
    for line in paf_lines:
        fields = line.split('\t')
        assert len(fields) == 12


def test_get_paf_all_empty_index_returns_empty():
    """Test that get_paf_all on a single-sequence index returns no lines."""
    idx = SequenceIndex(k=4)
    idx.add_sequence('only', 'ACGTACGTACGT')
    assert idx.get_paf_all() == []


def test_optimal_contig_order_unmatched_sorted_by_length_desc():
    """Unmatched contigs should be placed at the end sorted by descending length."""
    idx = SequenceIndex(k=4)
    # Sequences that share k-mers
    idx.add_sequence('seq_a', 'ACGTACGTACGTACGTACGT')   # 20 bp
    idx.add_sequence('seq_b', 'ACGTACGTACGTACGTACGT')   # 20 bp
    # Two unmatched sequences of different lengths (no ACGT k-mers → no matches)
    idx.add_sequence('long_unmatched', 'T' * 40)         # 40 bp
    idx.add_sequence('short_unmatched', 'T' * 10)        # 10 bp

    q_sorted, _ = idx.optimal_contig_order(
        ['long_unmatched', 'short_unmatched', 'seq_a'],
        ['seq_b'],
    )
    # seq_a matches seq_b → sorted first
    assert q_sorted[0] == 'seq_a'
    # Unmatched sorted by length descending: long before short
    assert q_sorted[1] == 'long_unmatched'
    assert q_sorted[2] == 'short_unmatched'
