"""Tests for the SequenceIndex class."""

import pytest
import tempfile
import os

from rusty_dot._rusty_dot import SequenceIndex


@pytest.fixture
def simple_index():
    """Create a simple index with two sequences."""
    idx = SequenceIndex(k=4)
    idx.add_sequence("seq1", "ACGTACGTACGTACGT")
    idx.add_sequence("seq2", "TACGTACGTACGTACG")
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
    assert "seq1" in names
    assert "seq2" in names


def test_get_sequence_length(simple_index):
    """Test sequence length retrieval."""
    assert simple_index.get_sequence_length("seq1") == 16


def test_get_sequence_length_missing(simple_index):
    """Test that missing sequence raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.get_sequence_length("nonexistent")


def test_get_kmer_set(simple_index):
    """Test k-mer set retrieval."""
    kmers = simple_index.get_kmer_set("seq1")
    assert isinstance(kmers, set)
    assert "ACGT" in kmers
    assert len(kmers) > 0


def test_get_kmer_set_missing(simple_index):
    """Test that missing sequence raises KeyError."""
    with pytest.raises(KeyError):
        simple_index.get_kmer_set("nonexistent")


def test_compare_sequences_returns_coords(simple_index):
    """Test that compare_sequences returns coordinate tuples."""
    matches = simple_index.compare_sequences("seq1", "seq2")
    assert isinstance(matches, list)
    for match in matches:
        assert len(match) == 4
        qs, qe, ts, te = match
        assert qs < qe
        assert ts < te


def test_compare_sequences_no_overlap():
    """Test comparison with no shared k-mers."""
    idx = SequenceIndex(k=4)
    idx.add_sequence("all_a", "AAAAAAAAAAAAAAAA")
    idx.add_sequence("all_c", "CCCCCCCCCCCCCCCC")
    matches = idx.compare_sequences("all_a", "all_c")
    assert matches == []


def test_compare_sequences_caching(simple_index):
    """Test that results are cached (same result on second call)."""
    matches1 = simple_index.compare_sequences("seq1", "seq2")
    matches2 = simple_index.compare_sequences("seq1", "seq2")
    assert matches1 == matches2


def test_compare_sequences_unmerged(simple_index):
    """Test unmerged comparison returns individual k-mer hits."""
    merged = simple_index.compare_sequences("seq1", "seq2", merge=True)
    # Unmerged should have >= merged count
    unmerged_results = simple_index.compare_sequences("seq1", "seq2", merge=False)
    assert len(unmerged_results) >= len(merged)


def test_get_paf(simple_index):
    """Test PAF output format."""
    paf_lines = simple_index.get_paf("seq1", "seq2")
    for line in paf_lines:
        fields = line.split("\t")
        assert len(fields) == 12
        assert fields[0] == "seq1"
        assert fields[5] == "seq2"


def test_load_fasta(fasta_file):
    """Test loading sequences from FASTA file."""
    idx = SequenceIndex(k=4)
    names = idx.load_fasta(fasta_file)
    assert len(names) == 3
    assert "seq1" in names


def test_load_fasta_gzip(gzip_fasta_file):
    """Test loading sequences from gzipped FASTA file."""
    idx = SequenceIndex(k=4)
    names = idx.load_fasta(gzip_fasta_file)
    assert len(names) == 3


def test_save_and_load(simple_index, tmp_path):
    """Test saving and loading an index."""
    path = str(tmp_path / "index.bin")
    simple_index.save(path)
    assert os.path.exists(path)

    # Load into new index with same k
    new_idx = SequenceIndex(k=4)
    new_idx.load(path)
    assert len(new_idx) == 2
    assert set(new_idx.sequence_names()) == set(simple_index.sequence_names())


def test_load_wrong_k(simple_index, tmp_path):
    """Test that loading with wrong k raises ValueError."""
    path = str(tmp_path / "index.bin")
    simple_index.save(path)

    new_idx = SequenceIndex(k=5)
    with pytest.raises(ValueError):
        new_idx.load(path)


def test_precompute_all_pairs():
    """Test precomputing all pairwise comparisons."""
    idx = SequenceIndex(k=4)
    idx.add_sequence("a", "ACGTACGTACGT")
    idx.add_sequence("b", "TACGTACGTACG")
    idx.add_sequence("c", "CGCGCGCGCGCG")
    pairs = idx.precompute_all_pairs()
    # Should compute n*(n-1) = 3*2 = 6 pairs
    assert len(pairs) == 6


def test_repr(simple_index):
    """Test string representation."""
    r = repr(simple_index)
    assert "SequenceIndex" in r
    assert "k=4" in r
