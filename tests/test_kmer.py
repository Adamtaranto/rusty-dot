"""Tests for k-mer set building and coordinate lookup."""

import pytest

from rusty_dot._rusty_dot import py_build_kmer_set, py_find_kmer_coords


def test_build_kmer_set_basic():
    """Test building a k-mer set for a simple sequence."""
    seq = 'ACGTACGT'
    kmers = py_build_kmer_set(seq, 4)
    assert 'ACGT' in kmers
    assert 'CGTA' in kmers
    assert 'GTAC' in kmers
    assert 'TACG' in kmers


def test_build_kmer_set_k_larger_than_seq():
    """Test that k larger than seq length returns empty set."""
    seq = 'ACGT'
    kmers = py_build_kmer_set(seq, 10)
    assert len(kmers) == 0


def test_build_kmer_set_invalid_k():
    """Test that k=0 raises ValueError."""
    with pytest.raises(ValueError):
        py_build_kmer_set('ACGT', 0)


def test_build_kmer_set_excludes_n():
    """Test that k-mers with N are excluded."""
    seq = 'ACGTNACGT'
    kmers = py_build_kmer_set(seq, 4)
    for kmer in kmers:
        assert 'N' not in kmer


def test_build_kmer_set_uniqueness():
    """Test that the k-mer set contains unique entries."""
    seq = 'ACGTACGTACGT'
    kmers = py_build_kmer_set(seq, 4)
    assert len(kmers) == len(set(kmers))


def test_find_kmer_coords_basic():
    """Test finding positions of k-mers in a sequence."""
    seq = 'ACGTACGTACGT'
    kmers = ['ACGT', 'CGTA']
    coords = py_find_kmer_coords(seq, kmers)
    assert 'ACGT' in coords
    # ACGT appears at positions 0, 4, 8 in "ACGTACGTACGT"
    assert 0 in coords['ACGT']
    assert 4 in coords['ACGT']
    assert 8 in coords['ACGT']


def test_find_kmer_coords_no_match():
    """Test that k-mers not in sequence return empty or missing entry."""
    seq = 'ACGTACGT'
    coords = py_find_kmer_coords(seq, ['TTTT'])
    assert 'TTTT' not in coords or coords.get('TTTT', []) == []


def test_find_kmer_coords_empty_kmers():
    """Test with empty k-mer list returns empty dict."""
    coords = py_find_kmer_coords('ACGTACGT', [])
    assert coords == {}
