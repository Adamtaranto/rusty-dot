"""Tests for FASTA reading functionality."""

import pytest

from rusty_dot._rusty_dot import py_read_fasta


def test_read_fasta_plain(fasta_file):
    """Test reading a plain FASTA file."""
    seqs = py_read_fasta(fasta_file)
    assert 'seq1' in seqs
    assert 'seq2' in seqs
    assert 'seq3' in seqs
    assert seqs['seq1'] == 'ACGTACGTACGTACGTACGT'
    assert seqs['seq2'] == 'TACGTACGTACGTACGTACG'


def test_read_fasta_gzip(gzip_fasta_file):
    """Test reading a gzipped FASTA file."""
    seqs = py_read_fasta(gzip_fasta_file)
    assert 'seq1' in seqs
    assert 'seq2' in seqs
    assert seqs['seq1'] == 'ACGTACGTACGTACGTACGT'


def test_read_fasta_missing_file():
    """Test that a missing file raises ValueError."""
    with pytest.raises(ValueError):
        py_read_fasta('/nonexistent/path/file.fasta')


def test_read_fasta_returns_uppercase(fasta_file):
    """Test that sequences are returned in uppercase."""
    seqs = py_read_fasta(fasta_file)
    for seq in seqs.values():
        assert seq == seq.upper()
