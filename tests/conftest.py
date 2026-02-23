"""Pytest configuration and shared fixtures."""

import gzip

import pytest

from tests.test_data import FASTA_CONTENT


@pytest.fixture
def fasta_file(tmp_path):
    """Write a plain FASTA file and return its path."""
    path = tmp_path / 'test.fasta'
    path.write_text(FASTA_CONTENT)
    return str(path)


@pytest.fixture
def gzip_fasta_file(tmp_path):
    """Write a gzipped FASTA file and return its path."""
    path = tmp_path / 'test.fasta.gz'
    with gzip.open(str(path), 'wt') as f:
        f.write(FASTA_CONTENT)
    return str(path)
