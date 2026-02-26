"""Pytest configuration and shared fixtures."""

import gzip
import os

import pytest

# Ensure a non-interactive backend is used for tests running in headless
# environments (e.g. CI).  This must be set before pyplot is imported.
os.environ.setdefault('MPLBACKEND', 'Agg')

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
