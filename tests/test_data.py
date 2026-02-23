"""Shared test data and fixtures for rusty-dot tests."""

import gzip
import pytest


# Short test sequences with known shared k-mers
SEQ_A = "ACGTACGTACGTACGTACGT"
SEQ_B = "TACGTACGTACGTACGTACG"
SEQ_C = "GCGCGCGCGCGCGCGCGCGC"

FASTA_CONTENT = """>seq1
ACGTACGTACGTACGTACGT
>seq2
TACGTACGTACGTACGTACG
>seq3
GCGCGCGCGCGCGCGCGCGC
"""


@pytest.fixture
def fasta_file(tmp_path):
    """Write a plain FASTA file and return its path."""
    path = tmp_path / "test.fasta"
    path.write_text(FASTA_CONTENT)
    return str(path)


@pytest.fixture
def gzip_fasta_file(tmp_path):
    """Write a gzipped FASTA file and return its path."""
    path = tmp_path / "test.fasta.gz"
    with gzip.open(str(path), "wt") as f:
        f.write(FASTA_CONTENT)
    return str(path)
