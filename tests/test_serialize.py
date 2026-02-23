"""Tests for FM-index serialization."""

import os
import pytest

from rusty_dot._rusty_dot import py_save_index, py_load_index


def test_save_and_load_index(tmp_path):
    """Test saving and loading an index collection."""
    path = str(tmp_path / "index.bin")
    sequences = {
        "seq1": "ACGTACGTACGT",
        "seq2": "TACGTACGTACG",
    }
    py_save_index(path, sequences, k=4)
    assert os.path.exists(path)

    kmer_sets, k = py_load_index(path)
    assert k == 4
    assert "seq1" in kmer_sets
    assert "seq2" in kmer_sets


def test_load_nonexistent_file():
    """Test that loading a nonexistent file raises ValueError."""
    with pytest.raises(ValueError):
        py_load_index("/nonexistent/path/index.bin")
