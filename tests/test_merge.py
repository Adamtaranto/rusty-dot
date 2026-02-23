"""Tests for k-mer run merging."""

import pytest

from rusty_dot._rusty_dot import py_merge_kmer_runs


def test_merge_consecutive_runs():
    """Test that consecutive k-mer hits are merged into one block."""
    # k-mer "ACGT" at query pos 0, target pos 10
    # k-mer "CGTA" at query pos 1, target pos 11 (consecutive)
    k = 4
    kmer_coords = {"ACGT": [10], "CGTA": [11]}
    query_positions = {"ACGT": [0], "CGTA": [1]}
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    # Should merge into a single run
    assert len(merged) == 1
    qs, qe, ts, te = merged[0]
    assert qs == 0
    assert qe == 1 + k  # = 5
    assert ts == 10
    assert te == 11 + k  # = 15


def test_merge_non_consecutive_stays_separate():
    """Test that non-consecutive hits remain separate."""
    k = 4
    kmer_coords = {"ACGT": [0], "TTTT": [20]}
    query_positions = {"ACGT": [0], "TTTT": [5]}
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    assert len(merged) == 2


def test_merge_empty_input():
    """Test that empty input returns empty list."""
    merged = py_merge_kmer_runs({}, {}, 4)
    assert merged == []


def test_merge_multiple_target_positions():
    """Test handling of k-mers with multiple target positions."""
    k = 4
    # ACGT appears at target positions 0 and 100
    kmer_coords = {"ACGT": [0, 100]}
    query_positions = {"ACGT": [0]}
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    # Two separate matches (same query pos, different target pos)
    assert len(merged) == 2
