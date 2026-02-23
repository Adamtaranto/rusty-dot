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


def test_merge_two_parallel_diagonals():
    """Test that hits on two parallel diagonals are each merged independently.

    Regression test for the diagonal-sort fix: naive (query, target) sorting
    would interleave hits from diagonal t-q=5 and diagonal t-q=10 and fail to
    merge each pair correctly.

    Diagonal t-q=5:  (q=0,t=5) + (q=1,t=6)  → merged to (qs=0, qe=5, ts=5, te=10)
    Diagonal t-q=10: (q=0,t=10) + (q=1,t=11) → merged to (qs=0, qe=5, ts=10, te=15)
    """
    k = 4
    kmer_coords = {"AAAA": [5, 10], "CCCC": [6, 11]}
    query_positions = {"AAAA": [0], "CCCC": [1]}
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    assert len(merged) == 2, f"Expected 2 merged runs, got {len(merged)}: {merged}"
    merged_set = set(merged)
    # Diagonal t-q=5: q=[0,5), t=[5,10)
    assert (0, 1 + k, 5, 6 + k) in merged_set, f"Missing diagonal-5 run in {merged_set}"
    # Diagonal t-q=10: q=[0,5), t=[10,15)
    assert (0, 1 + k, 10, 11 + k) in merged_set, f"Missing diagonal-10 run in {merged_set}"


def test_merge_longer_run():
    """Test merging a run of five consecutive k-mer hits into one block."""
    k = 3
    # 5 consecutive (q, t) pairs starting at (0, 10)
    kmer_coords = {
        "ACG": [10],
        "CGT": [11],
        "GTA": [12],
        "TAC": [13],
        "ACT": [14],
    }
    query_positions = {
        "ACG": [0],
        "CGT": [1],
        "GTA": [2],
        "TAC": [3],
        "ACT": [4],
    }
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    assert len(merged) == 1
    qs, qe, ts, te = merged[0]
    assert qs == 0
    assert qe == 4 + k  # = 7
    assert ts == 10
    assert te == 14 + k  # = 17


def test_merge_single_kmer():
    """Test that a single k-mer hit produces a single run of length k."""
    k = 5
    merged = py_merge_kmer_runs({"ACGTA": [3]}, {"ACGTA": [7]}, k)
    assert len(merged) == 1
    qs, qe, ts, te = merged[0]
    assert qs == 7
    assert qe == 7 + k
    assert ts == 3
    assert te == 3 + k
