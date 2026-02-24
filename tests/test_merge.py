"""Tests for k-mer run merging."""

from rusty_dot._rusty_dot import py_merge_kmer_runs, py_merge_rev_runs


def test_merge_consecutive_runs():
    """Test that consecutive k-mer hits are merged into one block."""
    # k-mer "ACGT" at query pos 0, target pos 10
    # k-mer "CGTA" at query pos 1, target pos 11 (consecutive)
    k = 4
    kmer_coords = {'ACGT': [10], 'CGTA': [11]}
    query_positions = {'ACGT': [0], 'CGTA': [1]}
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
    kmer_coords = {'ACGT': [0], 'TTTT': [20]}
    query_positions = {'ACGT': [0], 'TTTT': [5]}
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
    kmer_coords = {'ACGT': [0, 100]}
    query_positions = {'ACGT': [0]}
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
    kmer_coords = {'AAAA': [5, 10], 'CCCC': [6, 11]}
    query_positions = {'AAAA': [0], 'CCCC': [1]}
    merged = py_merge_kmer_runs(kmer_coords, query_positions, k)
    assert len(merged) == 2, f'Expected 2 merged runs, got {len(merged)}: {merged}'
    merged_set = set(merged)
    # Diagonal t-q=5: q=[0,5), t=[5,10)
    assert (0, 1 + k, 5, 6 + k) in merged_set, f'Missing diagonal-5 run in {merged_set}'
    # Diagonal t-q=10: q=[0,5), t=[10,15)
    assert (0, 1 + k, 10, 11 + k) in merged_set, (
        f'Missing diagonal-10 run in {merged_set}'
    )


def test_merge_longer_run():
    """Test merging a run of five consecutive k-mer hits into one block."""
    k = 3
    # 5 consecutive (q, t) pairs starting at (0, 10)
    kmer_coords = {
        'ACG': [10],
        'CGT': [11],
        'GTA': [12],
        'TAC': [13],
        'ACT': [14],
    }
    query_positions = {
        'ACG': [0],
        'CGT': [1],
        'GTA': [2],
        'TAC': [3],
        'ACT': [4],
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
    merged = py_merge_kmer_runs({'ACGTA': [3]}, {'ACGTA': [7]}, k)
    assert len(merged) == 1
    qs, qe, ts, te = merged[0]
    assert qs == 7
    assert qe == 7 + k
    assert ts == 3
    assert te == 3 + k


# ---------------------------------------------------------------------------
# py_merge_rev_runs
# ---------------------------------------------------------------------------


def test_merge_rev_single_kmer():
    """Single RC k-mer hit produces one block of length k."""
    k = 4
    # RC of "AAAC" = "GTTT" found at target pos 0
    merged = py_merge_rev_runs({'AAAC': [0]}, {'AAAC': [0]}, k)
    assert len(merged) == 1
    qs, qe, ts, te = merged[0]
    assert qs == 0
    assert qe == k
    assert ts == 0
    assert te == k


def test_merge_rev_consecutive_merge():
    """Consecutive anti-diagonal RC k-mer hits merge into one block.

    query='AAACCC' (k=3), target='GGGTTT' (= RC of query).
    RC positions in target:
      "AAA" at q=0 → RC "TTT" at t=3
      "AAC" at q=1 → RC "GTT" at t=2
      "ACC" at q=2 → RC "GGT" at t=1
      "CCC" at q=3 → RC "GGG" at t=0
    All share anti-diagonal q+t=3 → merge into (0, 6, 0, 6).
    """
    k = 3
    target_rev = {'AAA': [3], 'AAC': [2], 'ACC': [1], 'CCC': [0]}
    query_pos = {'AAA': [0], 'AAC': [1], 'ACC': [2], 'CCC': [3]}
    merged = py_merge_rev_runs(target_rev, query_pos, k)
    assert len(merged) == 1, f'expected 1 merged block, got {merged}'
    qs, qe, ts, te = merged[0]
    assert qs == 0
    assert qe == 6   # q_prev=3, 3+k=6
    assert ts == 0   # min RC position in target
    assert te == 6   # max RC position + k = 3+3=6


def test_merge_rev_non_consecutive_stays_separate():
    """Non-consecutive RC hits on different anti-diagonals remain separate."""
    k = 4
    # (q=0, t=0): q+t=0  and  (q=0, t=5): q+t=5 → different anti-diagonals
    merged = py_merge_rev_runs({'AAAC': [0], 'CCCC': [5]}, {'AAAC': [0], 'CCCC': [0]}, k)
    assert len(merged) == 2


def test_merge_rev_parallel_antidiagonals():
    """Two parallel RC anti-diagonals each merge independently.

    anti-diag q+t=5: (q=0,t=5), (q=1,t=4) → block (0, 5, 4, 9)
    anti-diag q+t=10: (q=0,t=10), (q=1,t=9) → block (0, 5, 9, 14)
    """
    k = 4
    target_rev = {'AAAA': [5, 10], 'CCCC': [4, 9]}
    query_pos = {'AAAA': [0], 'CCCC': [1]}
    merged = py_merge_rev_runs(target_rev, query_pos, k)
    assert len(merged) == 2, f'expected 2 blocks, got {merged}'
    result_set = set(merged)
    assert (0, 1 + k, 4, 5 + k) in result_set, f'missing anti-diag-5 block in {result_set}'
    assert (0, 1 + k, 9, 10 + k) in result_set, f'missing anti-diag-10 block in {result_set}'


def test_merge_rev_longer_run():
    """Five consecutive RC k-mer hits merge into one block."""
    k = 3
    # 5-kmer RC alignment: q_positions 0-4, t_positions 8-4 (decreasing)
    target_rev = {'K0': [8], 'K1': [7], 'K2': [6], 'K3': [5], 'K4': [4]}
    query_pos = {'K0': [0], 'K1': [1], 'K2': [2], 'K3': [3], 'K4': [4]}
    merged = py_merge_rev_runs(target_rev, query_pos, k)
    assert len(merged) == 1, f'expected 1 merged block, got {merged}'
    qs, qe, ts, te = merged[0]
    assert qs == 0
    assert qe == 4 + k   # q_prev=4, 4+3=7
    assert ts == 4        # min RC position
    assert te == 8 + k   # max RC position + k = 8+3=11


def test_merge_rev_empty_input():
    """Empty input returns empty list."""
    assert py_merge_rev_runs({}, {}, 4) == []
