//! Functions for merging sequential k-mer coordinate runs.
//!
//! When k-mers are consecutive in both query and target sequences,
//! they can be merged into a single match region, reducing the
//! number of reported features.
// pyo3 pyfunction return types trigger a false-positive useless_conversion lint.
#![allow(clippy::useless_conversion)]
//!
//! Two merging algorithms are provided:
//! - [`merge_fwd_runs`]: for forward (`+`) strand co-linear hits where
//!   consecutive pairs satisfy `q[i+1] = q[i] + 1`, `t[i+1] = t[i] + 1`.
//! - [`merge_rev_runs`]: for reverse-complement (`-`) strand anti-diagonal
//!   hits where consecutive pairs satisfy `q[i+1] = q[i] + 1`, `t[i+1] = t[i] - 1`.

use crate::strand::{STRAND_FWD, STRAND_REV};
use pyo3::prelude::*;
use std::collections::HashMap;

/// A matched coordinate pair between query and target positions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordPair {
    /// Start position in the query sequence (0-based, inclusive).
    pub query_start: usize,
    /// End position in the query sequence (0-based, exclusive).
    pub query_end: usize,
    /// Start position in the target sequence (0-based, inclusive).
    pub target_start: usize,
    /// End position in the target sequence (0-based, exclusive).
    pub target_end: usize,
    /// Strand of the match: `b'+'` (forward) or `b'-'` (reverse complement).
    pub strand: u8,
}

/// Merge forward-strand (`+`) consecutive k-mer matches into contiguous runs.
///
/// Collects all `(query_pos, target_pos)` pairs, groups them by *diagonal*
/// (`target_pos as i64 − query_pos as i64`), and within each diagonal merges
/// runs where `query_pos[i+1] == query_pos[i] + 1`.
///
/// Sorting by `(diagonal, query_pos)` is essential to avoid interleaving hits
/// from parallel diagonals that share the same query position.
///
/// # Arguments
///
/// * `kmer_coords` - Map from k-mer string to list of target positions.
/// * `query_kmer_positions` - Map from k-mer string to list of query positions.
/// * `k` - The k-mer length (used to compute end positions).
///
/// # Returns
///
/// A `Vec<CoordPair>` of merged `+`-strand match regions.
pub fn merge_fwd_runs(
    kmer_coords: &HashMap<String, Vec<usize>>,
    query_kmer_positions: &HashMap<String, Vec<usize>>,
    k: usize,
) -> Vec<CoordPair> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    for (kmer, target_positions) in kmer_coords {
        if let Some(query_positions) = query_kmer_positions.get(kmer) {
            for &qp in query_positions {
                for &tp in target_positions {
                    pairs.push((qp, tp));
                }
            }
        }
    }

    if pairs.is_empty() {
        return Vec::new();
    }

    pairs.sort_unstable_by_key(|&(q, t)| (t as i64 - q as i64, q as i64));
    pairs.dedup();

    let mut merged: Vec<CoordPair> = Vec::new();
    let mut iter = pairs.iter().peekable();

    if let Some(&(first_q, first_t)) = iter.next() {
        let mut q_start = first_q;
        let mut t_start = first_t;
        let mut q_prev = first_q;
        let mut t_prev = first_t;

        for &(q, t) in iter {
            let same_diag = (t as i64 - q as i64) == (t_prev as i64 - q_prev as i64);
            if same_diag && q == q_prev + 1 {
                q_prev = q;
                t_prev = t;
            } else {
                merged.push(CoordPair {
                    query_start: q_start,
                    query_end: q_prev + k,
                    target_start: t_start,
                    target_end: t_prev + k,
                    strand: STRAND_FWD,
                });
                q_start = q;
                t_start = t;
                q_prev = q;
                t_prev = t;
            }
        }
        merged.push(CoordPair {
            query_start: q_start,
            query_end: q_prev + k,
            target_start: t_start,
            target_end: t_prev + k,
            strand: STRAND_FWD,
        });
    }

    merged
}

/// Merge reverse-complement (`-`) strand k-mer matches into contiguous anti-diagonal runs.
///
/// `target_rev_coords` maps each query k-mer to the positions in the target where the
/// *reverse complement* of that k-mer was found.  As the query position advances by 1,
/// the corresponding RC position decreases by 1 (anti-diagonal: `q + t = constant`).
///
/// Target coordinates bracket the RC match region on the *forward* target strand:
/// `target_start = min(t)` in the run, `target_end = max(t) + k`.
///
/// # Arguments
///
/// * `target_rev_coords` - Map from query k-mer to positions of its RC in the target.
/// * `query_kmer_positions` - Map from query k-mer to positions in the query.
/// * `k` - The k-mer length.
///
/// # Returns
///
/// A `Vec<CoordPair>` of merged `-`-strand match regions.
pub fn merge_rev_runs(
    target_rev_coords: &HashMap<String, Vec<usize>>,
    query_kmer_positions: &HashMap<String, Vec<usize>>,
    k: usize,
) -> Vec<CoordPair> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    for (kmer, target_positions) in target_rev_coords {
        if let Some(query_positions) = query_kmer_positions.get(kmer) {
            for &qp in query_positions {
                for &tp in target_positions {
                    pairs.push((qp, tp));
                }
            }
        }
    }

    if pairs.is_empty() {
        return Vec::new();
    }

    // Sort by anti-diagonal (q + t = constant), then by ascending query position.
    // Within an anti-diagonal run, q increases and t decreases as we iterate.
    pairs.sort_unstable_by_key(|&(q, t)| (q as i64 + t as i64, q as i64));
    pairs.dedup();

    let mut merged: Vec<CoordPair> = Vec::new();
    let mut iter = pairs.iter().peekable();

    if let Some(&(first_q, first_t)) = iter.next() {
        // q_start/t_start track the first element (q_min, t_max) of the run.
        let mut q_start = first_q;
        let mut t_start = first_t;
        let mut q_prev = first_q;
        let mut t_prev = first_t;

        for &(q, t) in iter {
            let same_antidiag = (q as i64 + t as i64) == (q_prev as i64 + t_prev as i64);
            // Guard t_prev > 0 before subtraction to prevent usize underflow.
            let consecutive = same_antidiag && q == q_prev + 1 && t_prev > 0 && t == t_prev - 1;
            if consecutive {
                q_prev = q;
                t_prev = t;
            } else {
                // t_prev = min t in run (target_start); t_start + k = target_end.
                merged.push(CoordPair {
                    query_start: q_start,
                    query_end: q_prev + k,
                    target_start: t_prev,
                    target_end: t_start + k,
                    strand: STRAND_REV,
                });
                q_start = q;
                t_start = t;
                q_prev = q;
                t_prev = t;
            }
        }
        merged.push(CoordPair {
            query_start: q_start,
            query_end: q_prev + k,
            target_start: t_prev,
            target_end: t_start + k,
            strand: STRAND_REV,
        });
    }

    merged
}

/// Merge reverse-complement (`-`) strand k-mer matches that are co-linear on a
/// forward diagonal.
///
/// This function handles the case where sequential k-mers in the query have
/// their reverse complements at **increasing** positions in the target —
/// i.e., as the query position advances by 1, the forward-strand target
/// position of the RC match also advances by 1.  This pattern is distinct
/// from the anti-diagonal pattern handled by [`merge_rev_runs`] and arises
/// in inverted-repeat contexts where the repeat arms are in the same
/// left-to-right direction (both strands share a common direction of traversal).
///
/// The merging logic is identical to [`merge_fwd_runs`] — grouping by forward
/// diagonal (`target_pos − query_pos = constant`) — but results are labelled
/// with [`STRAND_REV`] because the underlying k-mer matches are still
/// reverse-complement matches.
///
/// # Arguments
///
/// * `target_rev_coords` - Map from query k-mer to positions of its RC in the target
///   (as returned by `find_rev_coords_in_index`).
/// * `query_kmer_positions` - Map from query k-mer to positions in the query.
/// * `k` - The k-mer length.
///
/// # Returns
///
/// A `Vec<CoordPair>` of merged `-`-strand match regions where consecutive
/// RC target positions advance by one per step (forward diagonal).
pub fn merge_rev_fwd_runs(
    target_rev_coords: &HashMap<String, Vec<usize>>,
    query_kmer_positions: &HashMap<String, Vec<usize>>,
    k: usize,
) -> Vec<CoordPair> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();

    for (kmer, target_positions) in target_rev_coords {
        if let Some(query_positions) = query_kmer_positions.get(kmer) {
            for &qp in query_positions {
                for &tp in target_positions {
                    pairs.push((qp, tp));
                }
            }
        }
    }

    if pairs.is_empty() {
        return Vec::new();
    }

    // Sort by forward diagonal (t − q = constant), then by ascending query position.
    pairs.sort_unstable_by_key(|&(q, t)| (t as i64 - q as i64, q as i64));
    pairs.dedup();

    let mut merged: Vec<CoordPair> = Vec::new();
    let mut iter = pairs.iter().peekable();

    if let Some(&(first_q, first_t)) = iter.next() {
        let mut q_start = first_q;
        let mut t_start = first_t;
        let mut q_prev = first_q;
        let mut t_prev = first_t;

        for &(q, t) in iter {
            let same_diag = (t as i64 - q as i64) == (t_prev as i64 - q_prev as i64);
            if same_diag && q == q_prev + 1 {
                q_prev = q;
                t_prev = t;
            } else {
                merged.push(CoordPair {
                    query_start: q_start,
                    query_end: q_prev + k,
                    target_start: t_start,
                    target_end: t_prev + k,
                    strand: STRAND_REV,
                });
                q_start = q;
                t_start = t;
                q_prev = q;
                t_prev = t;
            }
        }
        merged.push(CoordPair {
            query_start: q_start,
            query_end: q_prev + k,
            target_start: t_start,
            target_end: t_prev + k,
            strand: STRAND_REV,
        });
    }

    merged
}

///
/// Equivalent to [`merge_fwd_runs`].
pub fn merge_kmer_runs(
    kmer_coords: &HashMap<String, Vec<usize>>,
    query_kmer_positions: &HashMap<String, Vec<usize>>,
    k: usize,
) -> Vec<CoordPair> {
    merge_fwd_runs(kmer_coords, query_kmer_positions, k)
}

/// Python binding: merge sequential k-mer coordinate runs (forward strand).
///
/// Parameters
/// ----------
/// kmer_coords : dict[str, list[int]]
///     Mapping of k-mer to list of target start positions (0-based).
/// query_kmer_positions : dict[str, list[int]]
///     Mapping of k-mer to list of query start positions (0-based).
/// k : int
///     The k-mer length.
///
/// Returns
/// -------
/// list[tuple[int, int, int, int]]
///     List of (query_start, query_end, target_start, target_end) tuples.
///     Coordinates are 0-based, with end positions exclusive.
#[pyfunction]
pub fn py_merge_kmer_runs(
    kmer_coords: HashMap<String, Vec<usize>>,
    query_kmer_positions: HashMap<String, Vec<usize>>,
    k: usize,
) -> PyResult<Vec<(usize, usize, usize, usize)>> {
    let merged = merge_fwd_runs(&kmer_coords, &query_kmer_positions, k);
    Ok(merged
        .into_iter()
        .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
        .collect())
}

/// Python binding: merge reverse-complement (`-` strand) k-mer coordinate runs.
///
/// For each k-mer, `target_rev_coords` should contain the positions in the
/// target where the **reverse complement** of that k-mer was found (as returned
/// by ``find_rev_coords_in_index``).  Consecutive anti-diagonal pairs —
/// where `query_pos` advances by 1 and the corresponding RC target position
/// decreases by 1 — are merged into a single `CoordPair`.
///
/// Parameters
/// ----------
/// target_rev_coords : dict[str, list[int]]
///     Mapping of k-mer to the 0-based start positions of its reverse complement
///     in the target sequence.
/// query_kmer_positions : dict[str, list[int]]
///     Mapping of k-mer to its 0-based start positions in the query sequence.
/// k : int
///     The k-mer length.
///
/// Returns
/// -------
/// list[tuple[int, int, int, int]]
///     List of ``(query_start, query_end, target_start, target_end)`` tuples
///     representing merged ``-``-strand match regions.  Coordinates are
///     0-based with end positions exclusive.  ``target_start`` and
///     ``target_end`` are the forward-strand boundaries of the RC match region
///     on the target.
#[pyfunction]
pub fn py_merge_rev_runs(
    target_rev_coords: HashMap<String, Vec<usize>>,
    query_kmer_positions: HashMap<String, Vec<usize>>,
    k: usize,
) -> PyResult<Vec<(usize, usize, usize, usize)>> {
    let merged = merge_rev_runs(&target_rev_coords, &query_kmer_positions, k);
    Ok(merged
        .into_iter()
        .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
        .collect())
}

/// Python binding: merge reverse-complement (`-` strand) k-mer coordinate runs
/// that are co-linear on a **forward diagonal** (RC positions increase as query
/// position increases).
///
/// This is the complement of :func:`py_merge_rev_runs`, which handles the
/// anti-diagonal case.  Together they cover all possible orientations of
/// reverse-complement k-mer matches:
///
/// * ``py_merge_rev_runs`` — anti-diagonal: ``q`` advances +1, ``t_rc`` decreases by 1
///   (standard inverted-repeat / reverse-complement alignment).
/// * ``py_merge_rev_fwd_runs`` — forward diagonal: ``q`` advances +1, ``t_rc`` also
///   advances +1 (inverted-repeat case where both arms run in the same direction).
///
/// Parameters
/// ----------
/// target_rev_coords : dict[str, list[int]]
///     Mapping of k-mer to the 0-based start positions of its reverse complement
///     in the target sequence (as returned by ``find_rev_coords_in_index``).
/// query_kmer_positions : dict[str, list[int]]
///     Mapping of k-mer to its 0-based start positions in the query sequence.
/// k : int
///     The k-mer length.
///
/// Returns
/// -------
/// list[tuple[int, int, int, int]]
///     List of ``(query_start, query_end, target_start, target_end)`` tuples
///     representing merged ``-``-strand match regions where RC target positions
///     advance together with query positions.  Coordinates are 0-based with
///     end positions exclusive.
#[pyfunction]
pub fn py_merge_rev_fwd_runs(
    target_rev_coords: HashMap<String, Vec<usize>>,
    query_kmer_positions: HashMap<String, Vec<usize>>,
    k: usize,
) -> PyResult<Vec<(usize, usize, usize, usize)>> {
    let merged = merge_rev_fwd_runs(&target_rev_coords, &query_kmer_positions, k);
    Ok(merged
        .into_iter()
        .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
        .collect())
}

/// Python binding: unified merge for both strand orientations.
///
/// A single entry-point for merging k-mer coordinate runs on either the
/// forward (``"+"```) or reverse-complement (``"-"``) strand.  The function
/// dispatches to :func:`merge_fwd_runs` for the forward strand.  For the
/// reverse strand it applies **both** the anti-diagonal algorithm
/// (:func:`merge_rev_runs`) and the co-diagonal algorithm
/// (:func:`merge_rev_fwd_runs`), combining the results and deduplicating
/// any identical blocks.
///
/// For the **forward** strand, ``kmer_coords`` contains the positions of each
/// k-mer in the *target* sequence and ``query_kmer_positions`` contains the
/// positions of the same k-mers in the *query* sequence.
///
/// For the **reverse** strand, ``kmer_coords`` should contain the positions
/// of the **reverse complement** of each query k-mer in the target (as
/// returned by ``find_rev_coords_in_index``), and ``query_kmer_positions``
/// contains the positions of the original k-mers in the query.
///
/// Parameters
/// ----------
/// kmer_coords : dict[str, list[int]]
///     Mapping of k-mer to 0-based target positions.  For ``strand="-"``,
///     these are the positions of the RC of each k-mer in the target.
/// query_kmer_positions : dict[str, list[int]]
///     Mapping of k-mer to 0-based query positions.
/// k : int
///     The k-mer length.
/// strand : str
///     Orientation of the match: ``"+"`` for forward (co-linear diagonal)
///     or ``"-"`` for reverse-complement (both anti-diagonal and co-diagonal
///     patterns are merged and returned together).
///
/// Returns
/// -------
/// list[tuple[int, int, int, int, str]]
///     List of ``(query_start, query_end, target_start, target_end, strand)``
///     5-tuples.  Coordinates are 0-based; end positions are exclusive.
///     ``strand`` echoes the input argument so callers can mix results from
///     multiple calls without losing orientation information.
///
/// Raises
/// ------
/// ValueError
///     If ``strand`` is neither ``"+"`` nor ``"-"``.
#[pyfunction]
pub fn py_merge_runs(
    kmer_coords: HashMap<String, Vec<usize>>,
    query_kmer_positions: HashMap<String, Vec<usize>>,
    k: usize,
    strand: &str,
) -> PyResult<Vec<(usize, usize, usize, usize, String)>> {
    let merged: Vec<CoordPair> = match strand {
        "+" => merge_fwd_runs(&kmer_coords, &query_kmer_positions, k),
        "-" => {
            // Apply both anti-diagonal and co-diagonal merging; deduplicate
            // identical blocks that can arise when a single RC pair has no
            // neighbours on either diagonal.
            let anti = merge_rev_runs(&kmer_coords, &query_kmer_positions, k);
            let co = merge_rev_fwd_runs(&kmer_coords, &query_kmer_positions, k);
            let mut seen = std::collections::HashSet::new();
            let mut combined: Vec<CoordPair> = Vec::new();
            for block in anti.into_iter().chain(co.into_iter()) {
                let key = (
                    block.query_start,
                    block.query_end,
                    block.target_start,
                    block.target_end,
                );
                if seen.insert(key) {
                    combined.push(block);
                }
            }
            combined
        }
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "strand must be '+' or '-', got {:?}",
                other
            )))
        }
    };
    Ok(merged
        .into_iter()
        .map(|c| {
            (
                c.query_start,
                c.query_end,
                c.target_start,
                c.target_end,
                (c.strand as char).to_string(),
            )
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_map(entries: &[(&str, Vec<usize>)]) -> HashMap<String, Vec<usize>> {
        entries
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect()
    }

    // --- merge_fwd_runs ---

    #[test]
    fn test_fwd_empty() {
        assert!(merge_fwd_runs(&HashMap::new(), &HashMap::new(), 4).is_empty());
    }

    #[test]
    fn test_fwd_single_kmer() {
        let t = make_map(&[("ACGT", vec![10])]);
        let q = make_map(&[("ACGT", vec![3])]);
        let result = merge_fwd_runs(&t, &q, 4);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].query_start, 3);
        assert_eq!(result[0].query_end, 7);
        assert_eq!(result[0].target_start, 10);
        assert_eq!(result[0].target_end, 14);
        assert_eq!(result[0].strand, STRAND_FWD);
    }

    #[test]
    fn test_fwd_consecutive_merge() {
        // (q=0,t=10), (q=1,t=11) → one merged block
        let t = make_map(&[("ACGT", vec![10]), ("CGTA", vec![11])]);
        let q = make_map(&[("ACGT", vec![0]), ("CGTA", vec![1])]);
        let result = merge_fwd_runs(&t, &q, 4);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].query_start, 0);
        assert_eq!(result[0].query_end, 5); // 1 + 4
        assert_eq!(result[0].target_start, 10);
        assert_eq!(result[0].target_end, 15); // 11 + 4
    }

    #[test]
    fn test_fwd_non_consecutive_stays_separate() {
        let t = make_map(&[("ACGT", vec![0]), ("TTTT", vec![20])]);
        let q = make_map(&[("ACGT", vec![0]), ("TTTT", vec![5])]);
        let result = merge_fwd_runs(&t, &q, 4);
        assert_eq!(result.len(), 2);
    }

    // --- merge_rev_runs ---

    #[test]
    fn test_rev_empty() {
        assert!(merge_rev_runs(&HashMap::new(), &HashMap::new(), 4).is_empty());
    }

    #[test]
    fn test_rev_single_kmer() {
        // RC of "AAAC" = "GTTT" found at target pos 0
        let t_rev = make_map(&[("AAAC", vec![0])]);
        let q = make_map(&[("AAAC", vec![0])]);
        let result = merge_rev_runs(&t_rev, &q, 4);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].query_start, 0);
        assert_eq!(result[0].query_end, 4);
        assert_eq!(result[0].target_start, 0);
        assert_eq!(result[0].target_end, 4);
        assert_eq!(result[0].strand, STRAND_REV);
    }

    #[test]
    fn test_rev_consecutive_merge() {
        // query='AAACCC' (k=3), target='GGGTTT' (RC of query)
        // k-mer positions in query and their RC positions in target:
        //   "AAA" at q=0 → RC "TTT" at t=3
        //   "AAC" at q=1 → RC "GTT" at t=2
        //   "ACC" at q=2 → RC "GGT" at t=1
        //   "CCC" at q=3 → RC "GGG" at t=0
        // All on anti-diagonal q+t=3 → merge into (0, 6, 0, 6)
        let t_rev = make_map(&[
            ("AAA", vec![3]),
            ("AAC", vec![2]),
            ("ACC", vec![1]),
            ("CCC", vec![0]),
        ]);
        let q = make_map(&[
            ("AAA", vec![0]),
            ("AAC", vec![1]),
            ("ACC", vec![2]),
            ("CCC", vec![3]),
        ]);
        let result = merge_rev_runs(&t_rev, &q, 3);
        assert_eq!(
            result.len(),
            1,
            "expected one merged RC block, got {:?}",
            result
        );
        assert_eq!(result[0].query_start, 0);
        assert_eq!(result[0].query_end, 6); // q_prev=3, 3+3=6
        assert_eq!(result[0].target_start, 0); // t_prev=0 (min t)
        assert_eq!(result[0].target_end, 6); // t_start=3, 3+3=6
        assert_eq!(result[0].strand, STRAND_REV);
    }

    #[test]
    fn test_rev_non_consecutive_stays_separate() {
        // Two RC hits on different anti-diagonals
        let t_rev = make_map(&[("AAAC", vec![0]), ("CCCC", vec![5])]);
        let q = make_map(&[("AAAC", vec![0]), ("CCCC", vec![0])]);
        let result = merge_rev_runs(&t_rev, &q, 4);
        // q+t=0 and q+t=5 → two separate blocks
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_fwd_does_not_merge_antidiagonal_pairs() {
        // Anti-diagonal pairs: consecutive query positions with DECREASING target positions.
        // These belong to the reverse (-) strand; merge_fwd_runs must NOT merge them.
        //
        // (q=0, t=5): forward diagonal t-q = 5
        // (q=1, t=4): forward diagonal t-q = 3  ← different diagonal
        //
        // They sit on the SAME reverse anti-diagonal (q+t = 5), but that is irrelevant
        // here: merge_fwd_runs groups by t-q and the two pairs land on different
        // forward diagonals, so each becomes its own block.
        let t = make_map(&[("AAAA", vec![5]), ("CCCC", vec![4])]);
        let q = make_map(&[("AAAA", vec![0]), ("CCCC", vec![1])]);
        let result = merge_fwd_runs(&t, &q, 4);
        assert_eq!(
            result.len(),
            2,
            "merge_fwd_runs must not merge anti-diagonal (opposite-strand) pairs: {:?}",
            result
        );
        for block in &result {
            assert_eq!(block.strand, STRAND_FWD);
        }
    }

    #[test]
    fn test_rev_does_not_merge_diagonal_pairs() {
        // Forward-diagonal pairs: consecutive query positions with INCREASING target positions.
        // These belong to the forward (+) strand; merge_rev_runs must NOT merge them.
        //
        // (q=0, t=5): reverse anti-diagonal q+t = 5
        // (q=1, t=6): reverse anti-diagonal q+t = 7  ← different anti-diagonal
        //
        // They sit on the SAME forward diagonal (t-q = 5), but that is irrelevant
        // here: merge_rev_runs groups by q+t and the two pairs land on different
        // anti-diagonals, so each becomes its own block.
        let t_rev = make_map(&[("AAAA", vec![5]), ("CCCC", vec![6])]);
        let q = make_map(&[("AAAA", vec![0]), ("CCCC", vec![1])]);
        let result = merge_rev_runs(&t_rev, &q, 4);
        assert_eq!(
            result.len(),
            2,
            "merge_rev_runs must not merge forward-diagonal (opposite-strand) pairs: {:?}",
            result
        );
        for block in &result {
            assert_eq!(block.strand, STRAND_REV);
        }
    }

    #[test]
    fn test_rev_parallel_antidiagonals_stay_separate() {
        // Two separate RC alignments on different anti-diagonals:
        //   anti-diag q+t=5: (0,5), (1,4)
        //   anti-diag q+t=10: (0,10), (1,9)
        let t_rev = make_map(&[("AAAA", vec![5, 10]), ("CCCC", vec![4, 9])]);
        let q = make_map(&[("AAAA", vec![0]), ("CCCC", vec![1])]);
        let result = merge_rev_runs(&t_rev, &q, 4);
        assert_eq!(
            result.len(),
            2,
            "expected 2 separate merged RC blocks, got {:?}",
            result
        );
        let result_set: std::collections::HashSet<(usize, usize, usize, usize)> = result
            .iter()
            .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
            .collect();
        // anti-diag 5: q=[0,5), t=[4,9) (t_prev=4, t_start+k=5+4=9)
        assert!(
            result_set.contains(&(0, 1 + 4, 4, 5 + 4)),
            "missing anti-diag-5 block in {:?}",
            result_set
        );
        // anti-diag 10: q=[0,5), t=[9,14) (t_prev=9, t_start+k=10+4=14)
        assert!(
            result_set.contains(&(0, 1 + 4, 9, 10 + 4)),
            "missing anti-diag-10 block in {:?}",
            result_set
        );
    }

    // --- merge_rev_fwd_runs ---

    #[test]
    fn test_rev_fwd_empty() {
        assert!(merge_rev_fwd_runs(&HashMap::new(), &HashMap::new(), 4).is_empty());
    }

    #[test]
    fn test_rev_fwd_single_kmer() {
        // Single co-diagonal RC hit: should produce one block of length k, labelled -.
        let t_rev = make_map(&[("AAAC", vec![5])]);
        let q = make_map(&[("AAAC", vec![0])]);
        let result = merge_rev_fwd_runs(&t_rev, &q, 4);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].query_start, 0);
        assert_eq!(result[0].query_end, 4);
        assert_eq!(result[0].target_start, 5);
        assert_eq!(result[0].target_end, 9);
        assert_eq!(result[0].strand, STRAND_REV);
    }

    #[test]
    fn test_rev_fwd_consecutive_merge() {
        // Co-diagonal RC case: query='ACACAC' (k=3), target='GTGTGT'.
        // RC(ACA)=TGT at t=1,3; RC(CAC)=GTG at t=0,2.
        // Forward diagonal t-q=1: pairs (q=0,t=1),(q=1,t=2),(q=2,t=3) → merge
        // Forward diagonal t-q=-1: pairs (q=1,t=0),(q=2,t=1),(q=3,t=2) → merge
        let t_rev = make_map(&[("ACA", vec![1, 3]), ("CAC", vec![0, 2])]);
        let q = make_map(&[("ACA", vec![0, 2]), ("CAC", vec![1, 3])]);
        let result = merge_rev_fwd_runs(&t_rev, &q, 3);
        // All forward-diagonal groups; expect merged blocks (not all 8 singletons)
        assert!(
            result.len() < 8,
            "expected merging to reduce blocks, got {:?}",
            result
        );
        for block in &result {
            assert_eq!(block.strand, STRAND_REV, "all blocks must be STRAND_REV");
        }
        let result_set: std::collections::HashSet<(usize, usize, usize, usize)> = result
            .iter()
            .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
            .collect();
        // Diagonal t-q=1: (0,1),(1,2),(2,3) → (qs=0, qe=2+3=5, ts=1, te=3+3=6)
        assert!(
            result_set.contains(&(0, 5, 1, 6)),
            "expected co-diagonal block (0,5,1,6) in {:?}",
            result_set
        );
        // Diagonal t-q=-1: (1,0),(2,1),(3,2) → (qs=1, qe=3+3=6, ts=0, te=2+3=5)
        assert!(
            result_set.contains(&(1, 6, 0, 5)),
            "expected co-diagonal block (1,6,0,5) in {:?}",
            result_set
        );
    }

    #[test]
    fn test_rev_fwd_does_not_merge_antidiagonal_pairs() {
        // Anti-diagonal pairs (t decreases as q increases) should NOT be merged
        // by merge_rev_fwd_runs (they sit on different forward diagonals).
        let t_rev = make_map(&[("AAAA", vec![5]), ("CCCC", vec![4])]);
        let q = make_map(&[("AAAA", vec![0]), ("CCCC", vec![1])]);
        let result = merge_rev_fwd_runs(&t_rev, &q, 4);
        assert_eq!(
            result.len(),
            2,
            "merge_rev_fwd_runs must not merge anti-diagonal pairs: {:?}",
            result
        );
        for block in &result {
            assert_eq!(block.strand, STRAND_REV);
        }
    }

    #[test]
    fn test_rev_fwd_parallel_codiagonals_stay_separate() {
        // Two separate co-diagonal RC runs on different forward diagonals.
        // Diagonal t-q=5: (0,5),(1,6)
        // Diagonal t-q=10: (0,10),(1,11)
        let t_rev = make_map(&[("AAAA", vec![5, 10]), ("CCCC", vec![6, 11])]);
        let q = make_map(&[("AAAA", vec![0]), ("CCCC", vec![1])]);
        let result = merge_rev_fwd_runs(&t_rev, &q, 4);
        assert_eq!(
            result.len(),
            2,
            "expected 2 merged co-diagonal blocks, got {:?}",
            result
        );
        let result_set: std::collections::HashSet<(usize, usize, usize, usize)> = result
            .iter()
            .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
            .collect();
        assert!(
            result_set.contains(&(0, 1 + 4, 5, 6 + 4)),
            "missing diagonal-5 block in {:?}",
            result_set
        );
        assert!(
            result_set.contains(&(0, 1 + 4, 10, 11 + 4)),
            "missing diagonal-10 block in {:?}",
            result_set
        );
    }
}
