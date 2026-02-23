//! Functions for merging sequential k-mer coordinate runs.
//!
//! When k-mers are consecutive in both query and target sequences,
//! they can be merged into a single match region, reducing the
//! number of reported features.
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

/// Public alias kept for backward compatibility: forward-strand merge only.
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
