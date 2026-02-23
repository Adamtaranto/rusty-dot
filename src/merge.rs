//! Functions for merging sequential k-mer coordinate runs.
//!
//! When k-mers are consecutive in both query and target sequences,
//! they can be merged into a single match region, reducing the
//! number of reported features.

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
}

/// Merge consecutive k-mer matches into contiguous coordinate runs.
///
/// For each query position, there may be multiple target positions.
/// This function collects all (query_pos, target_pos) pairs for a given
/// k-mer length and merges runs where:
/// - `query_pos[i+1] == query_pos[i] + 1` (consecutive in query)
/// - `target_pos[i+1] == target_pos[i] + 1` (consecutive in target)
///
/// Each run is reported as a `CoordPair` with the merged start/end.
///
/// # Arguments
///
/// * `kmer_coords` - Map from k-mer string to list of target positions.
/// * `query_kmer_positions` - Map from k-mer string to list of query positions.
/// * `k` - The k-mer length (used to compute end positions).
///
/// # Returns
///
/// A sorted `Vec<CoordPair>` of merged match regions.
pub fn merge_kmer_runs(
    kmer_coords: &HashMap<String, Vec<usize>>,
    query_kmer_positions: &HashMap<String, Vec<usize>>,
    k: usize,
) -> Vec<CoordPair> {
    // Collect all (query_pos, target_pos) tuples
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

    // Sort by query position then target position
    pairs.sort_unstable();
    pairs.dedup();

    // Merge consecutive runs
    let mut merged: Vec<CoordPair> = Vec::new();
    let mut iter = pairs.iter().peekable();

    if let Some(&(mut q_start, mut t_start)) = iter.next() {
        let mut q_prev = q_start;
        let mut t_prev = t_start;

        while let Some(&(q, t)) = iter.next() {
            // Check if this pair extends the current run
            if q == q_prev + 1 && t == t_prev + 1 {
                q_prev = q;
                t_prev = t;
            } else {
                // Save the current run
                merged.push(CoordPair {
                    query_start: q_start,
                    query_end: q_prev + k,
                    target_start: t_start,
                    target_end: t_prev + k,
                });
                q_start = q;
                t_start = t;
                q_prev = q;
                t_prev = t;
            }
        }
        // Push the last run
        merged.push(CoordPair {
            query_start: q_start,
            query_end: q_prev + k,
            target_start: t_start,
            target_end: t_prev + k,
        });
    }

    merged
}

/// Python binding: merge sequential k-mer coordinate runs.
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
    let merged = merge_kmer_runs(&kmer_coords, &query_kmer_positions, k);
    Ok(merged
        .into_iter()
        .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
        .collect())
}
