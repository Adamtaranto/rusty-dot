//! PAF (Pairwise mApping Format) output utilities.
//!
//! PAF format columns (tab-separated):
//! 1.  Query sequence name
//! 2.  Query length
//! 3.  Query start (0-based)
//! 4.  Query end (exclusive)
//! 5.  Strand (+/-)
//! 6.  Target sequence name
//! 7.  Target length
//! 8.  Target start (0-based)
//! 9.  Target end (exclusive)
//! 10. Number of residue matches
//! 11. Alignment block length
//! 12. Mapping quality (255 = missing)

use crate::merge::CoordPair;
use pyo3::prelude::*;

/// Convert a list of `CoordPair` matches to PAF format strings.
///
/// # Arguments
///
/// * `matches` - Slice of `CoordPair` structs representing merged k-mer runs.
/// * `query_name` - Name of the query sequence.
/// * `query_len` - Total length of the query sequence.
/// * `target_name` - Name of the target sequence.
/// * `target_len` - Total length of the target sequence.
///
/// # Returns
///
/// A `Vec<String>` where each element is one PAF line.
pub fn coords_to_paf(
    matches: &[CoordPair],
    query_name: &str,
    query_len: usize,
    target_name: &str,
    target_len: usize,
) -> Vec<String> {
    matches
        .iter()
        .map(|m| {
            let match_len = m.query_end - m.query_start;
            let strand = m.strand as char;
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255",
                query_name,
                query_len,
                m.query_start,
                m.query_end,
                strand,
                target_name,
                target_len,
                m.target_start,
                m.target_end,
                match_len,
                match_len,
            )
        })
        .collect()
}

/// Python binding: convert coordinate tuples to PAF format strings.
///
/// Parameters
/// ----------
/// matches : list[tuple[int, int, int, int]]
///     List of (query_start, query_end, target_start, target_end) tuples.
/// query_name : str
///     Name of the query sequence.
/// query_len : int
///     Total length of the query sequence.
/// target_name : str
///     Name of the target sequence.
/// target_len : int
///     Total length of the target sequence.
///
/// Returns
/// -------
/// list[str]
///     List of PAF format lines (tab-separated).
#[pyfunction]
pub fn py_coords_to_paf(
    matches: Vec<(usize, usize, usize, usize)>,
    query_name: &str,
    query_len: usize,
    target_name: &str,
    target_len: usize,
) -> PyResult<Vec<String>> {
    let coord_pairs: Vec<CoordPair> = matches
        .into_iter()
        .map(|(qs, qe, ts, te)| CoordPair {
            query_start: qs,
            query_end: qe,
            target_start: ts,
            target_end: te,
            strand: crate::strand::STRAND_FWD,
        })
        .collect();
    Ok(coords_to_paf(
        &coord_pairs,
        query_name,
        query_len,
        target_name,
        target_len,
    ))
}
