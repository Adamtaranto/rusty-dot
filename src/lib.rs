//! rusty-dot: Fast dot plot comparisons of DNA sequences using an FM-Index.
//!
//! This library provides:
//! - FASTA/gzipped FASTA parsing via needletail
//! - FM-index construction via rust-bio
//! - K-mer set building and lookup
//! - Sequential k-mer run merging
//! - PAF format output
//! - FM-index serialization/deserialization via serde
//! - Python bindings via PyO3

use pyo3::prelude::*;

pub mod error;
pub mod fasta;
pub mod index;
pub mod kmer;
pub mod merge;
pub mod paf;
pub mod serialize;
pub mod strand;

use index::SequenceIndex;

/// Python module for rusty-dot.
///
/// Exposes all public functions and the SequenceIndex class to Python.
#[pymodule]
fn _rusty_dot(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SequenceIndex>()?;
    m.add_function(wrap_pyfunction!(fasta::py_read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::py_build_kmer_set, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::py_find_kmer_coords, m)?)?;
    m.add_function(wrap_pyfunction!(merge::py_merge_kmer_runs, m)?)?;
    m.add_function(wrap_pyfunction!(paf::py_coords_to_paf, m)?)?;
    m.add_function(wrap_pyfunction!(serialize::py_save_index, m)?)?;
    m.add_function(wrap_pyfunction!(serialize::py_load_index, m)?)?;
    Ok(())
}
