//! rusty-dot: Fast dot plot comparisons of DNA sequences using a canonical ntHash k-mer index
//!
//! This library provides:
//! - FASTA/gzipped FASTA parsing via needletail
//! - Canonical-hash ntHash k-mer index construction
//! - K-mer set building and lookup
//! - Sequential k-mer run merging
//! - PAF format output
//! - Index serialization/deserialization via serde (sequence bytes; indexes rebuilt on load)
//! - Python bindings via PyO3

use pyo3::prelude::*;

pub mod error;
pub mod fasta;
pub mod index;
pub mod kmer;
pub mod kmer_hash;
pub mod merge;
pub mod paf;
pub mod serialize;
pub mod strand;

use index::SequenceIndex;

/// Python module for rusty-dot.
///
/// Exposes all public functions and the SequenceIndex class to Python.
///
/// `gil_used = false` declares the module safe for free-threaded CPython
/// (3.13t/3.14t): the interpreter leaves the GIL disabled instead of switching
/// it back on when this module is imported.  pyo3 0.29 already defaults to
/// this, so the attribute is written out to pin the behaviour rather than
/// inherit it -- the published `cp314t` wheel would otherwise be free-threaded
/// or not depending on a pyo3 upgrade.  `python-tests-freethreaded` in CI
/// asserts `sys._is_gil_enabled()` is false after import.
///
/// The declaration is sound because nothing here holds GIL-bound state: no
/// `static mut`, no `GILOnceCell`/`GILProtected`, no `thread_local`, no
/// `unsafe`, and no `#[pyclass(unsendable)]`.  `SequenceIndex` owns its data
/// outright (plain `HashMap`s, hence `Send + Sync`) and the long-running Rust
/// loops already release via the GIL-agnostic `py.detach(...)`.
///
/// Its mutating methods take `&mut self`, so PyO3's runtime borrow checker
/// raises `RuntimeError: Already borrowed` if two threads mutate one index
/// concurrently, rather than blocking.  That is intended: sharing a single
/// mutable index across threads is a caller error, not something to lock
/// around.  Anything added here that does hold cross-call state must be made
/// thread-safe, or this attribute set back to `true`.
#[pymodule(gil_used = false)]
fn _rusty_dot(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SequenceIndex>()?;
    #[cfg(feature = "fasta")]
    m.add_function(wrap_pyfunction!(fasta::py_read_fasta, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::py_build_kmer_set, m)?)?;
    m.add_function(wrap_pyfunction!(kmer::py_find_kmer_coords, m)?)?;
    m.add_function(wrap_pyfunction!(merge::py_merge_kmer_runs, m)?)?;
    m.add_function(wrap_pyfunction!(merge::py_merge_rev_runs, m)?)?;
    m.add_function(wrap_pyfunction!(merge::py_merge_rev_fwd_runs, m)?)?;
    m.add_function(wrap_pyfunction!(merge::py_merge_runs, m)?)?;
    m.add_function(wrap_pyfunction!(paf::py_coords_to_paf, m)?)?;
    m.add_function(wrap_pyfunction!(serialize::py_save_index, m)?)?;
    m.add_function(wrap_pyfunction!(serialize::py_load_index, m)?)?;
    Ok(())
}
