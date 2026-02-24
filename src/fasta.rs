//! FASTA/gzipped FASTA parsing utilities using needletail.
//!
//! This module is only compiled when the `fasta` feature is enabled (the
//! default for native builds).  Wasm builds omit it because needletail
//! cannot link against the `wasm32-unknown-emscripten` target.
// pyo3 pyfunction return types trigger a false-positive useless_conversion lint.
#![allow(clippy::useless_conversion)]

#[cfg(feature = "fasta")]
use crate::error::RustyDotError;
#[cfg(feature = "fasta")]
use needletail::{parse_fastx_file, parse_fastx_stdin};
#[cfg(feature = "fasta")]
use pyo3::prelude::*;
#[cfg(feature = "fasta")]
use std::collections::HashMap;
#[cfg(feature = "fasta")]
use std::path::Path;

/// Read sequences from a FASTA (or gzipped FASTA) file.
///
/// Uses needletail for fast parsing. Gzip decompression is handled
/// automatically when the file extension ends in `.gz`.
///
/// Sequences are returned in the order they appear in the file.
///
/// # Arguments
///
/// * `path` - Path to the FASTA or FASTA.gz file.
///
/// # Returns
///
/// A `Vec` of `(name, sequence)` pairs in file order.
///
/// # Errors
///
/// Returns a `RustyDotError` if the file cannot be opened, parsed, or
/// contains duplicate sequence names.
#[cfg(feature = "fasta")]
pub fn read_fasta(path: &str) -> Result<Vec<(String, String)>, RustyDotError> {
    let mut seqs: Vec<(String, String)> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();

    let mut reader =
        parse_fastx_file(Path::new(path)).map_err(|e| RustyDotError::FastaParse(e.to_string()))?;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| RustyDotError::FastaParse(e.to_string()))?;
        let name = String::from_utf8_lossy(record.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = String::from_utf8_lossy(&record.seq()).to_uppercase();
        if !seen.insert(name.clone()) {
            return Err(RustyDotError::FastaParse(format!(
                "duplicate sequence name '{name}' in FASTA file '{path}'"
            )));
        }
        seqs.push((name, seq));
    }

    Ok(seqs)
}

/// Read sequences from stdin in FASTA format.
///
/// # Returns
///
/// A `HashMap` mapping sequence names (String) to their sequences (String).
///
/// # Errors
///
/// Returns a `RustyDotError` if stdin cannot be parsed.
#[cfg(feature = "fasta")]
pub fn read_fasta_stdin() -> Result<HashMap<String, String>, RustyDotError> {
    let mut seqs: HashMap<String, String> = HashMap::new();

    let mut reader = parse_fastx_stdin().map_err(|e| RustyDotError::FastaParse(e.to_string()))?;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| RustyDotError::FastaParse(e.to_string()))?;
        let name = String::from_utf8_lossy(record.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = String::from_utf8_lossy(&record.seq()).to_uppercase();
        seqs.insert(name, seq);
    }

    Ok(seqs)
}

/// Python binding: read sequences from a FASTA or gzipped FASTA file.
///
/// Parameters
/// ----------
/// path : str
///     Path to the FASTA or FASTA.gz file.
///
/// Returns
/// -------
/// dict[str, str]
///     Dictionary mapping sequence name to sequence string.
///
/// Raises
/// ------
/// ValueError
///     If the file cannot be opened, parsed, or contains duplicate sequence
///     names.
#[cfg(feature = "fasta")]
#[pyfunction]
pub fn py_read_fasta(path: &str) -> PyResult<HashMap<String, String>> {
    let seqs = read_fasta(path).map_err(|e| -> PyErr { e.into() })?;
    Ok(seqs.into_iter().collect())
}
