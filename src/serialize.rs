//! FM-index serialization and deserialization using serde + bincode.
//!
//! Rather than serializing the FM-index data structures directly (which may
//! not implement serde traits), we store only the original sequence text
//! and k-mer sets, then rebuild the FM-index in memory on load.

use crate::error::RustyDotError;
use crate::kmer::{build_kmer_set, sequence_to_index_text, FmIdx};
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// Serializable representation of a single indexed sequence.
///
/// Stores only the original sequence text so that the FM-index can be
/// reconstructed deterministically on load without needing to persist the
/// internal BWT/Occ structures.
#[derive(Serialize, Deserialize)]
pub struct SerializableSequence {
    /// The original sequence bytes (uppercase DNA, without the sentinel).
    pub seq_bytes: Vec<u8>,
}

/// A serializable collection of indexed sequences.
#[derive(Serialize, Deserialize)]
pub struct IndexCollection {
    /// Map from sequence name to its serializable representation.
    pub sequences: HashMap<String, SerializableSequence>,
    /// Map from sequence name to k-mer set.
    pub kmer_sets: HashMap<String, Vec<String>>,
    /// The k-mer length used.
    pub k: usize,
}

/// Save an `IndexCollection` to a binary file using bincode.
///
/// # Arguments
///
/// * `collection` - The index collection to serialize.
/// * `path` - Path to the output file.
///
/// # Errors
///
/// Returns `RustyDotError::Serialization` if serialization fails.
pub fn save_index(collection: &IndexCollection, path: &str) -> Result<(), RustyDotError> {
    let file = File::create(Path::new(path))
        .map_err(|e| RustyDotError::Serialization(e.to_string()))?;
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, collection)
        .map_err(|e| RustyDotError::Serialization(e.to_string()))?;
    Ok(())
}

/// Load an `IndexCollection` from a binary file using bincode.
///
/// # Arguments
///
/// * `path` - Path to the serialized index file.
///
/// # Returns
///
/// The deserialized `IndexCollection`.
///
/// # Errors
///
/// Returns `RustyDotError::Serialization` if deserialization fails.
pub fn load_index(path: &str) -> Result<IndexCollection, RustyDotError> {
    let file =
        File::open(Path::new(path)).map_err(|e| RustyDotError::Serialization(e.to_string()))?;
    let reader = BufReader::new(file);
    bincode::deserialize_from(reader)
        .map_err(|e| RustyDotError::Serialization(e.to_string()))
}

/// Rebuild an `FmIdx` from stored sequence bytes.
///
/// # Arguments
///
/// * `seq_bytes` - The original sequence bytes (without sentinel).
///
/// # Returns
///
/// A usable `FmIdx`.
///
/// # Errors
///
/// Returns `RustyDotError` if FM-index construction fails.
pub fn rebuild_fm_from_bytes(seq_bytes: &[u8]) -> Result<FmIdx, RustyDotError> {
    let seq_str = String::from_utf8_lossy(seq_bytes);
    let text = sequence_to_index_text(&seq_str);
    FmIdx::new(text)
}

/// Python binding: save an index collection to a file.
///
/// Parameters
/// ----------
/// path : str
///     Path to save the serialized index.
/// sequences : dict[str, str]
///     Dictionary of sequence name to sequence string.
/// k : int
///     The k-mer length used to build the index.
///
/// Raises
/// ------
/// ValueError
///     If serialization fails.
#[pyfunction]
pub fn py_save_index(path: &str, sequences: HashMap<String, String>, k: usize) -> PyResult<()> {
    let mut seq_map = HashMap::new();
    let mut kmer_sets_map = HashMap::new();

    for (name, seq) in &sequences {
        let kmer_set = build_kmer_set(seq, k)
            .map_err(|e| -> pyo3::PyErr { e.into() })?;
        seq_map.insert(
            name.clone(),
            SerializableSequence {
                seq_bytes: seq.as_bytes().to_vec(),
            },
        );
        kmer_sets_map.insert(name.clone(), kmer_set.into_iter().collect());
    }

    let collection = IndexCollection {
        sequences: seq_map,
        kmer_sets: kmer_sets_map,
        k,
    };

    save_index(&collection, path).map_err(|e| -> pyo3::PyErr { e.into() })
}

/// Python binding: load an index collection from a file.
///
/// Parameters
/// ----------
/// path : str
///     Path to the serialized index file.
///
/// Returns
/// -------
/// tuple[dict[str, list[str]], int]
///     A tuple of (kmer_sets_dict, k) where kmer_sets_dict maps
///     sequence names to their k-mer lists.
///
/// Raises
/// ------
/// ValueError
///     If deserialization fails.
#[pyfunction]
pub fn py_load_index(path: &str) -> PyResult<(HashMap<String, Vec<String>>, usize)> {
    let collection = load_index(path).map_err(|e| -> pyo3::PyErr { e.into() })?;
    Ok((collection.kmer_sets, collection.k))
}
