//! K-mer set building and coordinate lookup using the FM-index.
// pyo3 pyfunction return types trigger a false-positive useless_conversion lint.
#![allow(clippy::useless_conversion)]

use crate::error::RustyDotError;
use crate::strand::revcomp;
use ahash::AHashSet;
use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use pyo3::prelude::*;
use std::collections::{HashMap, HashSet};

/// A simple FM-index wrapper holding the necessary data structures.
///
/// The suffix array is computed once at construction time and stored so that
/// `find` queries can be answered in O(k + m) time (where `k` is the pattern
/// length and `m` is the number of occurrences) without rebuilding the SA.
pub struct FmIdx {
    /// The BWT of the text.
    pub bwt_data: Vec<u8>,
    /// The less array.
    pub less_data: Vec<usize>,
    /// The Occ table.
    pub occ_data: Occ,
    /// Precomputed suffix array (never serialised; rebuilt from bytes on load).
    pub sa: Vec<usize>,
    /// Length of the indexed text including the sentinel character.
    pub text_len: usize,
    /// Sampling rate for the occurrence table (stored as u32 to match bio API).
    pub occ_sampling: u32,
}

impl FmIdx {
    /// Build an FM-index from a text sequence.
    ///
    /// The text must end with a sentinel character `$`.
    /// The suffix array is computed once and stored for reuse across queries.
    /// The alphabet used is the DNA alphabet (A, C, G, T, N, $).
    ///
    /// # Arguments
    ///
    /// * `text` - The sequence bytes (uppercase DNA + sentinel).  Consumed by
    ///   this function; only derived data structures are retained.
    ///
    /// # Returns
    ///
    /// A new `FmIdx` instance.
    ///
    /// # Errors
    ///
    /// Returns a `RustyDotError` if construction fails.
    pub fn new(text: Vec<u8>) -> Result<Self, RustyDotError> {
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&text);
        let bwt_data = bwt(&text, &sa);
        let less_data = less(&bwt_data, &alphabet);
        let occ_sampling: u32 = 3;
        let occ_data = Occ::new(&bwt_data, occ_sampling, &alphabet);
        let text_len = text.len();

        Ok(FmIdx {
            bwt_data,
            less_data,
            occ_data,
            sa,
            text_len,
            occ_sampling,
        })
    }

    /// Query the FM-index for all occurrences of a pattern.
    ///
    /// Uses the precomputed suffix array stored at construction time, so each
    /// call runs in O(k + m) rather than O(n log n + m).
    ///
    /// # Arguments
    ///
    /// * `pattern` - The query k-mer (bytes).
    ///
    /// # Returns
    ///
    /// A `Vec<usize>` of 0-based start positions in the original text.
    pub fn find(&self, pattern: &[u8]) -> Vec<usize> {
        let fm = FMIndex::new(&self.bwt_data, &self.less_data, &self.occ_data);
        match fm.backward_search(pattern.iter()) {
            BackwardSearchResult::Complete(interval) => interval
                .occ(&self.sa)
                .into_iter()
                .filter(|&pos| pos + pattern.len() <= self.text_len)
                .collect(),
            BackwardSearchResult::Partial(_, _) | BackwardSearchResult::Absent => Vec::new(),
        }
    }
}

/// Build a set of all unique k-mers present in a DNA sequence.
///
/// Scans the sequence with a sliding window of size `k` and collects
/// all k-mers that consist only of valid DNA characters (A, C, G, T).
///
/// # Arguments
///
/// * `seq` - The DNA sequence string.
/// * `k` - The k-mer length.
///
/// # Returns
///
/// A `AHashSet<String>` of unique k-mer strings.
///
/// # Errors
///
/// Returns `RustyDotError::InvalidKmerLength` if `k == 0`.
pub fn build_kmer_set(seq: &str, k: usize) -> Result<AHashSet<String>, RustyDotError> {
    if k == 0 {
        return Err(RustyDotError::InvalidKmerLength(k));
    }
    let bytes = seq.as_bytes();
    if k > bytes.len() {
        return Ok(AHashSet::new());
    }

    let mut set = AHashSet::new();
    for i in 0..=(bytes.len() - k) {
        let kmer = &bytes[i..i + k];
        if kmer.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            set.insert(String::from_utf8_lossy(kmer).into_owned());
        }
    }
    Ok(set)
}

/// Build an FM-index text for a single sequence.
///
/// Appends the sentinel byte `b'$'` to allow unambiguous indexing.
///
/// # Arguments
///
/// * `seq` - The DNA sequence string.
///
/// # Returns
///
/// The byte vector for FM-index construction (sequence + sentinel).
pub fn sequence_to_index_text(seq: &str) -> Vec<u8> {
    let mut text: Vec<u8> = seq.bytes().collect();
    text.push(b'$');
    text
}

/// Find all positions of each k-mer from a set in a sequence using the FM-index.
///
/// # Arguments
///
/// * `kmer_set` - The set of k-mers to look up.
/// * `fm` - The FM-index for the target sequence.
///
/// # Returns
///
/// A `HashMap` mapping each k-mer string to its list of 0-based start positions.
pub fn find_kmer_coords_in_index(
    kmer_set: &AHashSet<String>,
    fm: &FmIdx,
) -> HashMap<String, Vec<usize>> {
    let mut coords: HashMap<String, Vec<usize>> = HashMap::new();
    for kmer in kmer_set {
        let positions = fm.find(kmer.as_bytes());
        if !positions.is_empty() {
            coords.insert(kmer.clone(), positions);
        }
    }
    coords
}

/// Find positions of the reverse complement of each k-mer in the FM-index.
///
/// For each k-mer in `kmer_set`, searches the target FM-index for the reverse
/// complement of that k-mer.  This identifies positions where a k-mer in the
/// query matches the reverse-complement strand of the target sequence.
///
/// The returned map uses the *original* (non-complemented) k-mer as key, so
/// that callers can correlate results with query k-mer positions without needing
/// to recompute reverse complements.
///
/// # Arguments
///
/// * `kmer_set` - The set of query k-mers to look up by reverse complement.
/// * `fm` - The FM-index for the target sequence.
///
/// # Returns
///
/// A `HashMap` mapping each query k-mer string to the 0-based start positions
/// of its reverse complement in the target.  K-mers whose RC is not found are
/// omitted.
pub fn find_rev_coords_in_index(
    kmer_set: &AHashSet<String>,
    fm: &FmIdx,
) -> HashMap<String, Vec<usize>> {
    let mut coords: HashMap<String, Vec<usize>> = HashMap::new();
    for kmer in kmer_set {
        let rc = revcomp(kmer.as_bytes());
        let positions = fm.find(&rc);
        if !positions.is_empty() {
            coords.insert(kmer.clone(), positions);
        }
    }
    coords
}

/// Python binding: build the set of unique k-mers in a sequence.
///
/// Parameters
/// ----------
/// seq : str
///     The DNA sequence string (uppercase recommended).
/// k : int
///     The k-mer length.
///
/// Returns
/// -------
/// set[str]
///     Set of unique k-mer strings found in the sequence.
///
/// Raises
/// ------
/// ValueError
///     If k is 0 or the sequence is empty.
#[pyfunction]
pub fn py_build_kmer_set(seq: &str, k: usize) -> PyResult<HashSet<String>> {
    let set = build_kmer_set(seq, k).map_err(|e| -> pyo3::PyErr { e.into() })?;
    Ok(set.into_iter().collect())
}

/// Python binding: find all positions of each k-mer in a sequence.
///
/// Parameters
/// ----------
/// seq : str
///     The DNA sequence to search in.
/// kmers : list[str]
///     List of k-mer strings to search for.
///
/// Returns
/// -------
/// dict[str, list[int]]
///     Dictionary mapping each k-mer to its list of start positions (0-based).
///
/// Raises
/// ------
/// ValueError
///     If the sequence is invalid or k-mers are inconsistent.
#[pyfunction]
pub fn py_find_kmer_coords(seq: &str, kmers: Vec<String>) -> PyResult<HashMap<String, Vec<usize>>> {
    let text = sequence_to_index_text(seq);
    let fm = FmIdx::new(text).map_err(|e| -> pyo3::PyErr { e.into() })?;

    let kmer_set: AHashSet<String> = kmers.into_iter().collect();
    let coords = find_kmer_coords_in_index(&kmer_set, &fm);
    Ok(coords)
}
