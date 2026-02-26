//! The main `SequenceIndex` class providing the core functionality.
//!
//! Provides FM-index construction, k-mer lookup, sequence comparison,
//! and PAF output, with full Python bindings via PyO3.
// pyo3 pyfunction/pymethods return types trigger a false-positive useless_conversion lint.
#![allow(clippy::useless_conversion)]

use crate::kmer::{
    build_kmer_set, find_kmer_coords_in_index, find_rev_coords_in_index, sequence_to_index_text,
    FmIdx,
};
use crate::merge::{
    merge_fwd_runs, merge_kmer_runs, merge_rev_fwd_runs, merge_rev_runs, CoordPair,
};
use crate::paf::coords_to_paf;
use crate::serialize::{
    load_index, rebuild_fm_from_bytes, save_index, IndexCollection, SerializableSequence,
};
use crate::strand::{revcomp, STRAND_REV};
use ahash::AHashSet;
use pyo3::prelude::*;
use std::collections::HashMap;

/// Stranded match coordinates: (query_start, query_end, target_start, target_end, strand).
type StrandedMatch = (usize, usize, usize, usize, String);

/// In-memory store for a single sequence's index data.
struct SequenceData {
    /// FM-index for this sequence.
    fm: FmIdx,
    /// Set of unique k-mers in this sequence.
    kmer_set: AHashSet<String>,
    /// Original sequence bytes (without sentinel), for serialization.
    seq_bytes: Vec<u8>,
    /// Original sequence length.
    seq_len: usize,
}

/// PyO3-exposed class for building and querying FM-indexes for DNA sequences.
///
/// Each sequence added to the index receives its **own independent FM-index**
/// built by [rust-bio](https://docs.rs/bio).  The rust-bio FM-index cannot be
/// updated or extended after construction, so adding more sequences never
/// modifies an existing FM-index — it only creates a new one.
///
/// The index behaves as a **dictionary of per-sequence FM-indexes**:
///
/// * `add_sequence` / `load_fasta` — **add** new entries to the collection;
///   calling either method multiple times accumulates sequences rather than
///   replacing them.
/// * If `add_sequence` (or `load_fasta`) is called with a name that already
///   exists in the index, the existing entry is **silently overwritten** with a
///   new FM-index for the new sequence.
/// * Pairwise comparisons always operate on exactly two independent FM-indexes.
///
/// The `k` value is fixed at construction time and applies to all sequences.
#[pyclass]
pub struct SequenceIndex {
    /// Map from sequence name to index data.
    sequences: HashMap<String, SequenceData>,
    /// The k-mer length used for indexing.
    k: usize,
    /// Cache of pairwise comparison results: (query, target, merge) -> coord pairs.
    pair_cache: HashMap<(String, String, bool), Vec<CoordPair>>,
}

#[pymethods]
impl SequenceIndex {
    /// Create a new SequenceIndex.
    ///
    /// Parameters
    /// ----------
    /// k : int
    ///     The k-mer length to use for indexing and comparison.
    ///
    /// Returns
    /// -------
    /// SequenceIndex
    ///     A new empty index.
    #[new]
    pub fn new(k: usize) -> PyResult<Self> {
        if k == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "k must be greater than 0",
            ));
        }
        Ok(SequenceIndex {
            sequences: HashMap::new(),
            k,
            pair_cache: HashMap::new(),
        })
    }

    /// Add a single sequence to the index.
    ///
    /// Builds a **new independent FM-index** for `seq` using rust-bio and
    /// stores it alongside the k-mer set and raw sequence bytes.  The
    /// rust-bio FM-index is constructed once and cannot be extended; each
    /// call to `add_sequence` creates a separate FM-index for that sequence
    /// only.
    ///
    /// Calling `add_sequence` does **not** affect any other sequence already
    /// in the index — each sequence has its own isolated FM-index.
    ///
    /// If a sequence with `name` already exists in the index, a
    /// `UserWarning` is emitted and the existing entry is **overwritten**
    /// with a new FM-index for the new `seq`.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     A unique identifier for the sequence. Re-using an existing name
    ///     emits a warning and replaces the previous entry.
    /// seq : str
    ///     The DNA sequence string (uppercase recommended; lowercase is
    ///     accepted and treated as uppercase).
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the FM-index cannot be built (e.g., invalid characters).
    pub fn add_sequence(&mut self, py: Python<'_>, name: &str, seq: &str) -> PyResult<()> {
        if self.sequences.contains_key(name) {
            let warnings = py.import("warnings")?;
            warnings.call_method1(
                "warn",
                (
                    format!(
                        "Sequence name '{name}' already exists in the index and will be overwritten."
                    ),
                    py.get_type::<pyo3::exceptions::PyUserWarning>(),
                ),
            )?;
        }
        let text = sequence_to_index_text(seq);
        let fm = FmIdx::new(text).map_err(|e| -> pyo3::PyErr { e.into() })?;
        let kmer_set = build_kmer_set(seq, self.k).map_err(|e| -> pyo3::PyErr { e.into() })?;
        let seq_len = seq.len();
        let seq_bytes = seq.as_bytes().to_vec();
        self.sequences.insert(
            name.to_string(),
            SequenceData {
                fm,
                kmer_set,
                seq_bytes,
                seq_len,
            },
        );
        Ok(())
    }

    /// Load all sequences from a FASTA or gzipped FASTA file and add them to the index.
    ///
    /// Parses the file with needletail (automatic gzip detection) and
    /// builds a fresh **independent FM-index** for each record.
    ///
    /// Sequences already in the index are **preserved** — `load_fasta` only
    /// adds new entries (or overwrites entries whose name already exists).
    /// Calling `load_fasta` twice on two different files accumulates all
    /// sequences from both files in the same index.
    ///
    /// If the FASTA file contains **duplicate sequence names** (two records
    /// with the same identifier), a `ValueError` is raised before any sequences
    /// are added to the index.
    ///
    /// If a record's name **already exists in the index**, a `UserWarning` is
    /// emitted and the existing entry is overwritten.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to the FASTA or FASTA.gz file.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     List of sequence names that were added (in file order).
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the file cannot be read or parsed, or if it contains duplicate
    ///     sequence names.
    #[cfg(feature = "fasta")]
    pub fn load_fasta(&mut self, py: Python<'_>, path: &str) -> PyResult<Vec<String>> {
        use crate::error::RustyDotError;
        use needletail::parse_fastx_file;
        use std::collections::{HashMap, HashSet};
        use std::path::Path;

        let warnings = py.import("warnings")?;
        let mut names: Vec<String> = Vec::new();
        let mut seen_in_file: HashSet<String> = HashSet::new();
        let mut temp: HashMap<String, SequenceData> = HashMap::new();

        let mut reader = parse_fastx_file(Path::new(path))
            .map_err(|e| -> pyo3::PyErr { RustyDotError::FastaParse(e.to_string()).into() })?;

        while let Some(record) = reader.next() {
            let record = record
                .map_err(|e| -> pyo3::PyErr { RustyDotError::FastaParse(e.to_string()).into() })?;
            let name = String::from_utf8_lossy(record.id())
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let seq = String::from_utf8_lossy(&record.seq()).to_uppercase();

            if !seen_in_file.insert(name.clone()) {
                return Err(RustyDotError::FastaParse(format!(
                    "duplicate sequence name '{name}' in FASTA file '{path}'"
                ))
                .into());
            }
            let text = sequence_to_index_text(&seq);
            let fm = FmIdx::new(text).map_err(|e| -> pyo3::PyErr { e.into() })?;
            let kmer_set = build_kmer_set(&seq, self.k).map_err(|e| -> pyo3::PyErr { e.into() })?;
            let seq_len = seq.len();
            let seq_bytes = seq.as_bytes().to_vec();
            temp.insert(
                name.clone(),
                SequenceData {
                    fm,
                    kmer_set,
                    seq_bytes,
                    seq_len,
                },
            );
            names.push(name);
        }

        // No errors — merge the fully-validated batch into self.sequences.
        for name in &names {
            if self.sequences.contains_key(name.as_str()) {
                warnings.call_method1(
                    "warn",
                    (
                        format!(
                            "Sequence name '{name}' already exists in the index and will be overwritten."
                        ),
                        py.get_type::<pyo3::exceptions::PyUserWarning>(),
                    ),
                )?;
            }
            if let Some(data) = temp.remove(name) {
                self.sequences.insert(name.clone(), data);
            }
        }
        Ok(names)
    }

    /// Get the list of sequence names in the index.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     List of sequence names.
    pub fn sequence_names(&self) -> Vec<String> {
        self.sequences.keys().cloned().collect()
    }

    /// Get the k-mer set for a named sequence.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     The sequence name.
    ///
    /// Returns
    /// -------
    /// set[str]
    ///     The set of unique k-mers in this sequence.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If the sequence name is not found.
    pub fn get_kmer_set(&self, name: &str) -> PyResult<std::collections::HashSet<String>> {
        match self.sequences.get(name) {
            Some(data) => Ok(data.kmer_set.iter().cloned().collect()),
            None => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Sequence '{}' not found in index",
                name
            ))),
        }
    }

    /// Get the length of a named sequence.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     The sequence name.
    ///
    /// Returns
    /// -------
    /// int
    ///     The sequence length.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If the sequence name is not found.
    pub fn get_sequence_length(&self, name: &str) -> PyResult<usize> {
        match self.sequences.get(name) {
            Some(data) => Ok(data.seq_len),
            None => Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Sequence '{}' not found in index",
                name
            ))),
        }
    }

    /// Find the shared k-mers between two sequences and return their coordinates.
    ///
    /// Uses the smaller k-mer set for lookup efficiency.
    /// Results are cached for subsequent queries.
    ///
    /// Parameters
    /// ----------
    /// query_name : str
    ///     Name of the query sequence.
    /// target_name : str
    ///     Name of the target sequence.
    /// merge : bool, optional
    ///     Whether to merge sequential k-mer runs. Default is True.
    ///
    /// Returns
    /// -------
    /// list[tuple[int, int, int, int]]
    ///     List of (query_start, query_end, target_start, target_end) tuples.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If either sequence name is not found.
    #[pyo3(signature = (query_name, target_name, merge=true))]
    pub fn compare_sequences(
        &mut self,
        query_name: &str,
        target_name: &str,
        merge: bool,
    ) -> PyResult<Vec<(usize, usize, usize, usize)>> {
        let cache_key = (query_name.to_string(), target_name.to_string(), merge);

        if let Some(cached) = self.pair_cache.get(&cache_key) {
            return Ok(cached
                .iter()
                .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
                .collect());
        }

        // Check existence
        if !self.sequences.contains_key(query_name) {
            return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Query sequence '{}' not found",
                query_name
            )));
        }
        if !self.sequences.contains_key(target_name) {
            return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                "Target sequence '{}' not found",
                target_name
            )));
        }

        // Determine which k-mer set is smaller (use as probe set)
        let probe_is_query = {
            let q_len = self.sequences[query_name].kmer_set.len();
            let t_len = self.sequences[target_name].kmer_set.len();
            q_len <= t_len
        };

        // Get shared k-mers: probe ∩ other
        let shared_kmers: Vec<String> = {
            let (probe_name, other_name) = if probe_is_query {
                (query_name, target_name)
            } else {
                (target_name, query_name)
            };
            let probe_set = &self.sequences[probe_name].kmer_set;
            let other_set = &self.sequences[other_name].kmer_set;
            probe_set
                .iter()
                .filter(|k| other_set.contains(*k))
                .cloned()
                .collect()
        };

        if shared_kmers.is_empty() {
            self.pair_cache.insert(cache_key, Vec::new());
            return Ok(Vec::new());
        }

        let shared_set: AHashSet<String> = shared_kmers.into_iter().collect();

        // Find coords of shared k-mers in both query and target
        let query_coords = {
            let fm = &self.sequences[query_name].fm;
            find_kmer_coords_in_index(&shared_set, fm)
        };

        let target_coords = {
            let fm = &self.sequences[target_name].fm;
            find_kmer_coords_in_index(&shared_set, fm)
        };

        let result = if merge {
            let merged = merge_kmer_runs(&target_coords, &query_coords, self.k);
            let tuples: Vec<(usize, usize, usize, usize)> = merged
                .iter()
                .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
                .collect();
            self.pair_cache.insert(cache_key, merged);
            tuples
        } else {
            // Return individual k-mer hits (unmerged)
            use crate::strand::STRAND_FWD;
            let mut unmerged: Vec<CoordPair> = Vec::new();
            for (kmer, q_positions) in &query_coords {
                if let Some(t_positions) = target_coords.get(kmer) {
                    for &qp in q_positions {
                        for &tp in t_positions {
                            unmerged.push(CoordPair {
                                query_start: qp,
                                query_end: qp + self.k,
                                target_start: tp,
                                target_end: tp + self.k,
                                strand: STRAND_FWD,
                            });
                        }
                    }
                }
            }
            unmerged.sort_unstable_by_key(|c| (c.query_start, c.target_start));
            let tuples: Vec<(usize, usize, usize, usize)> = unmerged
                .iter()
                .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
                .collect();
            self.pair_cache.insert(cache_key, unmerged);
            tuples
        };

        Ok(result)
    }

    /// Find shared k-mer matches between two sequences, reporting both strands.
    ///
    /// Extends `compare_sequences` to also search for reverse-complement (`-` strand)
    /// matches: k-mers in the query whose reverse complement appears in the target.
    /// Results are cached per `(query, target, merge)` tuple in a separate stranded cache.
    ///
    /// Parameters
    /// ----------
    /// query_name : str
    ///     Name of the query sequence.
    /// target_name : str
    ///     Name of the target sequence.
    /// merge : bool, optional
    ///     Whether to merge co-linear k-mer runs. Default is True.
    ///
    /// Returns
    /// -------
    /// list[tuple[int, int, int, int, str]]
    ///     List of (query_start, query_end, target_start, target_end, strand) tuples.
    ///     Strand is ``"+"`` for forward matches and ``"-"`` for reverse-complement matches.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If either sequence name is not found.
    #[pyo3(signature = (query_name, target_name, merge=true))]
    pub fn compare_sequences_stranded(
        &mut self,
        query_name: &str,
        target_name: &str,
        merge: bool,
    ) -> PyResult<Vec<StrandedMatch>> {
        for name in [query_name, target_name] {
            if !self.sequences.contains_key(name) {
                return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                    "Sequence '{}' not found",
                    name
                )));
            }
        }

        // --- Forward (+ strand) shared k-mers ---
        let fwd_shared: AHashSet<String> = {
            let probe_is_query = self.sequences[query_name].kmer_set.len()
                <= self.sequences[target_name].kmer_set.len();
            let (probe, other) = if probe_is_query {
                (query_name, target_name)
            } else {
                (target_name, query_name)
            };
            let probe_set = &self.sequences[probe].kmer_set;
            let other_set = &self.sequences[other].kmer_set;
            probe_set
                .iter()
                .filter(|k| other_set.contains(*k))
                .cloned()
                .collect()
        };

        // --- Reverse (- strand) k-mers: query kmers whose RC is in target ---
        let rev_shared: AHashSet<String> = {
            let q_set = &self.sequences[query_name].kmer_set;
            let t_set = &self.sequences[target_name].kmer_set;
            q_set
                .iter()
                .filter(|k| {
                    let rc_bytes = revcomp(k.as_bytes());
                    if let Ok(rc_str) = std::str::from_utf8(&rc_bytes) {
                        t_set.contains(rc_str)
                    } else {
                        false
                    }
                })
                .cloned()
                .collect()
        };

        let mut all_pairs: Vec<CoordPair> = Vec::new();

        // + strand hits
        if !fwd_shared.is_empty() {
            let query_fwd = {
                let fm = &self.sequences[query_name].fm;
                find_kmer_coords_in_index(&fwd_shared, fm)
            };
            let target_fwd = {
                let fm = &self.sequences[target_name].fm;
                find_kmer_coords_in_index(&fwd_shared, fm)
            };
            if merge {
                all_pairs.extend(merge_fwd_runs(&target_fwd, &query_fwd, self.k));
            } else {
                use crate::strand::STRAND_FWD;
                for (kmer, q_pos) in &query_fwd {
                    if let Some(t_pos) = target_fwd.get(kmer) {
                        for &qp in q_pos {
                            for &tp in t_pos {
                                all_pairs.push(CoordPair {
                                    query_start: qp,
                                    query_end: qp + self.k,
                                    target_start: tp,
                                    target_end: tp + self.k,
                                    strand: STRAND_FWD,
                                });
                            }
                        }
                    }
                }
            }
        }

        // - strand hits
        if !rev_shared.is_empty() {
            let query_rev = {
                let fm = &self.sequences[query_name].fm;
                find_kmer_coords_in_index(&rev_shared, fm)
            };
            let target_rev = {
                let fm = &self.sequences[target_name].fm;
                find_rev_coords_in_index(&rev_shared, fm)
            };
            if merge {
                // Apply both anti-diagonal and co-diagonal merging for RC hits,
                // deduplicating identical blocks that arise when a single RC pair
                // has no neighbours on either diagonal.
                let anti = merge_rev_runs(&target_rev, &query_rev, self.k);
                let co = merge_rev_fwd_runs(&target_rev, &query_rev, self.k);
                let mut seen = std::collections::HashSet::new();
                for block in anti.into_iter().chain(co.into_iter()) {
                    let key = (
                        block.query_start,
                        block.query_end,
                        block.target_start,
                        block.target_end,
                    );
                    if seen.insert(key) {
                        all_pairs.push(block);
                    }
                }
            } else {
                for (kmer, q_pos) in &query_rev {
                    if let Some(t_pos) = target_rev.get(kmer) {
                        for &qp in q_pos {
                            for &tp in t_pos {
                                all_pairs.push(CoordPair {
                                    query_start: qp,
                                    query_end: qp + self.k,
                                    target_start: tp,
                                    target_end: tp + self.k,
                                    strand: STRAND_REV,
                                });
                            }
                        }
                    }
                }
            }
        }

        Ok(all_pairs
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

    /// Compute the optimal ordering of query and target contigs to maximise collinearity.
    ///
    /// Uses a gravity-centre algorithm inspired by d-genies: for each query contig the
    /// gravity centre is the weighted mean of the target position mid-points across all
    /// matches, normalised by the total target sequence length.  Contigs are then sorted
    /// by ascending gravity centre.
    ///
    /// This method calls ``compare_sequences_stranded`` for every ordered pair (which
    /// uses the internal cache for repeated calls) and then sorts by gravity.
    ///
    /// Parameters
    /// ----------
    /// query_names : list[str]
    ///     Names of query sequences to reorder.  Must all be present in the index.
    /// target_names : list[str]
    ///     Names of target sequences to use as the reference axis.  Must all be in the index.
    ///
    /// Returns
    /// -------
    /// tuple[list[str], list[str]]
    ///     ``(sorted_query_names, sorted_target_names)`` — both lists are reordered so that
    ///     the resulting dotplot shows maximum collinearity along the diagonal.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If any sequence name is not present in the index.
    pub fn optimal_contig_order(
        &mut self,
        query_names: Vec<String>,
        target_names: Vec<String>,
    ) -> PyResult<(Vec<String>, Vec<String>)> {
        for n in query_names.iter().chain(target_names.iter()) {
            if !self.sequences.contains_key(n.as_str()) {
                return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                    "Sequence '{}' not found",
                    n
                )));
            }
        }

        let total_target_len: usize = target_names
            .iter()
            .map(|n| self.sequences[n.as_str()].seq_len)
            .sum::<usize>()
            .max(1);

        // Build cumulative target offsets so we can express each target position as an
        // absolute coordinate across the concatenated target.
        let mut target_offsets: HashMap<String, usize> = HashMap::new();
        let mut offset = 0usize;
        for n in &target_names {
            target_offsets.insert(n.clone(), offset);
            offset += self.sequences[n.as_str()].seq_len;
        }

        // Gravity of each query contig = weighted mean of match mid-points on the target axis.
        let mut q_gravity: Vec<(f64, String)> = Vec::new();
        for q in &query_names {
            let mut weight_sum = 0.0f64;
            let mut weighted_pos = 0.0f64;
            for t in &target_names {
                let matches = self.compare_sequences_stranded(q.as_str(), t.as_str(), true)?;
                let t_offset = *target_offsets.get(t).unwrap_or(&0) as f64;
                for (_, _, ts, te, _) in &matches {
                    let size = (te - ts) as f64;
                    let mid = t_offset + (*ts as f64 + *te as f64) / 2.0;
                    weighted_pos += size * mid;
                    weight_sum += size;
                }
            }
            let gravity = if weight_sum > 0.0 {
                weighted_pos / weight_sum / total_target_len as f64
            } else {
                f64::MAX // unmatched contigs sort to the end
            };
            q_gravity.push((gravity, q.clone()));
        }
        // Sort: matched contigs by ascending gravity, then unmatched by descending length.
        let mut matched_q: Vec<(f64, String)> = q_gravity
            .iter()
            .filter(|(g, _)| *g < f64::MAX)
            .cloned()
            .collect();
        let mut unmatched_q: Vec<(usize, String)> = q_gravity
            .iter()
            .filter(|(g, _)| *g >= f64::MAX)
            .map(|(_, n)| (self.sequences[n.as_str()].seq_len, n.clone()))
            .collect();
        matched_q.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        unmatched_q.sort_by(|a, b| b.0.cmp(&a.0)); // descending length
        let q_gravity: Vec<String> = matched_q
            .into_iter()
            .map(|(_, n)| n)
            .chain(unmatched_q.into_iter().map(|(_, n)| n))
            .collect();

        let total_query_len: usize = query_names
            .iter()
            .map(|n| self.sequences[n.as_str()].seq_len)
            .sum::<usize>()
            .max(1);

        let mut query_offsets: HashMap<String, usize> = HashMap::new();
        let mut q_offset = 0usize;
        for n in &query_names {
            query_offsets.insert(n.clone(), q_offset);
            q_offset += self.sequences[n.as_str()].seq_len;
        }

        let mut t_gravity: Vec<(f64, String)> = Vec::new();
        for t in &target_names {
            let mut weight_sum = 0.0f64;
            let mut weighted_pos = 0.0f64;
            for q in &query_names {
                let matches = self.compare_sequences_stranded(q.as_str(), t.as_str(), true)?;
                let q_offset = *query_offsets.get(q).unwrap_or(&0) as f64;
                for (qs, qe, _, _, _) in &matches {
                    let size = (qe - qs) as f64;
                    let mid = q_offset + (*qs as f64 + *qe as f64) / 2.0;
                    weighted_pos += size * mid;
                    weight_sum += size;
                }
            }
            let gravity = if weight_sum > 0.0 {
                weighted_pos / weight_sum / total_query_len as f64
            } else {
                f64::MAX
            };
            t_gravity.push((gravity, t.clone()));
        }
        // Sort: matched contigs by ascending gravity, then unmatched by descending length.
        let mut matched_t: Vec<(f64, String)> = t_gravity
            .iter()
            .filter(|(g, _)| *g < f64::MAX)
            .cloned()
            .collect();
        let mut unmatched_t: Vec<(usize, String)> = t_gravity
            .iter()
            .filter(|(g, _)| *g >= f64::MAX)
            .map(|(_, n)| (self.sequences[n.as_str()].seq_len, n.clone()))
            .collect();
        matched_t.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        unmatched_t.sort_by(|a, b| b.0.cmp(&a.0)); // descending length
        let t_gravity: Vec<String> = matched_t
            .into_iter()
            .map(|(_, n)| n)
            .chain(unmatched_t.into_iter().map(|(_, n)| n))
            .collect();

        Ok((q_gravity, t_gravity))
    }

    /// Get PAF-formatted strings for a pair of sequences.
    ///
    /// Parameters
    /// ----------
    /// query_name : str
    ///     Name of the query sequence.
    /// target_name : str
    ///     Name of the target sequence.
    /// merge : bool, optional
    ///     Whether to merge sequential k-mer runs. Default is True.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     List of PAF format lines.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If either sequence name is not found.
    #[pyo3(signature = (query_name, target_name, merge=true))]
    pub fn get_paf(
        &mut self,
        query_name: &str,
        target_name: &str,
        merge: bool,
    ) -> PyResult<Vec<String>> {
        let _ = self.compare_sequences(query_name, target_name, merge)?;
        let cache_key = (query_name.to_string(), target_name.to_string(), merge);
        let matches = self
            .pair_cache
            .get(&cache_key)
            .expect("cache was just populated by compare_sequences");
        let query_len = self.sequences[query_name].seq_len;
        let target_len = self.sequences[target_name].seq_len;
        Ok(coords_to_paf(
            matches,
            query_name,
            query_len,
            target_name,
            target_len,
        ))
    }

    /// Return PAF-formatted strings for every ordered sequence pair in the index.
    ///
    /// Calls :meth:`get_paf` for every ``(i, j)`` pair where ``i != j``,
    /// populating the comparison cache as a side-effect.
    ///
    /// Parameters
    /// ----------
    /// merge : bool, optional
    ///     Whether to merge sequential k-mer runs before generating PAF lines.
    ///     Default is ``True``.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     All PAF lines for every pairwise comparison, one line per match.
    #[pyo3(signature = (merge=true))]
    pub fn get_paf_all(&mut self, merge: bool) -> PyResult<Vec<String>> {
        let names: Vec<String> = self.sequences.keys().cloned().collect();
        let mut all_paf: Vec<String> = Vec::new();
        for i in 0..names.len() {
            for j in 0..names.len() {
                if i != j {
                    let q = names[i].clone();
                    let t = names[j].clone();
                    let lines = self.get_paf(&q, &t, merge)?;
                    all_paf.extend(lines);
                }
            }
        }
        Ok(all_paf)
    }

    /// Pre-calculate pairwise comparisons for all sequence pairs.
    ///
    /// Parameters
    /// ----------
    /// merge : bool, optional
    ///     Whether to merge sequential k-mer runs. Default is True.
    ///
    /// Returns
    /// -------
    /// list[tuple[str, str]]
    ///     List of (query_name, target_name) pairs that were computed.
    #[pyo3(signature = (merge=true))]
    pub fn precompute_all_pairs(&mut self, merge: bool) -> PyResult<Vec<(String, String)>> {
        let names: Vec<String> = self.sequences.keys().cloned().collect();
        let mut pairs = Vec::new();
        for i in 0..names.len() {
            for j in 0..names.len() {
                if i != j {
                    let q = names[i].clone();
                    let t = names[j].clone();
                    self.compare_sequences(&q, &t, merge)?;
                    pairs.push((q, t));
                }
            }
        }
        Ok(pairs)
    }

    /// Save the current index to a binary file.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to the output file.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If serialization fails.
    pub fn save(&self, path: &str) -> PyResult<()> {
        let mut seq_map = HashMap::new();
        let mut kmer_sets_map = HashMap::new();

        for (name, data) in &self.sequences {
            seq_map.insert(
                name.clone(),
                SerializableSequence {
                    seq_bytes: data.seq_bytes.clone(),
                },
            );
            kmer_sets_map.insert(
                name.clone(),
                data.kmer_set.iter().cloned().collect::<Vec<_>>(),
            );
        }

        let collection = IndexCollection {
            sequences: seq_map,
            kmer_sets: kmer_sets_map,
            k: self.k,
        };

        save_index(&collection, path).map_err(|e| -> pyo3::PyErr { e.into() })
    }

    /// Load sequences from a previously saved index file.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to the serialized index file.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If deserialization fails or the k-mer length does not match.
    pub fn load(&mut self, path: &str) -> PyResult<()> {
        let collection = load_index(path).map_err(|e| -> pyo3::PyErr { e.into() })?;
        if collection.k != self.k {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Index k={} does not match current k={}",
                collection.k, self.k
            )));
        }
        for (name, serializable) in collection.sequences {
            let kmer_set: AHashSet<String> = collection
                .kmer_sets
                .get(&name)
                .map(|v| v.iter().cloned().collect())
                .unwrap_or_default();
            let seq_len = serializable.seq_bytes.len();
            let fm = rebuild_fm_from_bytes(&serializable.seq_bytes)
                .map_err(|e| -> pyo3::PyErr { e.into() })?;
            self.sequences.insert(
                name,
                SequenceData {
                    fm,
                    kmer_set,
                    seq_bytes: serializable.seq_bytes,
                    seq_len,
                },
            );
        }
        Ok(())
    }

    /// Get the k-mer length used for this index.
    ///
    /// Returns
    /// -------
    /// int
    ///     The k-mer length.
    #[getter]
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get the number of sequences in the index.
    ///
    /// Returns
    /// -------
    /// int
    ///     The number of indexed sequences.
    pub fn __len__(&self) -> usize {
        self.sequences.len()
    }

    /// String representation of the SequenceIndex.
    pub fn __repr__(&self) -> String {
        format!(
            "SequenceIndex(k={}, sequences={})",
            self.k,
            self.sequences.len()
        )
    }
}
