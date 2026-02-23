//! The main `SequenceIndex` class providing the core functionality.
//!
//! Provides FM-index construction, k-mer lookup, sequence comparison,
//! and PAF output, with full Python bindings via PyO3.

use crate::fasta::read_fasta;
use crate::kmer::{
    build_kmer_set, find_kmer_coords_in_index, find_rev_coords_in_index, sequence_to_index_text,
    FmIdx,
};
use crate::merge::{merge_fwd_runs, merge_kmer_runs, merge_rev_runs, CoordPair};
use crate::paf::coords_to_paf;
use crate::serialize::{
    load_index, rebuild_fm_from_bytes, save_index, IndexCollection, SerializableSequence,
};
use crate::strand::{revcomp, STRAND_REV};
use ahash::AHashSet;
use pyo3::prelude::*;
use std::collections::HashMap;

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
/// Consumes FASTA file paths or individual sequences, builds FM-indexes
/// and k-mer sets, and supports pairwise k-mer coordinate lookup with
/// optional merging of sequential runs.
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
    /// Parameters
    /// ----------
    /// name : str
    ///     The sequence identifier.
    /// seq : str
    ///     The DNA sequence string (uppercase recommended).
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the FM-index cannot be built.
    pub fn add_sequence(&mut self, name: &str, seq: &str) -> PyResult<()> {
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
    /// Parameters
    /// ----------
    /// path : str
    ///     Path to the FASTA or FASTA.gz file.
    ///
    /// Returns
    /// -------
    /// list[str]
    ///     List of sequence names that were added.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the file cannot be read or parsed.
    pub fn load_fasta(&mut self, path: &str) -> PyResult<Vec<String>> {
        let seqs = read_fasta(path).map_err(|e| -> pyo3::PyErr { e.into() })?;
        let mut names = Vec::new();
        for (name, seq) in &seqs {
            let text = sequence_to_index_text(seq);
            let fm = FmIdx::new(text).map_err(|e| -> pyo3::PyErr { e.into() })?;
            let kmer_set =
                build_kmer_set(seq, self.k).map_err(|e| -> pyo3::PyErr { e.into() })?;
            let seq_len = seq.len();
            let seq_bytes = seq.as_bytes().to_vec();
            self.sequences.insert(
                name.clone(),
                SequenceData {
                    fm,
                    kmer_set,
                    seq_bytes,
                    seq_len,
                },
            );
            names.push(name.clone());
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
    ) -> PyResult<Vec<(usize, usize, usize, usize, String)>> {
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
                all_pairs.extend(merge_rev_runs(&target_rev, &query_rev, self.k));
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
                let matches =
                    self.compare_sequences_stranded(q.as_str(), t.as_str(), true)?;
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
        q_gravity.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

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
                let matches =
                    self.compare_sequences_stranded(q.as_str(), t.as_str(), true)?;
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
        t_gravity.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        Ok((
            q_gravity.into_iter().map(|(_, n)| n).collect(),
            t_gravity.into_iter().map(|(_, n)| n).collect(),
        ))
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
