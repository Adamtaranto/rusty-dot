//! The main `SequenceIndex` class providing the core functionality.
//!
//! Provides rolling-hash k-mer index construction, k-mer lookup, sequence
//! comparison, and PAF output, with full Python bindings via PyO3.
// pyo3 pyfunction/pymethods return types trigger a false-positive useless_conversion lint.
#![allow(clippy::useless_conversion)]

use crate::error::RustyDotError;
use crate::kmer::build_kmer_set;
use crate::kmer_hash::{shared_fwd_pairs, shared_rev_pairs, KmerIndex};
use crate::merge::{merge_fwd_pairs, merge_rev_fwd_pairs, merge_rev_pairs, CoordPair};
use crate::paf::coords_to_paf;
use crate::serialize::{load_index, save_index, IndexCollection, SerializableSequence};
use crate::strand::STRAND_REV;
use pyo3::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::HashMap;

/// Stranded match coordinates: (query_start, query_end, target_start, target_end, strand).
type StrandedMatch = (usize, usize, usize, usize, String);

/// Squared d-genies gravity weight of a single match: `(1 + euclidean_length)^2`,
/// where the euclidean length spans both the query and target axes.  Large,
/// collinear blocks dominate, so a contig's best chromosome and sort position are
/// driven by its major alignments rather than short spurious hits.
fn match_weight(c: &CoordPair) -> f64 {
    let dx = c.target_end as f64 - c.target_start as f64;
    let dy = c.query_end as f64 - c.query_start as f64;
    let eucl = (dx * dx + dy * dy).sqrt();
    (1.0 + eucl).powi(2)
}

/// Order `items` along the concatenated `others` axis by gravity centre, following
/// the d-genies algorithm.
///
/// For each item we sum the squared match weights (`match_weight`) against every
/// "other" and assign the item to its best-matching other (argmax of the summed
/// weight).  Its sort position is the weighted mean of its match mid-points on the
/// concatenated other axis, using *only* the matches to that best other, normalised
/// by the total other length.  Matched items sort by ascending position (ties broken
/// by descending length then name); items with no matches sort last, by descending
/// length then name.
///
/// `self_is_query` selects the axis convention: the shared `match_map` is always
/// keyed `(query, target)`, so when ordering targets we swap the lookup key and read
/// the query-axis mid-points instead of the target-axis mid-points.
fn gravity_order(
    items: &[String],
    others: &[String],
    match_map: &HashMap<(String, String), Vec<CoordPair>>,
    seq_len: &HashMap<String, usize>,
    self_is_query: bool,
) -> Vec<String> {
    // Cumulative offsets so each match mid-point becomes an absolute coordinate
    // along the concatenated other axis (in the given `others` order).
    let mut other_offsets: HashMap<&str, f64> = HashMap::new();
    let mut acc = 0f64;
    for o in others {
        other_offsets.insert(o.as_str(), acc);
        acc += seq_len[o.as_str()] as f64;
    }
    let total_other = acc.max(1.0);

    let mid_of = |c: &CoordPair| -> f64 {
        if self_is_query {
            (c.target_start as f64 + c.target_end as f64) / 2.0
        } else {
            (c.query_start as f64 + c.query_end as f64) / 2.0
        }
    };
    let key_of = |it: &str, o: &str| -> (String, String) {
        if self_is_query {
            (it.to_string(), o.to_string())
        } else {
            (o.to_string(), it.to_string())
        }
    };

    let mut matched: Vec<(f64, usize, String)> = Vec::new();
    let mut unmatched: Vec<(usize, String)> = Vec::new();

    for it in items {
        // Total squared weight of this item against each other; pick the argmax.
        // The first other reaching the maximum wins (strict `>`), matching the
        // Python implementation's tie behaviour.
        let mut best_other: Option<&str> = None;
        let mut best_w = 0f64;
        for o in others {
            let ms = &match_map[&key_of(it, o)];
            if ms.is_empty() {
                continue;
            }
            let w: f64 = ms.iter().map(match_weight).sum();
            if w > best_w {
                best_w = w;
                best_other = Some(o.as_str());
            }
        }

        match best_other {
            None => unmatched.push((seq_len[it.as_str()], it.clone())),
            Some(bo) => {
                let ms = &match_map[&key_of(it, bo)];
                let offset = other_offsets[bo];
                let mut wsum = 0f64;
                let mut wpos = 0f64;
                for c in ms {
                    let w = match_weight(c);
                    wsum += w;
                    wpos += w * (offset + mid_of(c));
                }
                let pos = wpos / wsum / total_other;
                matched.push((pos, seq_len[it.as_str()], it.clone()));
            }
        }
    }

    // Matched: ascending gravity, then descending length, then name (determinism).
    matched.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.2.cmp(&b.2))
    });
    // Unmatched: descending length, then name.
    unmatched.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));

    matched
        .into_iter()
        .map(|(_, _, n)| n)
        .chain(unmatched.into_iter().map(|(_, n)| n))
        .collect()
}

/// Build the in-memory index data (rolling-hash k-mer index) for a single sequence.
///
/// This is the CPU-heavy per-sequence work (a single O(n) ntHash scan building
/// the forward and reverse-complement hash maps).  It is a free function taking
/// no `&self` so it can be called from a rayon parallel iterator while building
/// many sequences concurrently.
fn build_sequence_data(seq: &str, k: usize) -> Result<SequenceData, RustyDotError> {
    if k == 0 {
        return Err(RustyDotError::InvalidKmerLength(k));
    }
    let seq_bytes = seq.as_bytes().to_vec();
    let index = KmerIndex::build(&seq_bytes, k);
    let seq_len = seq.len();
    Ok(SequenceData {
        index,
        seq_bytes,
        seq_len,
    })
}

/// Build index data for many sequences, in parallel when the `parallel`
/// feature is enabled.
///
/// Returns the built `SequenceData` in the same order as `records`.  The caller
/// is responsible for releasing the GIL (via `Python::detach`) around
/// this call so the rayon worker threads can run without contention.
fn build_many_sequence_data(
    records: &[(String, String)],
    k: usize,
) -> Result<Vec<SequenceData>, RustyDotError> {
    #[cfg(feature = "parallel")]
    {
        records
            .par_iter()
            .map(|(_, seq)| build_sequence_data(seq, k))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        records
            .iter()
            .map(|(_, seq)| build_sequence_data(seq, k))
            .collect()
    }
}

/// In-memory store for a single sequence's index data.
struct SequenceData {
    /// Rolling-hash k-mer index (forward + reverse-complement hash maps).
    index: KmerIndex,
    /// Original sequence bytes (without sentinel), used for matching (byte
    /// verification) and serialization.
    seq_bytes: Vec<u8>,
    /// Original sequence length.
    seq_len: usize,
}

/// PyO3-exposed class for building and querying k-mer indexes for DNA sequences.
///
/// Each sequence added to the index receives its **own independent k-mer
/// index** — forward and reverse-complement ntHash maps built in a single O(n)
/// pass (see [`crate::kmer_hash`]).  Per-sequence indexes are independent, so
/// adding more sequences never modifies an existing one.
///
/// The index behaves as a **dictionary of per-sequence k-mer indexes**:
///
/// * `add_sequence` / `load_fasta` — **add** new entries to the collection;
///   calling either method multiple times accumulates sequences rather than
///   replacing them.
/// * If `add_sequence` (or `load_fasta`) is called with a name that already
///   exists in the index, the existing entry is **silently overwritten** with a
///   new index for the new sequence.
/// * Pairwise comparisons always operate on exactly two independent indexes.
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

/// Forward-strand comparison of two sequences, returning the coordinate pairs
/// that would be cached (merged runs, or sorted unmerged hits).
///
/// This is the cache-free, `&self`-free core of
/// [`SequenceIndex::compare_sequences`].  It takes only `&sequences` + `k` so it
/// can be invoked from a rayon parallel iterator without capturing the
/// `#[pyclass]`.  Callers must have already validated that the named sequences
/// exist.
fn compute_pair_fwd(
    sequences: &HashMap<String, SequenceData>,
    k: usize,
    query_name: &str,
    target_name: &str,
    merge: bool,
) -> Vec<CoordPair> {
    let query = &sequences[query_name];
    let target = &sequences[target_name];
    let hits = shared_fwd_pairs(
        &query.seq_bytes,
        &query.index,
        &target.seq_bytes,
        &target.index,
        k,
    );

    if hits.is_empty() {
        return Vec::new();
    }

    if merge {
        merge_fwd_pairs(hits, k)
    } else {
        use crate::strand::STRAND_FWD;
        let mut unmerged: Vec<CoordPair> = hits
            .into_iter()
            .map(|(qp, tp)| CoordPair {
                query_start: qp,
                query_end: qp + k,
                target_start: tp,
                target_end: tp + k,
                strand: STRAND_FWD,
            })
            .collect();
        unmerged.sort_unstable_by_key(|c| (c.query_start, c.target_start));
        unmerged
    }
}

/// Both-strand comparison of two sequences.
///
/// The cache-free, `&self`-free core of
/// [`SequenceIndex::compare_sequences_stranded`].
fn compute_pair_stranded(
    sequences: &HashMap<String, SequenceData>,
    k: usize,
    query_name: &str,
    target_name: &str,
    merge: bool,
) -> Vec<CoordPair> {
    let query = &sequences[query_name];
    let target = &sequences[target_name];

    let mut all_pairs: Vec<CoordPair> = Vec::new();

    // + strand hits.  The pair-based path skips the per-k-mer String maps
    // entirely — on repeat-rich assembly pairs those maps dominated the peak
    // memory of a comparison.
    let fwd_hits = shared_fwd_pairs(
        &query.seq_bytes,
        &query.index,
        &target.seq_bytes,
        &target.index,
        k,
    );
    if !fwd_hits.is_empty() {
        if merge {
            all_pairs.extend(merge_fwd_pairs(fwd_hits, k));
        } else {
            use crate::strand::STRAND_FWD;
            all_pairs.extend(fwd_hits.into_iter().map(|(qp, tp)| CoordPair {
                query_start: qp,
                query_end: qp + k,
                target_start: tp,
                target_end: tp + k,
                strand: STRAND_FWD,
            }));
        }
    }

    // - strand hits
    let rev_hits = shared_rev_pairs(
        &query.seq_bytes,
        &query.index,
        &target.seq_bytes,
        &target.index,
        k,
    );
    if !rev_hits.is_empty() {
        if merge {
            // Apply both anti-diagonal and co-diagonal merging for RC hits,
            // deduplicating identical blocks that arise when a single RC pair
            // has no neighbours on either diagonal.
            let anti = merge_rev_pairs(rev_hits.clone(), k);
            let co = merge_rev_fwd_pairs(rev_hits, k);
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
            all_pairs.extend(rev_hits.into_iter().map(|(qp, tp)| CoordPair {
                query_start: qp,
                query_end: qp + k,
                target_start: tp,
                target_end: tp + k,
                strand: STRAND_REV,
            }));
        }
    }

    all_pairs
}

/// Compute forward-strand coord pairs for a set of ordered pairs, in parallel
/// when the `parallel` feature is enabled.  Pure: reads only `sequences`.
fn compute_pairs_fwd_batch(
    sequences: &HashMap<String, SequenceData>,
    k: usize,
    pairs: &[(String, String)],
    merge: bool,
) -> Vec<Vec<CoordPair>> {
    #[cfg(feature = "parallel")]
    {
        pairs
            .par_iter()
            .map(|(q, t)| compute_pair_fwd(sequences, k, q, t, merge))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        pairs
            .iter()
            .map(|(q, t)| compute_pair_fwd(sequences, k, q, t, merge))
            .collect()
    }
}

/// Compute both-strand coord pairs for a set of ordered pairs, in parallel when
/// the `parallel` feature is enabled.  Pure: reads only `sequences`.
fn compute_pairs_stranded_batch(
    sequences: &HashMap<String, SequenceData>,
    k: usize,
    pairs: &[(String, String)],
    merge: bool,
) -> Vec<Vec<CoordPair>> {
    #[cfg(feature = "parallel")]
    {
        pairs
            .par_iter()
            .map(|(q, t)| compute_pair_stranded(sequences, k, q, t, merge))
            .collect()
    }
    #[cfg(not(feature = "parallel"))]
    {
        pairs
            .iter()
            .map(|(q, t)| compute_pair_stranded(sequences, k, q, t, merge))
            .collect()
    }
}

impl SequenceIndex {
    /// Enumerate every ordered `(query, target)` pair (`i != j`), compute the
    /// forward-strand coords for any not already cached — in parallel, with the
    /// GIL released — and populate `pair_cache`.  Returns all ordered pairs.
    fn ensure_fwd_pairs_cached(&mut self, py: Python<'_>, merge: bool) -> Vec<(String, String)> {
        let names: Vec<String> = self.sequences.keys().cloned().collect();
        let mut pairs: Vec<(String, String)> = Vec::with_capacity(names.len().saturating_mul(2));
        for i in 0..names.len() {
            for j in 0..names.len() {
                if i != j {
                    pairs.push((names[i].clone(), names[j].clone()));
                }
            }
        }

        let uncached: Vec<(String, String)> = pairs
            .iter()
            .filter(|(q, t)| !self.pair_cache.contains_key(&(q.clone(), t.clone(), merge)))
            .cloned()
            .collect();

        if !uncached.is_empty() {
            let k = self.k;
            let sequences = &self.sequences;
            // Release the GIL: the batch does no Python interaction.
            let results = py.detach(|| compute_pairs_fwd_batch(sequences, k, &uncached, merge));
            for ((q, t), coords) in uncached.into_iter().zip(results.into_iter()) {
                self.pair_cache.insert((q, t, merge), coords);
            }
        }
        pairs
    }
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
        let k = self.k;
        // Release the GIL for the CPU-heavy build (no Python interaction here).
        let data = py
            .detach(|| build_sequence_data(seq, k))
            .map_err(|e| -> pyo3::PyErr { e.into() })?;
        self.sequences.insert(name.to_string(), data);
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
        use needletail::parse_fastx_file;
        use std::collections::HashSet;
        use std::path::Path;

        let warnings = py.import("warnings")?;
        let mut seen_in_file: HashSet<String> = HashSet::new();

        // Pass 1: read every record into memory and validate names.  Parsing is
        // inherently sequential (single reader), so we do it up front and defer
        // the CPU-heavy index construction to a parallel pass.
        let mut records: Vec<(String, String)> = Vec::new();
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
            records.push((name, seq));
        }

        // Pass 2: build the FM-index + k-mer set for every record in parallel.
        // Release the GIL so the rayon worker threads run unhindered; this call
        // performs no Python interaction.
        let k = self.k;
        let built = py
            .detach(|| build_many_sequence_data(&records, k))
            .map_err(|e| -> pyo3::PyErr { e.into() })?;

        // Pass 3: merge the fully-validated batch into self.sequences (serial,
        // holds the GIL so the overwrite warnings can be emitted).
        let names: Vec<String> = records.into_iter().map(|(name, _)| name).collect();
        for (name, data) in names.iter().cloned().zip(built.into_iter()) {
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
            self.sequences.insert(name, data);
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
            // Regenerate the k-mer string set on demand from the retained bytes;
            // the matching path uses the rolling-hash index instead of storing
            // the (memory-heavy) string set.
            Some(data) => {
                let seq = String::from_utf8_lossy(&data.seq_bytes);
                let set = build_kmer_set(&seq, self.k).map_err(|e| -> pyo3::PyErr { e.into() })?;
                Ok(set.into_iter().collect())
            }
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

    /// Get the stored sequence for a named sequence.
    ///
    /// Returns the original sequence bases retained when the sequence was added,
    /// decoded as a UTF-8 string.  Used to write (optionally reordered and
    /// reoriented) FASTA output.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     The sequence name.
    ///
    /// Returns
    /// -------
    /// str
    ///     The sequence bases.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If the sequence name is not found.
    pub fn get_sequence(&self, name: &str) -> PyResult<String> {
        match self.sequences.get(name) {
            Some(data) => Ok(String::from_utf8_lossy(&data.seq_bytes).into_owned()),
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

        let coords = compute_pair_fwd(&self.sequences, self.k, query_name, target_name, merge);
        let tuples: Vec<(usize, usize, usize, usize)> = coords
            .iter()
            .map(|c| (c.query_start, c.query_end, c.target_start, c.target_end))
            .collect();
        self.pair_cache.insert(cache_key, coords);
        Ok(tuples)
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

        let all_pairs =
            compute_pair_stranded(&self.sequences, self.k, query_name, target_name, merge);
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

    /// Find shared k-mer matches for many (query, target) pairs in one call.
    ///
    /// Batched equivalent of `compare_sequences_stranded`: all pairs are
    /// computed inside a single native call with the GIL released, in parallel
    /// when the `parallel` feature is enabled.  This removes the per-pair
    /// Python-call overhead that dominates all-vs-all comparisons of
    /// many-contig assemblies (Q x T round-trips collapse into one).
    ///
    /// Parameters
    /// ----------
    /// pairs : list[tuple[str, str]]
    ///     Ordered ``(query_name, target_name)`` pairs to compare.  Names may
    ///     repeat; each entry produces one result list.
    /// merge : bool, optional
    ///     Whether to merge co-linear k-mer runs. Default is True.
    /// min_block_len : int, optional
    ///     Drop matches whose longest span (query or target) is shorter than
    ///     this many bases, before crossing back into Python.  Repeat-rich
    ///     genome pairs can produce millions of short blocks; filtering here
    ///     avoids materialising them as Python objects. Default is 0 (keep
    ///     all).
    ///
    /// Returns
    /// -------
    /// list[list[tuple[int, int, int, int, str]]]
    ///     One list per input pair, in input order, of
    ///     (query_start, query_end, target_start, target_end, strand) tuples.
    ///     Strand is ``"+"`` for forward and ``"-"`` for reverse-complement
    ///     matches.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If any sequence name is not found.
    #[pyo3(signature = (pairs, merge=true, min_block_len=0))]
    pub fn compare_pairs_stranded(
        &self,
        py: Python<'_>,
        pairs: Vec<(String, String)>,
        merge: bool,
        min_block_len: usize,
    ) -> PyResult<Vec<Vec<StrandedMatch>>> {
        for (q, t) in &pairs {
            for name in [q, t] {
                if !self.sequences.contains_key(name.as_str()) {
                    return Err(pyo3::exceptions::PyKeyError::new_err(format!(
                        "Sequence '{}' not found",
                        name
                    )));
                }
            }
        }

        let k = self.k;
        let sequences = &self.sequences;
        // Release the GIL: the batch does no Python interaction.
        let results = py.detach(|| compute_pairs_stranded_batch(sequences, k, &pairs, merge));
        Ok(results
            .into_iter()
            .map(|coords| {
                coords
                    .into_iter()
                    .filter(|c| {
                        min_block_len == 0
                            || (c.query_end - c.query_start).max(c.target_end - c.target_start)
                                >= min_block_len
                    })
                    .map(|c| {
                        (
                            c.query_start,
                            c.query_end,
                            c.target_start,
                            c.target_end,
                            (c.strand as char).to_string(),
                        )
                    })
                    .collect()
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
        py: Python<'_>,
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

        // Precompute the both-strand matches for every (query, target) pair once,
        // in parallel with the GIL released.  Both gravity passes below read from
        // this map, so each pair is compared a single time.
        let pair_list: Vec<(String, String)> = query_names
            .iter()
            .flat_map(|q| target_names.iter().map(move |t| (q.clone(), t.clone())))
            .collect();
        let match_map: HashMap<(String, String), Vec<CoordPair>> = {
            let k = self.k;
            let sequences = &self.sequences;
            let results =
                py.detach(|| compute_pairs_stranded_batch(sequences, k, &pair_list, true));
            pair_list.iter().cloned().zip(results.into_iter()).collect()
        };

        // Sequence lengths for every name involved (used for offsets and the
        // unmatched-contig length tie-break).
        let seq_len: HashMap<String, usize> = query_names
            .iter()
            .chain(target_names.iter())
            .map(|n| (n.clone(), self.sequences[n.as_str()].seq_len))
            .collect();

        // d-genies groups query contigs against the *displayed* target order, so we
        // order the targets first (using the queries in their input order for the
        // concatenated query axis) and then order the queries against the freshly
        // sorted target axis.  This asymmetry is intentional and deterministic.
        let sorted_t = gravity_order(&target_names, &query_names, &match_map, &seq_len, false);
        let sorted_q = gravity_order(&query_names, &sorted_t, &match_map, &seq_len, true);

        Ok((sorted_q, sorted_t))
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
    pub fn get_paf_all(&mut self, py: Python<'_>, merge: bool) -> PyResult<Vec<String>> {
        // Compute every pair in parallel (cache-filling), then format PAF lines
        // from the cached coords.
        let pairs = self.ensure_fwd_pairs_cached(py, merge);
        let mut all_paf: Vec<String> = Vec::new();
        for (q, t) in &pairs {
            let coords = self
                .pair_cache
                .get(&(q.clone(), t.clone(), merge))
                .expect("pair cache populated by ensure_fwd_pairs_cached");
            let query_len = self.sequences[q].seq_len;
            let target_len = self.sequences[t].seq_len;
            all_paf.extend(coords_to_paf(coords, q, query_len, t, target_len));
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
    pub fn precompute_all_pairs(
        &mut self,
        py: Python<'_>,
        merge: bool,
    ) -> PyResult<Vec<(String, String)>> {
        Ok(self.ensure_fwd_pairs_cached(py, merge))
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
            // Regenerate the k-mer string set from the retained bytes for the
            // serialized form (the in-memory index no longer stores it).
            let seq = String::from_utf8_lossy(&data.seq_bytes);
            let kmer_set = build_kmer_set(&seq, self.k).map_err(|e| -> pyo3::PyErr { e.into() })?;
            kmer_sets_map.insert(name.clone(), kmer_set.into_iter().collect::<Vec<_>>());
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
            let seq_len = serializable.seq_bytes.len();
            // Rebuild the rolling-hash index from the stored bytes.
            let index = KmerIndex::build(&serializable.seq_bytes, self.k);
            self.sequences.insert(
                name,
                SequenceData {
                    index,
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

#[cfg(test)]
mod gravity_tests {
    use super::*;

    /// Build a forward-strand `CoordPair` (strand is irrelevant to gravity).
    fn cp(qs: usize, qe: usize, ts: usize, te: usize) -> CoordPair {
        CoordPair {
            query_start: qs,
            query_end: qe,
            target_start: ts,
            target_end: te,
            strand: b'+',
        }
    }

    #[test]
    fn match_weight_is_one_plus_euclidean_squared() {
        // Target span 3, query span 4 -> euclidean length 5 -> (1 + 5)^2 = 36.
        let w = match_weight(&cp(0, 4, 0, 3));
        assert!((w - 36.0).abs() < 1e-9);
    }

    #[test]
    fn best_chromosome_argmax_ignores_small_hits() {
        // `q` has five tiny hits to `tA` but one big hit to `tB`; squared
        // weighting makes `tB` its best chromosome, so it sorts after `q2`
        // (which maps to the start of `tA`) on the concatenated [tA, tB] axis.
        let mut mm: HashMap<(String, String), Vec<CoordPair>> = HashMap::new();
        let small: Vec<CoordPair> = (0..5)
            .map(|i| cp(100 + i * 10, 110 + i * 10, 100 + i * 10, 110 + i * 10))
            .collect();
        mm.insert(("q".into(), "tA".into()), small);
        mm.insert(("q".into(), "tB".into()), vec![cp(300, 800, 400, 900)]);
        mm.insert(("q2".into(), "tA".into()), vec![cp(0, 400, 0, 400)]);
        mm.insert(("q2".into(), "tB".into()), vec![]);
        let seq_len: HashMap<String, usize> = [
            ("q".to_string(), 900usize),
            ("q2".to_string(), 500),
            ("tA".to_string(), 1000),
            ("tB".to_string(), 1000),
        ]
        .into_iter()
        .collect();
        let order = gravity_order(
            &["q".to_string(), "q2".to_string()],
            &["tA".to_string(), "tB".to_string()],
            &mm,
            &seq_len,
            true,
        );
        assert_eq!(order, vec!["q2".to_string(), "q".to_string()]);
    }

    #[test]
    fn unmatched_sorts_last_by_descending_length() {
        let mut mm: HashMap<(String, String), Vec<CoordPair>> = HashMap::new();
        mm.insert(("m".into(), "t".into()), vec![cp(0, 100, 0, 100)]);
        mm.insert(("u_long".into(), "t".into()), vec![]);
        mm.insert(("u_short".into(), "t".into()), vec![]);
        let seq_len: HashMap<String, usize> = [
            ("m".to_string(), 100usize),
            ("u_long".to_string(), 500),
            ("u_short".to_string(), 50),
            ("t".to_string(), 1000),
        ]
        .into_iter()
        .collect();
        let order = gravity_order(
            &["u_short".to_string(), "u_long".to_string(), "m".to_string()],
            &["t".to_string()],
            &mm,
            &seq_len,
            true,
        );
        assert_eq!(
            order,
            vec!["m".to_string(), "u_long".to_string(), "u_short".to_string()]
        );
    }
}
