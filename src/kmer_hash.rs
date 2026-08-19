//! Rolling-hash k-mer position index (ntHash) with a compact canonical layout.
//!
//! This module replaces the FM-index on the dot-plot *matching* path.  Each
//! sequence is scanned once with the [ntHash] rolling hash, which rolls both
//! the forward hash `fh` and reverse-complement hash `rh` of every k-mer
//! window in O(1) per position.  A key ntHash identity makes strand handling
//! allocation-free:
//!
//! ```text
//! rev_hash(Q) == fwd_hash(reverse_complement(Q))
//! ```
//!
//! # Memory layout
//!
//! Earlier revisions stored two hash maps per sequence (`fwd` and `rev`,
//! each `AHashMap<u64, SmallVec<u32>>`), costing ~100 bytes per base pair —
//! enough to abort the 2 GB WebAssembly heap while indexing real assemblies.
//! The index is now a single **canonical-hash CSR table**:
//!
//! * every window is keyed by `canon = min(fh, rh)`;
//! * a strand flag (the position's high bit) records *which* of the two
//!   hashes was the canonical one (`1` when `rh < fh`);
//! * entries are sorted by `(canon, position)` and compacted into three flat
//!   arrays (`hashes`, CSR `starts`, packed `entries`).
//!
//! This is ~12–16 bytes per base pair — both strands included — and lookups
//! intersect two indexes with a cache-friendly two-pointer walk instead of
//! hash probes.  Forward and reverse matches are both recovered from the one
//! table: windows sharing a canonical hash with **equal** flags have equal
//! forward hashes (same k-mer), while **opposite** flags mean one window is
//! the reverse complement of the other.
//!
//! Hashes are 64-bit, so collisions are astronomically rare but not
//! impossible.  Callers resolve them by byte-verifying one representative
//! k-mer per shared hash class against the retained sequence bytes, which
//! keeps matching output byte-identical to exact k-mer search.
//!
//! [ntHash]: https://github.com/bcgsc/ntHash

use ahash::AHashMap;

/// ntHash seed for each nucleotide (the canonical ntHash constants).
const H_A: u64 = 0x3c8b_fbb3_95c6_0474;
const H_C: u64 = 0x3193_c185_62a0_2b4c;
const H_G: u64 = 0x2032_3ed0_8257_2324;
const H_T: u64 = 0x2955_49f5_4be2_4456;

/// High bit of a packed entry: set when the canonical hash is the
/// reverse-complement hash (`rh < fh`).
const FLAG_BIT: u32 = 1 << 31;
/// Low 31 bits of a packed entry: the 0-based window start position.
/// Caps indexable sequences at 2^31 - 1 bases (~2.1 Gbp), far above any
/// single contig this library targets.
const POS_MASK: u32 = FLAG_BIT - 1;

/// Forward-strand seed for a base (`0` for anything non-ACGT).
#[inline]
fn h(b: u8) -> u64 {
    match b {
        b'A' => H_A,
        b'C' => H_C,
        b'G' => H_G,
        b'T' => H_T,
        _ => 0,
    }
}

/// Reverse-complement seed for a base: the forward seed of its complement.
#[inline]
fn rc(b: u8) -> u64 {
    match b {
        b'A' => H_T,
        b'C' => H_G,
        b'G' => H_C,
        b'T' => H_A,
        _ => 0,
    }
}

/// True for the four canonical uppercase DNA bases.
#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Compute the forward ntHash of a k-mer directly (non-rolling).
///
/// Used for unit testing the rolling implementation and for callers that need a
/// one-off hash.  The k-mer must be exactly `k` bytes.
pub fn forward_hash(kmer: &[u8]) -> u64 {
    let k = kmer.len();
    let mut fh = 0u64;
    for (i, &b) in kmer.iter().enumerate() {
        fh ^= h(b).rotate_left((k - 1 - i) as u32);
    }
    fh
}

/// Compute the reverse-complement ntHash of a k-mer directly (non-rolling).
///
/// Satisfies `reverse_hash(q) == forward_hash(revcomp(q))`.
pub fn reverse_hash(kmer: &[u8]) -> u64 {
    let mut rh = 0u64;
    for (i, &b) in kmer.iter().enumerate() {
        rh ^= rc(b).rotate_left(i as u32);
    }
    rh
}

/// A compact canonical-hash k-mer index for a single sequence.
///
/// Three flat arrays in CSR form: `hashes[i]` is the i-th distinct canonical
/// hash (ascending); its packed entries are
/// `entries[starts[i] as usize..starts[i + 1] as usize]`, each holding a
/// window start position in the low 31 bits and the strand flag in the high
/// bit.  Entries within a group are position-ascending.
#[derive(Default)]
pub struct KmerIndex {
    /// Distinct canonical hashes, ascending.
    hashes: Vec<u64>,
    /// CSR offsets into `entries` (length `hashes.len() + 1`).
    starts: Vec<u32>,
    /// Packed `(strand_flag << 31) | position` entries.
    entries: Vec<u32>,
}

impl KmerIndex {
    /// Build the canonical hash index for `seq` at k-mer length `k`.
    ///
    /// Only windows consisting entirely of the canonical bases A, C, G, T are
    /// indexed — windows containing `N` (or any other byte) are skipped, exactly
    /// matching the set of k-mers considered by the exact-search path.  Within
    /// each maximal ACGT run of length ≥ `k`, the first window's hashes are
    /// computed directly and every subsequent window is rolled in O(1).
    ///
    /// # Arguments
    ///
    /// * `seq` - The sequence bytes (uppercase DNA, no sentinel).
    /// * `k` - The k-mer length.
    ///
    /// # Returns
    ///
    /// A populated [`KmerIndex`]; empty when `k == 0` or `seq` is shorter than
    /// `k`.
    pub fn build(seq: &[u8], k: usize) -> KmerIndex {
        let n = seq.len();
        if k == 0 || n < k {
            return KmerIndex::default();
        }

        // Pass 1: collect (canonical hash, packed position) pairs.
        let mut pairs: Vec<(u64, u32)> = Vec::with_capacity(n.saturating_sub(k) + 1);
        let mut push = |fh: u64, rh: u64, pos: usize| {
            let (canon, flag) = if rh < fh { (rh, FLAG_BIT) } else { (fh, 0) };
            pairs.push((canon, pos as u32 | flag));
        };

        let mut start = 0usize;
        while start < n {
            // Skip non-ACGT bases to find the beginning of a maximal ACGT run.
            if !is_acgt(seq[start]) {
                start += 1;
                continue;
            }
            let mut end = start;
            while end < n && is_acgt(seq[end]) {
                end += 1;
            }

            // The run [start, end) yields windows [pos, pos+k) for
            // pos in start..=end-k, provided the run is at least k long.
            if end - start >= k {
                // First window in the run: compute both hashes directly.
                let mut fh = 0u64;
                let mut rh = 0u64;
                for j in 0..k {
                    fh ^= h(seq[start + j]).rotate_left((k - 1 - j) as u32);
                    rh ^= rc(seq[start + j]).rotate_left(j as u32);
                }
                push(fh, rh, start);

                // Roll across the remaining windows in this run.
                let last_pos = end - k;
                for pos in (start + 1)..=last_pos {
                    let out = seq[pos - 1];
                    let inb = seq[pos + k - 1];
                    fh = fh.rotate_left(1) ^ h(out).rotate_left(k as u32) ^ h(inb);
                    rh = rh.rotate_right(1)
                        ^ rc(out).rotate_right(1)
                        ^ rc(inb).rotate_left((k - 1) as u32);
                    push(fh, rh, pos);
                }
            }

            start = end;
        }

        // Pass 2: sort by (hash, position) — the flag lives above the
        // position bits, so sort on the masked position to keep groups
        // position-ascending regardless of strand.
        pairs.sort_unstable_by_key(|&(hash, packed)| (hash, packed & POS_MASK));

        // Pass 3: compact into CSR form.
        let mut idx = KmerIndex {
            hashes: Vec::new(),
            starts: Vec::new(),
            entries: Vec::with_capacity(pairs.len()),
        };
        for (hash, packed) in pairs {
            if idx.hashes.last() != Some(&hash) {
                idx.hashes.push(hash);
                idx.starts.push(idx.entries.len() as u32);
            }
            idx.entries.push(packed);
        }
        idx.starts.push(idx.entries.len() as u32);
        idx.hashes.shrink_to_fit();
        idx.entries.shrink_to_fit();
        idx
    }

    /// Whether the index contains no windows.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Iterate `(canonical_hash, packed_entries)` groups in ascending hash
    /// order.  Entry packing: high bit = strand flag, low 31 bits = position.
    pub fn iter_groups(&self) -> impl Iterator<Item = (u64, &[u32])> {
        self.hashes.iter().enumerate().map(move |(i, &hash)| {
            let lo = self.starts[i] as usize;
            let hi = self.starts[i + 1] as usize;
            (hash, &self.entries[lo..hi])
        })
    }
}

/// Unpack an entry's position.
#[inline]
fn pos_of(packed: u32) -> usize {
    (packed & POS_MASK) as usize
}

/// Whether an entry's canonical hash was the reverse-complement hash.
#[inline]
fn is_rc(packed: u32) -> bool {
    packed & FLAG_BIT != 0
}

/// Split a group's packed entries into (flag-0 positions, flag-1 positions),
/// preserving the ascending position order within each class.
fn split_by_flag(entries: &[u32]) -> (Vec<usize>, Vec<usize>) {
    let mut plain = Vec::new();
    let mut flagged = Vec::new();
    for &e in entries {
        if is_rc(e) {
            flagged.push(pos_of(e));
        } else {
            plain.push(pos_of(e));
        }
    }
    (plain, flagged)
}

/// Walk two indexes' ascending hash arrays and call `on_shared` for every
/// canonical hash present in both, with each side's packed entry group.
fn for_shared_groups(
    q_index: &KmerIndex,
    t_index: &KmerIndex,
    mut on_shared: impl FnMut(&[u32], &[u32]),
) {
    let mut qi = 0usize;
    let mut ti = 0usize;
    while qi < q_index.hashes.len() && ti < t_index.hashes.len() {
        let qh = q_index.hashes[qi];
        let th = t_index.hashes[ti];
        match qh.cmp(&th) {
            std::cmp::Ordering::Less => qi += 1,
            std::cmp::Ordering::Greater => ti += 1,
            std::cmp::Ordering::Equal => {
                let q_lo = q_index.starts[qi] as usize;
                let q_hi = q_index.starts[qi + 1] as usize;
                let t_lo = t_index.starts[ti] as usize;
                let t_hi = t_index.starts[ti + 1] as usize;
                on_shared(&q_index.entries[q_lo..q_hi], &t_index.entries[t_lo..t_hi]);
                qi += 1;
                ti += 1;
            }
        }
    }
}

/// Build the shared forward-strand coordinate maps for two sequences.
///
/// For every canonical hash present in both indexes, entries whose strand
/// flags are **equal** have equal forward hashes — i.e. the same k-mer occurs
/// in both sequences (subject to byte verification of one representative,
/// rejecting the astronomically rare 64-bit collision).  A single canonical
/// group can therefore contribute up to two k-mers (a k-mer and its reverse
/// complement, one per flag class).  The returned maps are keyed by the k-mer
/// string and are drop-in inputs for the `merge_*` functions — reproducing
/// exactly what exact k-mer search produced.
///
/// # Arguments
///
/// * `q_seq` / `t_seq` - Query and target sequence bytes (for byte verification).
/// * `q_index` / `t_index` - The corresponding [`KmerIndex`]es.
/// * `k` - The k-mer length.
///
/// # Returns
///
/// `(query_coords, target_coords)`, each mapping a shared k-mer string to its
/// ascending positions in the respective sequence.
pub fn shared_fwd_coords(
    q_seq: &[u8],
    q_index: &KmerIndex,
    t_seq: &[u8],
    t_index: &KmerIndex,
    k: usize,
) -> (AHashMap<String, Vec<usize>>, AHashMap<String, Vec<usize>>) {
    let mut query_coords: AHashMap<String, Vec<usize>> = AHashMap::new();
    let mut target_coords: AHashMap<String, Vec<usize>> = AHashMap::new();

    for_shared_groups(q_index, t_index, |q_entries, t_entries| {
        let (q_plain, q_flag) = split_by_flag(q_entries);
        let (t_plain, t_flag) = split_by_flag(t_entries);
        for (q_pos, t_pos) in [(q_plain, t_plain), (q_flag, t_flag)] {
            if q_pos.is_empty() || t_pos.is_empty() {
                continue;
            }
            // Verify one representative k-mer to reject hash collisions.
            let q_rep = &q_seq[q_pos[0]..q_pos[0] + k];
            let t_rep = &t_seq[t_pos[0]..t_pos[0] + k];
            if q_rep != t_rep {
                continue;
            }
            let kmer = String::from_utf8_lossy(q_rep).into_owned();
            query_coords.insert(kmer.clone(), q_pos);
            target_coords.insert(kmer, t_pos);
        }
    });

    (query_coords, target_coords)
}

/// Build the shared reverse-strand coordinate maps for two sequences.
///
/// Finds query k-mers whose reverse complement occurs in the target.  Under
/// the canonical layout these are entries sharing a canonical hash with
/// **opposite** strand flags (using the ntHash identity
/// `rev_hash(Q) == fwd_hash(RC(Q))`), plus the palindromic case
/// (`Q == RC(Q)`) where equal-flag forward matches are also reverse matches.
/// Returns `(target_rev_coords, query_rev_coords)` keyed by the *original*
/// query k-mer, matching the inputs expected by `merge_rev_*`.
pub fn shared_rev_coords(
    q_seq: &[u8],
    q_index: &KmerIndex,
    t_seq: &[u8],
    t_index: &KmerIndex,
    k: usize,
) -> (AHashMap<String, Vec<usize>>, AHashMap<String, Vec<usize>>) {
    let mut query_rev: AHashMap<String, Vec<usize>> = AHashMap::new();
    let mut target_rev: AHashMap<String, Vec<usize>> = AHashMap::new();

    for_shared_groups(q_index, t_index, |q_entries, t_entries| {
        let (q_plain, q_flag) = split_by_flag(q_entries);
        let (t_plain, t_flag) = split_by_flag(t_entries);
        // Opposite flag classes: the target windows are the reverse
        // complement of the query windows.
        for (q_pos, t_pos) in [(&q_plain, &t_flag), (&q_flag, &t_plain)] {
            if q_pos.is_empty() || t_pos.is_empty() {
                continue;
            }
            let q_rep = &q_seq[q_pos[0]..q_pos[0] + k];
            let t_rep = &t_seq[t_pos[0]..t_pos[0] + k];
            if crate::strand::revcomp(q_rep) != t_rep {
                continue;
            }
            let kmer = String::from_utf8_lossy(q_rep).into_owned();
            query_rev.insert(kmer.clone(), q_pos.clone());
            target_rev.insert(kmer, t_pos.clone());
        }
        // Palindromic k-mers (Q == RC(Q)) hash identically on both strands
        // (fh == rh, flag 0 on both sides): a forward co-occurrence is then
        // also a reverse-complement occurrence.
        if !q_plain.is_empty() && !t_plain.is_empty() {
            let q_rep = &q_seq[q_plain[0]..q_plain[0] + k];
            let t_rep = &t_seq[t_plain[0]..t_plain[0] + k];
            if q_rep == t_rep && crate::strand::revcomp(q_rep) == q_rep {
                let kmer = String::from_utf8_lossy(q_rep).into_owned();
                query_rev.insert(kmer.clone(), q_plain.clone());
                target_rev.insert(kmer, t_plain.clone());
            }
        }
    });

    (target_rev, query_rev)
}

/// Emit raw forward-strand hit pairs `(query_pos, target_pos)` for two
/// sequences, without materialising per-k-mer maps.
///
/// Semantically identical to cross-joining the maps returned by
/// [`shared_fwd_coords`], but the intermediate String-keyed maps are skipped
/// entirely — on repeat-rich assembly pairs those maps were the peak-memory
/// hot spot of a comparison.  Feed the result to
/// [`crate::merge::merge_fwd_pairs`].
pub fn shared_fwd_pairs(
    q_seq: &[u8],
    q_index: &KmerIndex,
    t_seq: &[u8],
    t_index: &KmerIndex,
    k: usize,
) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for_shared_groups(q_index, t_index, |q_entries, t_entries| {
        emit_fwd_group(q_seq, t_seq, k, q_entries, t_entries, &mut pairs);
    });
    pairs
}

/// Emit the forward-strand hit pairs for one shared canonical group.
///
/// Iterates the two flag classes in place (no per-group allocation — this
/// path runs once per shared hash and dominates comparison time on large
/// sequences).
#[inline]
fn emit_fwd_group(
    q_seq: &[u8],
    t_seq: &[u8],
    k: usize,
    q_entries: &[u32],
    t_entries: &[u32],
    pairs: &mut Vec<(usize, usize)>,
) {
    for flag in [0u32, FLAG_BIT] {
        let Some(&q_first) = q_entries.iter().find(|&&e| e & FLAG_BIT == flag) else {
            continue;
        };
        let Some(&t_first) = t_entries.iter().find(|&&e| e & FLAG_BIT == flag) else {
            continue;
        };
        // Verify one representative k-mer to reject hash collisions.
        let q0 = pos_of(q_first);
        let t0 = pos_of(t_first);
        if q_seq[q0..q0 + k] != t_seq[t0..t0 + k] {
            continue;
        }
        for &qe in q_entries.iter().filter(|&&e| e & FLAG_BIT == flag) {
            for &te in t_entries.iter().filter(|&&e| e & FLAG_BIT == flag) {
                pairs.push((pos_of(qe), pos_of(te)));
            }
        }
    }
}

/// Raw `(query_pos, target_pos)` hit pairs for one strand.
pub type HitPairs = Vec<(usize, usize)>;

/// Emit both strands' hit pairs in a **single** walk over the two indexes.
///
/// Equivalent to calling [`shared_fwd_pairs`] and [`shared_rev_pairs`] but
/// intersecting the sorted hash arrays only once — the walk itself is a
/// significant share of comparison time on large sequences.
pub fn shared_stranded_pairs(
    q_seq: &[u8],
    q_index: &KmerIndex,
    t_seq: &[u8],
    t_index: &KmerIndex,
    k: usize,
) -> (HitPairs, HitPairs) {
    let mut fwd: HitPairs = Vec::new();
    let mut rev: HitPairs = Vec::new();
    for_shared_groups(q_index, t_index, |q_entries, t_entries| {
        emit_fwd_group(q_seq, t_seq, k, q_entries, t_entries, &mut fwd);
        emit_rev_group(q_seq, t_seq, k, q_entries, t_entries, &mut rev);
    });
    (fwd, rev)
}

/// Emit raw reverse-strand hit pairs `(query_pos, rc_target_pos)` for two
/// sequences, without materialising per-k-mer maps.
///
/// Semantically identical to cross-joining the maps returned by
/// [`shared_rev_coords`] (including the palindromic case).  Feed the result
/// to [`crate::merge::merge_rev_pairs`] /
/// [`crate::merge::merge_rev_fwd_pairs`].
pub fn shared_rev_pairs(
    q_seq: &[u8],
    q_index: &KmerIndex,
    t_seq: &[u8],
    t_index: &KmerIndex,
    k: usize,
) -> Vec<(usize, usize)> {
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for_shared_groups(q_index, t_index, |q_entries, t_entries| {
        emit_rev_group(q_seq, t_seq, k, q_entries, t_entries, &mut pairs);
    });
    pairs
}

/// Emit the reverse-strand hit pairs for one shared canonical group.
#[inline]
fn emit_rev_group(
    q_seq: &[u8],
    t_seq: &[u8],
    k: usize,
    q_entries: &[u32],
    t_entries: &[u32],
    pairs: &mut Vec<(usize, usize)>,
) {
    {
        // Opposite flag classes (in place, no per-group allocation): the
        // target windows are the reverse complement of the query windows.
        for (q_flag, t_flag) in [(0u32, FLAG_BIT), (FLAG_BIT, 0u32)] {
            let Some(&q_first) = q_entries.iter().find(|&&e| e & FLAG_BIT == q_flag) else {
                continue;
            };
            let Some(&t_first) = t_entries.iter().find(|&&e| e & FLAG_BIT == t_flag) else {
                continue;
            };
            let q0 = pos_of(q_first);
            let t0 = pos_of(t_first);
            if crate::strand::revcomp(&q_seq[q0..q0 + k]) != t_seq[t0..t0 + k] {
                continue;
            }
            for &qe in q_entries.iter().filter(|&&e| e & FLAG_BIT == q_flag) {
                for &te in t_entries.iter().filter(|&&e| e & FLAG_BIT == t_flag) {
                    pairs.push((pos_of(qe), pos_of(te)));
                }
            }
        }
        // Palindromic k-mers: a forward co-occurrence is also an RC match.
        let Some(&q_first) = q_entries.iter().find(|&&e| e & FLAG_BIT == 0) else {
            return;
        };
        let Some(&t_first) = t_entries.iter().find(|&&e| e & FLAG_BIT == 0) else {
            return;
        };
        let q0 = pos_of(q_first);
        let t0 = pos_of(t_first);
        let q_rep = &q_seq[q0..q0 + k];
        if q_rep == &t_seq[t0..t0 + k] && crate::strand::revcomp(q_rep) == q_rep {
            for &qe in q_entries.iter().filter(|&&e| e & FLAG_BIT == 0) {
                for &te in t_entries.iter().filter(|&&e| e & FLAG_BIT == 0) {
                    pairs.push((pos_of(qe), pos_of(te)));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strand::revcomp;

    /// All indexed positions of an index, ascending.
    fn all_positions(idx: &KmerIndex) -> Vec<u32> {
        let mut v: Vec<u32> = idx
            .iter_groups()
            .flat_map(|(_, entries)| entries.iter().map(|&e| e & POS_MASK))
            .collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn reverse_hash_equals_forward_hash_of_revcomp() {
        let kmers: &[&[u8]] = &[b"ACGT", b"AAAA", b"ACGTGCA", b"TTTTAAAAC", b"GATTACA"];
        for kmer in kmers {
            let rc_bytes = revcomp(kmer);
            assert_eq!(
                reverse_hash(kmer),
                forward_hash(&rc_bytes),
                "reverse_hash({:?}) must equal forward_hash(revcomp)",
                std::str::from_utf8(kmer).unwrap()
            );
        }
    }

    #[test]
    fn rolling_canonical_matches_direct() {
        let seq = b"ACGTACGTTGCAACGTNNACGTACGT";
        for k in [4usize, 5] {
            let idx = KmerIndex::build(seq, k);
            for (hash, entries) in idx.iter_groups() {
                for &e in entries {
                    let pos = (e & POS_MASK) as usize;
                    let window = &seq[pos..pos + k];
                    let fh = forward_hash(window);
                    let rh = reverse_hash(window);
                    assert_eq!(hash, fh.min(rh), "canonical hash mismatch");
                    assert_eq!(e & FLAG_BIT != 0, rh < fh, "strand flag mismatch");
                }
            }
        }
    }

    #[test]
    fn groups_sorted_and_positions_ascending() {
        let seq: Vec<u8> = b"ACGTTGCA".iter().cloned().cycle().take(400).collect();
        let idx = KmerIndex::build(&seq, 6);
        let hashes: Vec<u64> = idx.iter_groups().map(|(h, _)| h).collect();
        let mut sorted = hashes.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(hashes, sorted, "groups must be ascending and distinct");
        for (_, entries) in idx.iter_groups() {
            let pos: Vec<u32> = entries.iter().map(|&e| e & POS_MASK).collect();
            let mut ps = pos.clone();
            ps.sort_unstable();
            assert_eq!(pos, ps, "positions within a group must ascend");
        }
    }

    #[test]
    fn skips_windows_containing_n() {
        // Only the two flanking ACGT runs of length 4 should be indexed.
        let seq = b"ACGTNACGT";
        let idx = KmerIndex::build(seq, 4);
        assert_eq!(all_positions(&idx), vec![0, 5]);
    }

    #[test]
    fn all_acgt_windows_indexed() {
        let seq = b"ACGTACGT";
        let idx = KmerIndex::build(seq, 4);
        // n - k + 1 = 5 windows at positions 0..=4.
        assert_eq!(all_positions(&idx), vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn repeated_kmer_groups_positions() {
        let seq = b"AAAAAA";
        let idx = KmerIndex::build(seq, 3);
        // "AAA" (canonical with "TTT") occurs at positions 0..=3, one group.
        assert_eq!(idx.iter_groups().count(), 1);
        assert_eq!(all_positions(&idx), vec![0, 1, 2, 3]);
    }

    #[test]
    fn empty_when_seq_shorter_than_k() {
        let idx = KmerIndex::build(b"ACG", 4);
        assert!(idx.is_empty());
    }

    #[test]
    fn k_equal_one_and_large_rotations() {
        // k=1 exercises the rotate_left((k-1)=0) and rotate_left(k=1) paths.
        let idx = KmerIndex::build(b"ACGT", 1);
        assert_eq!(all_positions(&idx), vec![0, 1, 2, 3]);
    }

    #[test]
    fn large_k_rotation_consistency() {
        // k > 64 exercises rotate amounts that wrap the 64-bit rotation.
        let seq: Vec<u8> = b"ACGT".iter().cloned().cycle().take(200).collect();
        let k = 121;
        let idx = KmerIndex::build(&seq, k);
        for (hash, entries) in idx.iter_groups() {
            for &e in entries {
                let pos = (e & POS_MASK) as usize;
                let window = &seq[pos..pos + k];
                assert_eq!(hash, forward_hash(window).min(reverse_hash(window)));
            }
        }
    }

    #[test]
    fn fwd_coords_finds_shared_kmers_both_flag_classes() {
        // GATTACA's canon may come from either strand; shared occurrences on
        // the same strand class must be reported as forward matches.
        let q = b"GATTACAGGGG";
        let t = b"CCCGATTACA";
        let k = 7;
        let (qc, tc) = shared_fwd_coords(q, &KmerIndex::build(q, k), t, &KmerIndex::build(t, k), k);
        assert_eq!(qc.get("GATTACA"), Some(&vec![0usize]));
        assert_eq!(tc.get("GATTACA"), Some(&vec![3usize]));
    }

    #[test]
    fn rev_coords_finds_reverse_complement_matches() {
        // TGTAATC = revcomp(GATTACA): a reverse-strand match, no forward one.
        let q = b"GATTACA";
        let t = b"TGTAATC";
        let k = 7;
        let (tr, qr) = shared_rev_coords(q, &KmerIndex::build(q, k), t, &KmerIndex::build(t, k), k);
        assert_eq!(qr.get("GATTACA"), Some(&vec![0usize]));
        assert_eq!(tr.get("GATTACA"), Some(&vec![0usize]));
        let (qc, _) = shared_fwd_coords(q, &KmerIndex::build(q, k), t, &KmerIndex::build(t, k), k);
        assert!(qc.is_empty(), "no forward match between RC-only sequences");
    }

    #[test]
    fn palindromic_kmer_matches_both_strands() {
        // ACGT == revcomp(ACGT): a shared palindrome is both a forward and a
        // reverse match (mirroring the old two-map behaviour).
        let q = b"AACGTA";
        let t = b"TACGTT";
        let k = 4;
        let (qc, _) = shared_fwd_coords(q, &KmerIndex::build(q, k), t, &KmerIndex::build(t, k), k);
        assert!(qc.contains_key("ACGT"), "palindrome forward match expected");
        let (tr, qr) = shared_rev_coords(q, &KmerIndex::build(q, k), t, &KmerIndex::build(t, k), k);
        assert!(qr.contains_key("ACGT"), "palindrome reverse match expected");
        assert_eq!(tr.get("ACGT"), Some(&vec![1usize]));
    }
}
