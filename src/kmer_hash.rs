//! Rolling-hash k-mer position index (ntHash).
//!
//! This module replaces the FM-index on the dot-plot *matching* path.  Instead
//! of building a suffix array per sequence and searching it for every observed
//! k-mer, each sequence is scanned once with the [ntHash] rolling hash to build
//! two maps:
//!
//! * `fwd` — forward-strand hash → sorted 0-based positions of that k-mer, and
//! * `rev` — reverse-complement hash → sorted positions.
//!
//! ntHash rolls both the forward and reverse-complement hash of a k-mer in O(1)
//! per position, so building the whole index is O(n) with a tiny constant — far
//! cheaper than suffix-array construction, especially for large single
//! sequences.  A key ntHash identity makes strand handling allocation-free:
//!
//! ```text
//! rev_hash(Q) == fwd_hash(reverse_complement(Q))
//! ```
//!
//! so reverse-strand matches between a query `Q` and a target `T` (i.e. the
//! reverse complement of `Q` occurs in `T`) are found by intersecting the
//! query's `rev` map with the target's `fwd` map — no per-k-mer reverse
//! complement is computed.
//!
//! Hashes are 64-bit, so collisions are astronomically rare but not impossible.
//! Callers resolve them by byte-verifying one representative k-mer per shared
//! hash against the retained sequence bytes, which keeps matching output
//! byte-identical to exact k-mer search.
//!
//! [ntHash]: https://github.com/bcgsc/ntHash

use ahash::AHashMap;
use smallvec::SmallVec;

/// ntHash seed for each nucleotide (the canonical ntHash constants).
const H_A: u64 = 0x3c8b_fbb3_95c6_0474;
const H_C: u64 = 0x3193_c185_62a0_2b4c;
const H_G: u64 = 0x2032_3ed0_8257_2324;
const H_T: u64 = 0x2955_49f5_4be2_4456;

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

/// Per-k-mer position lists.  Most k-mers occur once, so a `SmallVec` inline
/// capacity of one avoids a heap allocation in the common case.
pub type PosVec = SmallVec<[u32; 1]>;

/// A rolling-hash k-mer index for a single sequence.
///
/// `fwd` maps the forward ntHash of each ACGT-only k-mer window to its sorted
/// 0-based start positions; `rev` maps the reverse-complement ntHash likewise.
#[derive(Default)]
pub struct KmerIndex {
    /// Forward-strand hash → ascending positions.
    pub fwd: AHashMap<u64, PosVec>,
    /// Reverse-complement hash → ascending positions.
    pub rev: AHashMap<u64, PosVec>,
}

impl KmerIndex {
    /// Build the forward/reverse hash index for `seq` at k-mer length `k`.
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
        let mut idx = KmerIndex::default();
        let n = seq.len();
        if k == 0 || n < k {
            return idx;
        }

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
                idx.fwd.entry(fh).or_default().push(start as u32);
                idx.rev.entry(rh).or_default().push(start as u32);

                // Roll across the remaining windows in this run.
                let last_pos = end - k;
                for pos in (start + 1)..=last_pos {
                    let out = seq[pos - 1];
                    let inb = seq[pos + k - 1];
                    fh = fh.rotate_left(1) ^ h(out).rotate_left(k as u32) ^ h(inb);
                    rh = rh.rotate_right(1)
                        ^ rc(out).rotate_right(1)
                        ^ rc(inb).rotate_left((k - 1) as u32);
                    idx.fwd.entry(fh).or_default().push(pos as u32);
                    idx.rev.entry(rh).or_default().push(pos as u32);
                }
            }

            start = end;
        }
        idx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strand::revcomp;

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
    fn rolling_forward_matches_direct() {
        let seq = b"ACGTACGTTGCAACGTNNACGTACGT";
        let k = 4;
        let idx = KmerIndex::build(seq, k);
        // Every indexed position's rolled forward hash must equal the direct hash.
        for (hash, positions) in &idx.fwd {
            for &pos in positions {
                let window = &seq[pos as usize..pos as usize + k];
                assert_eq!(forward_hash(window), *hash);
            }
        }
    }

    #[test]
    fn rolling_reverse_matches_direct() {
        let seq = b"ACGTACGTTGCAACGTNNACGTACGT";
        let k = 5;
        let idx = KmerIndex::build(seq, k);
        for (hash, positions) in &idx.rev {
            for &pos in positions {
                let window = &seq[pos as usize..pos as usize + k];
                assert_eq!(reverse_hash(window), *hash);
            }
        }
    }

    #[test]
    fn skips_windows_containing_n() {
        // Only the two flanking ACGT runs of length 4 should be indexed.
        let seq = b"ACGTNACGT";
        let k = 4;
        let idx = KmerIndex::build(seq, k);
        // Positions 0 and 5 are the only fully-ACGT windows.
        let mut positions: Vec<u32> = idx.fwd.values().flatten().copied().collect();
        positions.sort_unstable();
        assert_eq!(positions, vec![0, 5]);
    }

    #[test]
    fn all_acgt_windows_indexed() {
        let seq = b"ACGTACGT";
        let k = 4;
        let idx = KmerIndex::build(seq, k);
        // n - k + 1 = 5 windows at positions 0..=4.
        let mut positions: Vec<u32> = idx.fwd.values().flatten().copied().collect();
        positions.sort_unstable();
        assert_eq!(positions, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn repeated_kmer_groups_positions() {
        let seq = b"AAAAAA";
        let k = 3;
        let idx = KmerIndex::build(seq, k);
        // "AAA" occurs at positions 0..=3, all under one hash.
        assert_eq!(idx.fwd.len(), 1);
        let positions = idx.fwd.values().next().unwrap();
        let mut v: Vec<u32> = positions.iter().copied().collect();
        v.sort_unstable();
        assert_eq!(v, vec![0, 1, 2, 3]);
    }

    #[test]
    fn empty_when_seq_shorter_than_k() {
        let idx = KmerIndex::build(b"ACG", 4);
        assert!(idx.fwd.is_empty());
        assert!(idx.rev.is_empty());
    }

    #[test]
    fn k_equal_one_and_large_rotations() {
        // k=1 exercises the rotate_left((k-1)=0) and rotate_left(k=1) paths.
        let seq = b"ACGT";
        let idx = KmerIndex::build(seq, 1);
        let positions: Vec<u32> = {
            let mut v: Vec<u32> = idx.fwd.values().flatten().copied().collect();
            v.sort_unstable();
            v
        };
        assert_eq!(positions, vec![0, 1, 2, 3]);
    }

    #[test]
    fn large_k_rotation_consistency() {
        // k > 64 exercises rotate amounts that wrap the 64-bit rotation.
        let seq: Vec<u8> = b"ACGT".iter().cloned().cycle().take(200).collect();
        let k = 121;
        let idx = KmerIndex::build(&seq, k);
        for (hash, positions) in &idx.fwd {
            for &pos in positions {
                let window = &seq[pos as usize..pos as usize + k];
                assert_eq!(forward_hash(window), *hash);
            }
        }
        for (hash, positions) in &idx.rev {
            for &pos in positions {
                let window = &seq[pos as usize..pos as usize + k];
                assert_eq!(reverse_hash(window), *hash);
            }
        }
    }
}
