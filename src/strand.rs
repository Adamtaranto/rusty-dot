//! Strand-related utilities: reverse complement and strand constants.

/// The forward (`+`) strand byte.
pub const STRAND_FWD: u8 = b'+';
/// The reverse (`-`) strand byte.
pub const STRAND_REV: u8 = b'-';

/// Compute the complement of a single uppercase DNA byte.
///
/// `A ↔ T`, `C ↔ G`, `N → N`.  All other characters are mapped to `N`.
///
/// # Arguments
///
/// * `b` - A single uppercase DNA base byte.
///
/// # Returns
///
/// The complementary base byte.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'N' => b'N',
        _ => b'N',
    }
}

/// Compute the reverse complement of a DNA byte slice.
///
/// The input is reversed and each byte is complemented.  Non-ACGT characters
/// are mapped to `N`.
///
/// # Arguments
///
/// * `seq` - A slice of uppercase DNA bytes.
///
/// # Returns
///
/// A new `Vec<u8>` containing the reverse complement.
///
/// # Examples
///
/// ```
/// use _rusty_dot::strand::revcomp;
/// assert_eq!(revcomp(b"ACGT"), b"ACGT");  // palindrome
/// assert_eq!(revcomp(b"AAAC"), b"GTTT");
/// ```
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_revcomp_palindrome() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
    }

    #[test]
    fn test_revcomp_basic() {
        assert_eq!(revcomp(b"AAAC"), b"GTTT");
        assert_eq!(revcomp(b"TTTT"), b"AAAA");
    }

    #[test]
    fn test_revcomp_n() {
        assert_eq!(revcomp(b"ACNT"), b"ANGT");
    }

    #[test]
    fn test_revcomp_roundtrip() {
        let seq = b"ACGTTGCA";
        assert_eq!(revcomp(&revcomp(seq)), seq);
    }

    #[test]
    fn test_complement_all_bases() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }
}
