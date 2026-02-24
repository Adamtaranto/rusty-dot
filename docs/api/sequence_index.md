# SequenceIndex

The `SequenceIndex` class is the primary interface for building and querying FM-index-backed sequence comparison data.
It is implemented in Rust (via [PyO3](https://pyo3.rs)) for maximum performance.

## How multiple sequences are stored

Each sequence added to a `SequenceIndex` receives its **own independent FM-index**, built by [rust-bio](https://docs.rs/bio).
The rust-bio FM-index is a read-only data structure that is constructed once for a given sequence and cannot be extended or updated after construction.

This means:

* **Adding sequences accumulates independent indexes.** Calling `add_sequence` or `load_fasta` multiple times grows the collection — each call creates a new, isolated FM-index for that sequence only and does not affect any existing FM-index.
* **`load_fasta` can be called multiple times.** Each call parses a file and adds its sequences to the collection, preserving all sequences that were added previously. Two calls on two separate FASTA files will leave the index containing all sequences from both files.
* **Re-using a name emits a warning and overwrites.** If `add_sequence` or a FASTA record uses a name that already exists in the index, a `UserWarning` is emitted and the existing entry (and its FM-index) is replaced with a new one for the new sequence.
* **Duplicate names within a FASTA file raise an error.** If a FASTA file contains two records with the same sequence name, `load_fasta` raises a `ValueError` before adding any sequences from that file to the index.
* **Pairwise comparisons use two independent FM-indexes.** `compare_sequences` and `compare_sequences_stranded` look up the two named sequences from the dictionary and compare their individual FM-indexes — no combined or merged FM-index is ever created.

```python
from rusty_dot import SequenceIndex
import warnings

idx = SequenceIndex(k=15)

# Each call adds a new independent FM-index entry
idx.add_sequence("contig1", "ACGT" * 100)
idx.add_sequence("contig2", "TTTT" * 100)
print(idx.sequence_names())   # ['contig1', 'contig2']

# load_fasta accumulates — sequences from both files are kept
idx.load_fasta("assembly_a.fasta")
idx.load_fasta("assembly_b.fasta")
print(len(idx.sequence_names()))  # total from both files plus the two above

# Re-using a name emits a UserWarning then replaces the entry
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    idx.add_sequence("contig1", "GGGG" * 100)   # overwrites old contig1
assert len(w) == 1                               # exactly one warning emitted
assert issubclass(w[0].category, UserWarning)    # warned before overwriting
print(idx.get_sequence_length("contig1"))        # 400 (the new sequence)

# load_fasta raises ValueError if the FASTA file itself has duplicate names
try:
    idx.load_fasta("file_with_dup_names.fasta")
except ValueError as e:
    print(e)   # "duplicate sequence name 'seq1' in FASTA file '...'"
```

## Class

::: rusty_dot._rusty_dot.SequenceIndex
