# Low-Level Functions

These functions are implemented in Rust and exposed via PyO3.
They provide fine-grained access to the underlying FM-index and k-mer machinery.
For most use cases, the [`SequenceIndex`](sequence_index.md) class is more convenient.

## FASTA I/O

::: rusty_dot._rusty_dot.py_read_fasta

## K-mer Operations

::: rusty_dot._rusty_dot.py_build_kmer_set

::: rusty_dot._rusty_dot.py_find_kmer_coords

::: rusty_dot._rusty_dot.py_merge_kmer_runs

## PAF Formatting

::: rusty_dot._rusty_dot.py_coords_to_paf

## Index Serialization

::: rusty_dot._rusty_dot.py_save_index

::: rusty_dot._rusty_dot.py_load_index
