# Low-Level Functions

These functions are implemented in Rust and exposed via PyO3.
They provide fine-grained access to the underlying FM-index and k-mer machinery.
For most use cases, the [`SequenceIndex`](sequence_index.md) class is more convenient.

## FASTA I/O

::: rusty_dot._rusty_dot.py_read_fasta

## K-mer Operations

::: rusty_dot._rusty_dot.py_build_kmer_set

::: rusty_dot._rusty_dot.py_find_kmer_coords

## Merging K-mer Runs

rusty-dot provides four merge functions covering all k-mer alignment orientations.
`py_merge_runs` is the recommended entry-point for new code; the strand-specific
functions are available for lower-level control.

### Unified entry-point

::: rusty_dot._rusty_dot.py_merge_runs

### Forward-strand merge

::: rusty_dot._rusty_dot.py_merge_kmer_runs

### Reverse-complement merges

Two complementary algorithms cover all reverse-complement alignment patterns:

| Pattern | When to use |
|---------|-------------|
| `py_merge_rev_runs` | RC target positions *decrease* as query advances (query +1, target −1 per step) — standard inverted-repeat alignment where the two arms face each other |
| `py_merge_rev_fwd_runs` | RC target positions *increase* as query advances (query +1, target +1 per step) — both repeat arms run in the same left-to-right direction |

`py_merge_runs(strand="-")` calls both and deduplicates the results automatically.

::: rusty_dot._rusty_dot.py_merge_rev_runs

::: rusty_dot._rusty_dot.py_merge_rev_fwd_runs

## PAF Formatting

::: rusty_dot._rusty_dot.py_coords_to_paf

## Index Serialization

::: rusty_dot._rusty_dot.py_save_index

::: rusty_dot._rusty_dot.py_load_index
