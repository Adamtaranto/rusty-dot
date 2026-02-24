# API Reference

rusty-dot exposes its functionality through the following classes and functions.

## Classes

| Class | Module | Description |
|-------|--------|-------------|
| [`SequenceIndex`](sequence_index.md) | `rusty_dot` | Rust-backed FM-index for sequence comparison |
| [`DotPlotter`](dotplot.md) | `rusty_dot.dotplot` | All-vs-all dotplot visualisation |
| [`CrossIdx`](paf_io.md#rusty_dot.paf_io.CrossIdx) | `rusty_dot.paf_io` | Multi-group cross-index for cross-group pairwise comparisons; DotPlotter-compatible |
| [`PafRecord`](paf_io.md#rusty_dot.paf_io.PafRecord) | `rusty_dot.paf_io` | Single PAF alignment record |
| [`PafAlignment`](paf_io.md#rusty_dot.paf_io.PafAlignment) | `rusty_dot.paf_io` | Collection of PAF records with reordering utilities |

## Functions

| Function | Module | Description |
|----------|--------|-------------|
| [`py_read_fasta`](functions.md#rusty_dot._rusty_dot.py_read_fasta) | `rusty_dot` | Read a FASTA or gzipped FASTA file |
| [`py_build_kmer_set`](functions.md#rusty_dot._rusty_dot.py_build_kmer_set) | `rusty_dot` | Build the k-mer set for a sequence |
| [`py_find_kmer_coords`](functions.md#rusty_dot._rusty_dot.py_find_kmer_coords) | `rusty_dot` | Find k-mer positions in a sequence via FM-index |
| [`py_merge_runs`](functions.md#rusty_dot._rusty_dot.py_merge_runs) | `rusty_dot` | Unified strand-aware merge: forward and both RC patterns |
| [`py_merge_kmer_runs`](functions.md#rusty_dot._rusty_dot.py_merge_kmer_runs) | `rusty_dot` | Merge forward-strand (`+`) co-linear k-mer hits into blocks |
| [`py_merge_rev_runs`](functions.md#rusty_dot._rusty_dot.py_merge_rev_runs) | `rusty_dot` | Merge RC anti-diagonal k-mer hits (standard inverted repeat) |
| [`py_merge_rev_fwd_runs`](functions.md#rusty_dot._rusty_dot.py_merge_rev_fwd_runs) | `rusty_dot` | Merge RC co-diagonal k-mer hits (both arms same direction) |
| [`py_coords_to_paf`](functions.md#rusty_dot._rusty_dot.py_coords_to_paf) | `rusty_dot` | Convert coordinate tuples to PAF lines |
| [`py_save_index`](functions.md#rusty_dot._rusty_dot.py_save_index) | `rusty_dot` | Serialise an index collection to disk |
| [`py_load_index`](functions.md#rusty_dot._rusty_dot.py_load_index) | `rusty_dot` | Load a serialised index from disk |
| [`parse_paf_file`](paf_io.md#rusty_dot.paf_io.parse_paf_file) | `rusty_dot.paf_io` | Yield PAF records from a file |
| [`compute_gravity_contigs`](paf_io.md#rusty_dot.paf_io.compute_gravity_contigs) | `rusty_dot.paf_io` | Sort contigs by gravity centre |
