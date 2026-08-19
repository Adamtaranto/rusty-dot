[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# rusty-dot

Fast dot plot comparisons of DNA sequences using a canonical ntHash k-mer
index. Written in Rust with PyO3 python bindings.

## Browser app

[![Web App](https://img.shields.io/badge/Web%20App-live-teal)](https://adamtaranto.github.io/rusty-dot/app/)

Try rusty-dot without installing anything:
**<https://adamtaranto.github.io/rusty-dot/app/>**.
The app runs entirely in your browser (WebAssembly + Pyodide) — uploaded
assemblies never leave your machine. Align with the k-mer engine, minimap2,
LASTZ, or nucmer (or import a PAF file, or align an assembly to itself),
overlay GFF3 annotations with per-type toggles and colours, reconfigure the
dotplot without recomputing, cancel long aligner runs mid-flight, and
download SVG/PDF plots, PAF alignments, and a reordered/reoriented query
FASTA. The k-mer method accepts up to ~80 Mb of combined input in-browser;
larger genomes should use minimap2 or a local install. See
[docs/webapp.md](docs/webapp.md) for capabilities and limits, and
[app/README.md](app/README.md) to run or build it locally.

## Features

- Read FASTA / gzipped FASTA files via [needletail](https://docs.rs/needletail)
- Build a rolling-hash [ntHash](https://github.com/bcgsc/ntHash) k-mer index per
  sequence in a single O(n) pass, in parallel across sequences (via
  [rayon](https://docs.rs/rayon))
- Compact canonical-hash index (both strands in one sorted CSR table,
  ~12–16 bytes/bp) intersected with a two-pointer walk and byte-verified for
  exact matching
- **Both-strand k-mer matching**: forward (`+`) and reverse-complement (`-`)
hits detected via `compare_sequences_stranded`, with a batched
`compare_pairs_stranded(pairs, merge, min_block_len)` for whole Q×T grids
(one native call; short blocks can be filtered before they reach Python)
- Merge sequential k-mer runs into contiguous match blocks for both orientations:
  - Forward-strand co-linear diagonal merging (`py_merge_kmer_runs`)
  - RC anti-diagonal merging — standard inverted repeats (`py_merge_rev_runs`)
  - RC co-diagonal merging — both arms run in same direction (`py_merge_rev_fwd_runs`)
  - Unified strand-aware entry-point (`py_merge_runs`)
- PAF format output for alignment records
- Index serialization/deserialization with [serde](https://docs.rs/serde) +
  postcard (stores sequence bytes; the k-mer index is rebuilt on load)
- All-vs-all dotplot visualization with matplotlib:
  - Forward hits drawn in **blue** (configurable via `dot_color`)
  - Reverse-complement hits drawn in **red** (configurable via `rc_color`)
  - Sequence names rendered once — at the bottom of each column and left of each row
  - **SVG vector output** in addition to PNG/PDF via the `format` parameter
  - **Minimum alignment length filter** (`min_length`) to suppress
  short/spurious hits before rendering
- Cross-index comparisons between two sequence sets (e.g. two genome assemblies)
- Relative sequence scaling in dotplot subpanels
- d-genies-style gravity contig ordering for maximum collinearity: matches are
  weighted by `(1 + euclidean_length)²`, each contig is assigned to its single
  best-matching chromosome, and reverse-oriented contigs are detected and can be
  rendered flipped along the main diagonal (`DotPlotter.plot(reverse_contigs=…)`)
- Export the collinearity layout to FASTA (`CrossIndex.write_fasta()`): contigs
  written in the reordered order with reverse-oriented contigs reverse-complemented
- **`PafAlignment.filter_by_min_length()`** — discard short alignment records from a loaded PAF file
- **Identity colouring** for PAF-backed plots —
  `plot(color_by_identity=True, identity_palette=…, identity_colorbar=True)`
  shades match segments by percent identity with an optional colour key
- **GFF3 annotation overlays** — `GffAnnotation.from_file/from_text/from_bytes`
  (gzip auto-detected), per-type colours, diagonal shading behind alignments,
  and lane-packed side tracks with strand arrows on focused single-pair plots
- **Interactive HTML dotplot reports** (`DotPlotter.to_html()`, or an `.html`
  output path): single self-contained file with click-to-focus sub-panels,
  scroll zoom, a click-a-match detail bar, and clickable annotation features
- Opt-in **Nature-journal plot style** via the `nature_style()` context manager
  in `rusty_dot.style`
- `plot(hide_internal_axes=True)` for continuous grid plots without internal
  axis clutter
- **Plot-time contig ordering** — `plot(contig_order='length'|'colinearity')`
  with `auto_reverse=True` to flip reverse-oriented contigs automatically
- Full Python bindings via [PyO3](https://pyo3.rs)

## Installation

Requirements:

- Rust: See [rust-lang.org](https://rust-lang.org/tools/install/)
- Python >=3.9 <3.15

```bash
# Clone this project repo
git clone https://github.com/Adamtaranto/rusty-dot.git && cd rusty-dot

# Install maturin build tool
pip install maturin

# Build and install the python package
maturin develop --release
```

## Performance

rusty-dot is built for large sequences and many-contig genomes:

- **Always build with `--release`.** Debug builds are dramatically slower; the
  release profile additionally enables link-time optimisation.
- **Parallel index construction.** Building the per-sequence k-mer index for the
  records of a FASTA file, and computing all-vs-all pairwise comparisons, run
  across all available CPU cores (rayon) with the Python GIL released.
- **Rolling-hash matching on a compact canonical index.** Each sequence is
  scanned once with ntHash (forward and reverse-complement hashes rolled
  together) into a single canonical-hash table: every k-mer window is keyed by
  `min(forward, revcomp)` hash with a strand flag packed into its position,
  sorted into three flat CSR arrays. That covers **both strands in
  ~12–16 bytes per base pair** (the previous two-hash-map layout cost
  ~100 B/bp). Matching intersects two indexes with a cache-friendly
  two-pointer walk and byte-verifies representatives — no hash probing and no
  suffix-array construction on the comparison path.
- **Vector-first, high-resolution plotting.** Dotplot matches are built as NumPy
  segment arrays and drawn with one `LineCollection` per panel/strand. By default
  (`rasterized='auto'`) the match layer is **true vector** — infinitely zoomable in
  SVG/PDF — until a panel exceeds `rasterization_threshold` segments, above which
  only that layer is rasterised at `dpi` to bound file size (axes and labels
  always stay vector; raise `dpi` for a higher-resolution PNG). Enable
  `chain_gap=<bp>` to chain co-linear matches
  across gaps into a few long lines, cutting render time and file size for dense,
  genome-scale plots.

## Quick Start — single multi-FASTA index

Each sequence added to a `SequenceIndex` gets its **own independent k-mer
index**.

Calling `add_sequence` or `load_fasta` multiple times **accumulates** sequences
— it never merges or replaces the existing collection.

Re-using an existing sequence name emits a `UserWarning` and **overwrites** that
entry.

If a FASTA file contains duplicate sequence names, `load_fasta` raises a
`ValueError` before adding any sequences.

```python
from rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter

# Build an index from a multi-sequence FASTA file
# Each sequence in the file gets its own independent k-mer index entry
idx = SequenceIndex(k=15)
names = idx.load_fasta("assembly.fasta")

# load_fasta accumulates: calling it again adds more sequences, keeps existing ones
# idx.load_fasta("more_sequences.fasta")   # would add to the same index

# List the sequences now held in the index
print(idx.sequence_names())   # ['contig1', 'contig2', 'contig3', ...]

# Print all pairwise PAF lines (every i ≠ j combination)
for line in idx.get_paf_all():
    print(line)

# Print PAF lines for one specific pair
for line in idx.get_paf("contig1", "contig2"):
    print(line)

# All-vs-all dotplot
# Forward (+) hits are drawn in blue, reverse-complement (-) hits in red.
# Sequence names appear once per column (bottom) and once per row (left).
plotter = DotPlotter(idx)
plotter.plot(output_path="all_vs_all.png", title="All vs All")

# Save as an SVG vector image instead of PNG
plotter.plot(output_path="all_vs_all.svg", title="All vs All")

# Filter out short alignments (< 500 bp) before plotting
plotter.plot(output_path="filtered.png", min_length=500)

# Single pairwise dotplot
plotter.plot_single("contig1", "contig2", output_path="pair.png")
```

## All-vs-All Dotplot Between Two Genomes

Compare sequences from two separate FASTA files (e.g. two genome assemblies) and
plot an all-vs-all grid with subpanels scaled by relative sequence length.

```python
from rusty_dot.dotplot import DotPlotter
from rusty_dot.paf_io import CrossIndex, PafAlignment, PafRecord

# --- Build a cross-index for two assemblies ---
cross = CrossIndex(k=15)
cross.load_fasta("genome_a.fasta", group="a")   # query sequences (rows)
cross.load_fasta("genome_b.fasta", group="b")   # target sequences (columns)

# Compute and cache the cross-group matches.  This is the primary
# computation step and is required before any reordering; pass
# min_block_len= to drop short match blocks at compute time.
cross.compute_matches()

# --- Sort contigs for maximum collinearity ---
# Option 1: via CrossIndex (the d-genies gravity algorithm; its result
# matches SequenceIndex.optimal_contig_order by construction).  Each query
# contig is assigned to its best-matching target chromosome and ordered by
# its gravity centre there; reverse-oriented contigs are detected and
# exposed via cross.reversed_contigs("a").
q_sorted, t_sorted = cross.reorder_contigs()
reversed_a = cross.reversed_contigs("a")  # names to render flipped

# Option 2: via PafAlignment gravity-centre algorithm
# Retrieve all cross-group PAF lines
paf_lines = cross.get_paf_all()

records = [PafRecord.from_line(line) for line in paf_lines]
aln = PafAlignment.from_records(records)
q_sorted, t_sorted = aln.reorder_contigs(
    query_names=cross.query_names,
    target_names=cross.target_names,
)
# Unmatched contigs are placed at the end, sorted by descending length.

# --- Plot with relative scaling ---
# For a CrossIndex, select each axis by group; contig_order='colinearity'
# applies the gravity ordering at plot time, and reverse-oriented query
# contigs are auto-detected (pass reverse_contigs= to override).

plotter = DotPlotter(cross)
plotter.plot(
    query_group="a",
    target_group="b",
    contig_order="colinearity",
    output_path="cross_dotplot.png",
    scale_sequences=True,   # subplot size proportional to sequence length
    title="Genome A vs Genome B",
    # Render reverse-oriented query contigs flipped so they read along the main
    # diagonal.  Pass an explicit set, or omit to auto-pull the detected set.
    reverse_contigs=reversed_a,
)

# Save as SVG vector image for publication-quality output
plotter.plot(
    query_group="a",
    target_group="b",
    contig_order="colinearity",
    output_path="cross_dotplot.svg",
    scale_sequences=True,
    title="Genome A vs Genome B",
)

# Suppress short alignments (e.g. < 500 bp) from the plot
plotter.plot(
    query_group="a",
    target_group="b",
    contig_order="colinearity",
    output_path="cross_dotplot_filtered.png",
    scale_sequences=True,
    min_length=500,
    title="Genome A vs Genome B (≥500 bp alignments)",
)

# --- Persist the collinearity layout as FASTA ---
# Reorder + reorient assembly B against a FIXED assembly A, then write both.
# Reverse-oriented B contigs are reverse-complemented on write; A is untouched.
# (compute_matches must have been called for the (query, target) pair first.)
cross.compute_matches(query_group="b", target_group="a")
cross.reorder_for_colinearity("b", "a", reorder_target=False)
cross.write_fasta("assembly_a.sorted.fasta", "a")  # forward reference, unchanged
cross.write_fasta("assembly_b.sorted.fasta", "b")  # reordered + reoriented
```

## GFF Annotation Overlays

Overlay GFF3 features on dotplots. Feature types are auto-coloured from a
qualitative palette (override per type), and a colour legend is added
automatically.

```python
from rusty_dot import DotPlotter, SequenceIndex
from rusty_dot.annotation import GffAnnotation

# From a file (gzip detected automatically), text, or raw bytes.
ann = GffAnnotation.from_file("features.gff3.gz")
ann = ann.keep_feature_types(["gene", "CDS", "repeat_region"])
ann.set_colors({"gene": "#2c7fb8"})

plotter = DotPlotter(idx)

# Diagonal squares: features shade self-vs-self panels behind the alignments.
plotter.plot(annotation=ann, output_path="annotated_grid.png")

# Focused single-pair view with side annotation tracks: lane-packed feature
# shapes left of the y axis and below the x axis, direction arrows for
# stranded types (gene/mRNA/exon/CDS/ORF), and connector lines joining
# multi-part CDS groups.
plotter.plot(
    query_names=["chr1"],
    target_names=["chr2"],
    annotation_query=ann,
    annotation_target=ann,
    annotation_tracks=True,
    output_path="annotated_pair.png",
)

# Interactive HTML reports make diagonal features clickable (name, type,
# coordinates, strand, parent shown in the detail bar).
plotter.to_html("report.html", annotation=ann)
```

See the dotplot tutorial for rendered examples.

## Filtering PAF Alignments by Length

Use `PafAlignment.filter_by_min_length` to remove short alignment records after
loading a PAF file. This is particularly useful for cleaned-up visualisations
when alignments have been merged from k-mer runs (which can be longer than the
k-mer size) or when working with a pre-computed PAF file.

```python
from rusty_dot.paf_io import PafAlignment

aln = PafAlignment.from_file("alignments.paf")

# Keep only alignments of at least 500 bp on the query
aln_long = aln.filter_by_min_length(500)
print(f"Records before: {len(aln)}, after: {len(aln_long)}")
```

## Writing PAF Lines to a File

```python
# All pairwise alignments within a single index
paf_lines = idx.get_paf_all()

# Or one specific pair
paf_lines = idx.get_paf("contig1", "contig2", merge=True)

with open("alignments.paf", "w") as f:
    for line in paf_lines:
        f.write(line + "\n")
```

## Saving and Loading Indexes

```python
# Save the current index to a compact binary file
idx.save("my_index.bin")

# Load into a new index (k must match the saved index)
idx2 = SequenceIndex(k=15)
idx2.load("my_index.bin")
```
