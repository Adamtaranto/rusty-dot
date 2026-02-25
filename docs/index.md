# rusty-dot

**rusty-dot** is a Rust + PyO3 tool for making fast dot plot comparisons of DNA sequences using a Rust FM-Index.

## Overview

rusty-dot provides a high-performance toolkit for pairwise DNA sequence comparison and visualisation.
At its core, it builds an [FM-index](https://en.wikipedia.org/wiki/FM-index) (via [rust-bio](https://docs.rs/bio)) for each sequence and uses k-mer set intersection to efficiently find shared subsequences between any two sequences in the collection.

### Key Features

- **Fast FM-index construction** via Rust + PyO3 bindings
- **Read FASTA / gzipped FASTA files** via [needletail](https://docs.rs/needletail)
- **Build FM-indexes** per sequence using [rust-bio](https://docs.rs/bio)
- **K-mer set intersection** for efficient shared k-mer lookup
- **Both-strand k-mer matching**: forward (`+`) and reverse-complement (`-`) hits via `compare_sequences_stranded`
- **Complete RC hit coverage**: two patterns merged independently — anti-diagonal (standard inverted repeat) and co-diagonal (both arms same direction)
- **Unified merge API** (`py_merge_runs`) handles all orientation cases with a single call
- **PAF format output** for alignment records
- **FM-index serialization/deserialization** with [serde](https://docs.rs/serde) + bincode
- **All-vs-all dotplot visualization** with matplotlib: forward hits in blue, RC hits in red; edge-only axis labels in grid plots; subpanels scaled by sequence length by default (`scale_sequences=True`)
- **SVG vector output** via the `format` parameter (`format='svg'`) or by using a `.svg` file extension — suitable for publication-quality figures
- **Minimum alignment length filter** (`min_length`) on `DotPlotter.plot()` / `plot_single()` — suppresses short or spurious alignment hits before rendering
- **`CrossIdx`** multi-group cross-index: N arbitrary sequence groups, configurable group pairs for alignment, per-group contig ordering (insertion order, length, or collinearity), `run_merge` to update cached PAF records, compatible with `DotPlotter`
- **`PafAlignment.filter_by_min_length()`** — discard short alignment records from a loaded PAF file; filters on query aligned length
- **Full Python bindings** via [PyO3](https://pyo3.rs)

## Quick Start

```python
from rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter

# Build index for two sequences
idx = SequenceIndex(k=15)
idx.load_fasta("genome1.fasta")
idx.load_fasta("genome2.fasta")

# Get PAF-format alignments (forward strand only)
for line in idx.get_paf("seq1", "seq2"):
    print(line)

# Stranded comparison: forward (+) and reverse-complement (-) hits
hits = idx.compare_sequences_stranded("seq1", "seq2", merge=True)
for qs, qe, ts, te, strand in hits:
    print(f"{strand}  q[{qs}:{qe}]  t[{ts}:{te}]")

# Generate dotplot — forward hits blue, RC hits red
plotter = DotPlotter(idx)
plotter.plot(output_path="dotplot.png")

# Save as SVG vector image
plotter.plot(output_path="dotplot.svg")

# Filter short alignments (< 200 bp) before plotting
plotter.plot(output_path="dotplot_filtered.png", min_length=200)
```

## Documentation Sections

- **[Installation](installation.md)** — how to install rusty-dot and its dependencies.
- **[Tutorials](tutorials/quickstart.ipynb)** — step-by-step Jupyter notebook walkthroughs.
- **[API Reference](api/index.md)** — full documentation for all classes and functions.
