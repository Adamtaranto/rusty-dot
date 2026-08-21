# rusty-dot

**rusty-dot** is a Rust + PyO3 tool for making fast dot plot comparisons of DNA sequences.

## Overview

rusty-dot provides a high-performance toolkit for pairwise DNA sequence comparison and visualisation.
At its core, it builds a compact canonical-hash [ntHash](https://github.com/bcgsc/ntHash) k-mer index for each sequence — both strands in one sorted table, built in a single O(n) pass and in parallel across sequences — and intersects two indexes with a cache-friendly two-pointer walk to efficiently find shared subsequences between any two sequences in the collection.

### Key Features

- **Fast, parallel k-mer index construction** via Rust + PyO3 bindings
- **Read FASTA / gzipped FASTA files** via [needletail](https://docs.rs/needletail)
- **Rolling-hash k-mer index** per sequence ([ntHash](https://github.com/bcgsc/ntHash)); records of a file are indexed in parallel across CPU cores ([rayon](https://docs.rs/rayon))
- **Compact canonical index** — both strands share one sorted CSR table (~12–16 bytes/bp); shared k-mers are found by a two-pointer walk and byte-verified for exact matching
- **PAF format output** for alignment records
- **Index serialization/deserialization** with [serde](https://docs.rs/serde) + postcard (sequence bytes are stored; the k-mer index is rebuilt on load)
- **All-vs-all dotplot visualization** with matplotlib: forward hits in blue, RC hits in red; edge-only axis labels in grid plots; subpanels scaled by sequence length by default (`scale_sequences=True`)
- **SVG vector output** via the `format` parameter (`format='svg'`) or by using a `.svg` file extension — suitable for publication-quality figures
- **GFF3 annotation overlays** — `GffAnnotation` (from a file, text, or raw bytes; gzip auto-detected) shades features behind self-vs-self panels and draws lane-packed side tracks with strand arrows on focused single-pair plots; features are clickable in HTML reports
- **Interactive HTML dotplot reports** (`DotPlotter.to_html()`, or an `.html` output path) — single self-contained file with click-to-focus sub-panels, scroll zoom, and a click-a-match detail bar
- **Plot-time contig ordering** — `plot(contig_order='length'|'colinearity')` with `auto_reverse=True` to flip reverse-oriented contigs automatically

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

# Colour alignments by identity from a PAF file
from rusty_dot.paf_io import PafAlignment
aln = PafAlignment.from_file("alignments.paf")
plotter = DotPlotter(idx, paf_alignment=aln)
plotter.plot(output_path="identity_dotplot.png", color_by_identity=True, identity_palette="viridis")
plotter.plot_identity_colorbar(palette="viridis", output_path="colorbar.png")
```

## Documentation Sections

- **[Web App](webapp.md)** — try rusty-dot in your browser, no installation
  required: a fully client-side assembly-comparison app (files never leave
  your machine).
- **[Installation](installation.md)** — how to install rusty-dot and its dependencies.
- **[Tutorials](tutorials/quickstart.ipynb)** — step-by-step Jupyter notebook walkthroughs.
- **[API Reference](api/index.md)** — full documentation for all classes and functions.
