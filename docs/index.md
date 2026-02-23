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
- **Merge sequential k-mer runs** into contiguous match blocks
- **PAF format output** for alignment records
- **FM-index serialization/deserialization** with [serde](https://docs.rs/serde) + bincode
- **All-vs-all dotplot visualization** with matplotlib
- **Full Python bindings** via [PyO3](https://pyo3.rs)

## Quick Start

```python
from rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter

# Build index for two sequences
idx = SequenceIndex(k=15)
idx.load_fasta("genome1.fasta")
idx.load_fasta("genome2.fasta")

# Get PAF-format alignments
for line in idx.get_paf("seq1", "seq2"):
    print(line)

# Generate dotplot
plotter = DotPlotter(idx)
plotter.plot(output_path="dotplot.png")
```

## Documentation Sections

- **[Installation](installation.md)** — how to install rusty-dot and its dependencies.
- **[Tutorials](tutorials/quickstart.ipynb)** — step-by-step Jupyter notebook walkthroughs.
- **[API Reference](api/index.md)** — full documentation for all classes and functions.
