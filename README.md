[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# rusty-dot

Fast dot plot comparisons of DNA sequences using an FM-Index.
Written in Rust with PyO3 python bindings.

## Features

- Read FASTA / gzipped FASTA files via [needletail](https://docs.rs/needletail)
- Build FM-indexes per sequence using [rust-bio](https://docs.rs/bio)
- K-mer set intersection for efficient shared k-mer lookup
- Merge sequential k-mer runs into contiguous match blocks
- PAF format output for alignment records
- FM-index serialization/deserialization with [serde](https://docs.rs/serde) + bincode
- All-vs-all dotplot visualization with matplotlib
- Full Python bindings via [PyO3](https://pyo3.rs)

## Installation

```bash
pip install maturin
maturin develop --release
```

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

## Loading from FASTA files

```python
from rusty_dot import SequenceIndex, py_read_fasta

# Low-level: read FASTA (plain or gzipped)
seqs = py_read_fasta("sequences.fasta.gz")
print(seqs.keys())  # sequence names

# High-level: build index directly from FASTA
idx = SequenceIndex(k=12)
names = idx.load_fasta("sequences.fasta")
print(f"Indexed {len(names)} sequences")
```

## Saving and Loading Indexes

```python
# Save
idx.save("my_index.bin")

# Load into a new index (k must match)
idx2 = SequenceIndex(k=12)
idx2.load("my_index.bin")
```

## PAF Output

```python
paf_lines = idx.get_paf("query_seq", "target_seq", merge=True)
with open("alignments.paf", "w") as f:
    for line in paf_lines:
        f.write(line + "\n")
```

## All-vs-All Dotplot

```python
from rusty_dot.dotplot import DotPlotter

plotter = DotPlotter(idx)

# All sequences vs all sequences
plotter.plot(output_path="all_vs_all.png", title="All vs All Dotplot")

# Specific query vs target sets
plotter.plot(
    query_names=["seq1", "seq2"],
    target_names=["seq3", "seq4"],
    output_path="subset.png",
)

# Single pair
plotter.plot_single("seq1", "seq2", output_path="pair.png")
```

## Running Tests

```bash
pip install maturin pytest matplotlib numpy
maturin develop --extras dev,docs
pytest tests/ -v
```
