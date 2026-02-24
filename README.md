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
- Cross-index comparisons between two sequence sets (e.g. two genome assemblies)
- Relative sequence scaling in dotplot subpanels
- Gravity-centre contig ordering for maximum collinearity
- Full Python bindings via [PyO3](https://pyo3.rs)

## Installation

Requirements:

- Rust: See [rust-lang.org](https://rust-lang.org/tools/install/)
- Python >=3.9 <3.14

```bash
# Clone this project repo
git clone https://github.com/Adamtaranto/rusty-dot.git && cd rusty-dot

# Install maturin build tool
pip install maturin

# Build and install the python package
maturin develop --release
```

## Quick Start — single multi-FASTA index

```python
from rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter

# Build an index from a multi-sequence FASTA file
idx = SequenceIndex(k=15)
names = idx.load_fasta("assembly.fasta")

# List the sequences now held in the index (also available as idx.sequence_names())
print(idx.sequence_names())   # ['contig1', 'contig2', 'contig3', ...]

# Print all pairwise PAF lines (every i ≠ j combination)
for line in idx.get_paf_all():
    print(line)

# Print PAF lines for one specific pair
for line in idx.get_paf("contig1", "contig2"):
    print(line)

# All-vs-all dotplot
plotter = DotPlotter(idx)
plotter.plot(output_path="all_vs_all.png", title="All vs All")

# Single pairwise dotplot
plotter.plot_single("contig1", "contig2", output_path="pair.png")
```

## All-vs-All Dotplot Between Two Indexes

Compare sequences from two separate FASTA files (e.g. two genome assemblies) and
plot an all-vs-all grid with subpanels scaled by relative sequence length.

```python
from rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter
from rusty_dot.paf_io import CrossIndexPaf, PafAlignment, PafRecord

# --- Build a cross-index for two assemblies ---
cross = CrossIndexPaf(k=15)
cross.load_fasta("genome_a.fasta", group="a")   # query sequences (rows)
cross.load_fasta("genome_b.fasta", group="b")   # target sequences (columns)

# Retrieve all cross-group PAF lines
paf_lines = cross.get_paf_all()

# --- Sort contigs for maximum collinearity ---
# Option 1: via CrossIndexPaf (delegates to SequenceIndex.optimal_contig_order)
q_sorted, t_sorted = cross.reorder_contigs()

# Option 2: via PafAlignment gravity-centre algorithm
records = [PafRecord.from_line(line) for line in paf_lines]
aln = PafAlignment.from_records(records)
q_sorted, t_sorted = aln.reorder_contigs(
    query_names=cross.query_names,
    target_names=cross.target_names,
)
# Unmatched contigs are placed at the end, sorted by descending length.

# --- Plot with relative scaling ---
# Build a combined SequenceIndex containing sequences from both assemblies
idx = SequenceIndex(k=15)
idx.load_fasta("genome_a.fasta")
idx.load_fasta("genome_b.fasta")

plotter = DotPlotter(idx)
plotter.plot(
    query_names=q_sorted,
    target_names=t_sorted,
    output_path="cross_dotplot.png",
    scale_sequences=True,   # subplot size proportional to sequence length
    title="Genome A vs Genome B",
)
```

## Saving and Loading Indexes

```python
# Save the current index to a compact binary file
idx.save("my_index.bin")

# Load into a new index (k must match the saved index)
idx2 = SequenceIndex(k=15)
idx2.load("my_index.bin")
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