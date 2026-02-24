# PAF I/O

This module provides classes and helpers for reading, writing, and reordering
[PAF (Pairwise mApping Format)](https://github.com/lh3/miniasm/blob/master/PAF.md)
alignment records.

## CrossIdx — Multi-group cross-index

`CrossIdx` (formerly `CrossIndexPaf`) manages sequences divided into named
groups and computes cross-group pairwise comparisons.  It is compatible with
:class:`~rusty_dot.dotplot.DotPlotter`.

### Alignment scope by number of groups

* **2 groups** — alignments between the two groups only.
* **3+ groups** — all non-self ordered pairs of groups.
  Use the `group_pairs` argument of `get_paf` to restrict to specific pairs.

### DotPlotter usage

```python
from rusty_dot.paf_io import CrossIdx
from rusty_dot.dotplot import DotPlotter

cross = CrossIdx(k=15)
cross.load_fasta("assembly_a.fasta", group="a")
cross.load_fasta("assembly_b.fasta", group="b")

plotter = DotPlotter(cross)
plotter.plot(
    query_names=cross.sequence_names(group="a"),
    target_names=cross.sequence_names(group="b"),
    output_path="cross_plot.png",
)
```

::: rusty_dot.paf_io.CrossIdx

## PafAlignment — Alignment record collection

::: rusty_dot.paf_io.PafRecord

::: rusty_dot.paf_io.PafAlignment

## Functions

::: rusty_dot.paf_io.parse_paf_file

::: rusty_dot.paf_io.compute_gravity_contigs
