# PAF I/O

This module provides classes and helpers for reading, writing, and reordering
[PAF (Pairwise mApping Format)](https://github.com/lh3/miniasm/blob/master/PAF.md)
alignment records.

## PafAlignment — Alignment record collection

`PafAlignment` wraps a list of :class:`~rusty_dot.paf_io.PafRecord` objects
and provides filtering, contig reordering, and sequence-length lookup
utilities.  It can be passed directly to
:class:`~rusty_dot.dotplot.DotPlotter` — no
:class:`~rusty_dot.SequenceIndex` is required:

```python
from rusty_dot.paf_io import PafAlignment
from rusty_dot.dotplot import DotPlotter

aln = PafAlignment.from_file("alignments.paf")
q_order, t_order = aln.reorder_contigs()

plotter = DotPlotter(aln)
plotter.plot(
    query_names=q_order,
    target_names=t_order,
    output_path="dotplot.png",
    scale_sequences=True,
)
```

::: rusty_dot.paf_io.PafRecord

::: rusty_dot.paf_io.PafAlignment

## Functions

::: rusty_dot.paf_io.parse_paf_file

::: rusty_dot.paf_io.compute_gravity_contigs
