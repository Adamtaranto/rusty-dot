# CrossIndex

`CrossIndex` manages sequences divided into named groups and computes
cross-group pairwise comparisons.  It is compatible with
:class:`~rusty_dot.dotplot.DotPlotter`.

## Alignment scope by number of groups

* **2 groups** — alignments between the two groups only.
* **3+ groups** — all non-self ordered pairs of groups.
  Use the `group_pairs` argument of `get_paf` to restrict to specific pairs.

## DotPlotter usage

```python
from rusty_dot.paf_io import CrossIndex
from rusty_dot.dotplot import DotPlotter

cross = CrossIndex(k=15)
cross.load_fasta("assembly_a.fasta", group="a")
cross.load_fasta("assembly_b.fasta", group="b")

plotter = DotPlotter(cross)
plotter.plot(
    query_names=cross.sequence_names(group="a"),
    target_names=cross.sequence_names(group="b"),
    output_path="cross_plot.png",
)
```

## Class

::: rusty_dot.paf_io.CrossIndex
