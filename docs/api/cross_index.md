# CrossIndex

`CrossIndex` manages sequences divided into named groups and computes
cross-group pairwise comparisons.  It is compatible with
:class:`~rusty_dot.dotplot.DotPlotter`.

## Workflow

Loading sequences and computing matches are **separate explicit steps**:

1. Load sequences with `add_sequence()` or `load_fasta()`.
2. Call `compute_matches()` to compute k-mer matches between groups.
3. Call `reorder_contigs()` or `reorder_for_colinearity()` (requires step 2).

Progress is logged at `INFO` level for each loading and computation step.
Warnings are emitted when a sequence name already exists in the same or
another group.

## Alignment scope by number of groups

* **2 groups** — `compute_matches()` compares the two groups.
* **3+ groups** — `compute_matches()` computes all non-self ordered pairs by
  default.  Use the `query_group` / `target_group` arguments to restrict to a
  specific pair.

## Quick start

```python
from rusty_dot.paf_io import CrossIndex
from rusty_dot.dotplot import DotPlotter

cross = CrossIndex(k=15)
cross.load_fasta("assembly_a.fasta", group="a")
cross.load_fasta("assembly_b.fasta", group="b")

# Explicitly compute k-mer matches (required before reorder_contigs)
cross.compute_matches()
print("Computed pairs:", cross.computed_group_pairs)

# Sort contigs for maximum collinearity
q_sorted, t_sorted = cross.reorder_contigs()

# Plot
plotter = DotPlotter(cross)
plotter.plot(
    query_names=cross.sequence_names(group="a"),
    target_names=cross.sequence_names(group="b"),
    output_path="cross_plot.png",
)
```

## Custom group names

```python
cross = CrossIndex(k=15)
cross.load_fasta("genome_a.fasta", group="Group_A")
cross.load_fasta("genome_b.fasta", group="Group_B")

cross.compute_matches()  # auto-detects the two groups
q_sorted, t_sorted = cross.reorder_contigs()

# Or rename groups and use explicitly
cross.rename_group("Group_A", "query")
cross.rename_group("Group_B", "target")
cross.compute_matches(query_group="query", target_group="target")
q_sorted, t_sorted = cross.reorder_contigs(query_group="query", target_group="target")
```

## Class

::: rusty_dot.paf_io.CrossIndex
