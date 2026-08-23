# CrossIndex

`CrossIndex` manages sequences divided into named groups and computes
cross-group pairwise comparisons.  It is compatible with
[`DotPlotter`](dotplot.md#rusty_dot.dotplot.DotPlotter).

## Workflow

Loading sequences and computing matches are **separate explicit steps**:

1. Load sequences with `add_sequence()` or `load_fasta()`.
2. Call `compute_matches()` to compute k-mer matches between groups.
3. Call `reorder_contigs()` or `reorder_for_colinearity()` (requires step 2).
4. Plot with `DotPlotter` using `query_group` / `target_group` to specify groups.
5. Optionally persist the layout with `write_fasta()` — contigs are written in
   the reordered order with reverse-oriented contigs reverse-complemented.

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

# Sort contigs for maximum collinearity.  Each query contig is assigned to its
# best-matching target chromosome (argmax of the squared match weights) and
# ordered by its gravity centre there; reverse-oriented contigs are detected.
q_sorted, t_sorted = cross.reorder_contigs()
reversed_a = cross.reversed_contigs("a")  # query contigs to render flipped

# Plot directly from CrossIndex — sequence names resolved via group params.
# Reverse-oriented contigs are rendered flipped so they read along the main
# diagonal.  Omit reverse_contigs to auto-pull the detected set for query_group.
plotter = DotPlotter(cross)
plotter.plot(
    query_group="a",   # sequences from group 'a' as rows
    target_group="b",  # sequences from group 'b' as columns
    output_path="cross_plot.png",
    reverse_contigs=reversed_a,
)
```

## Writing reordered / reoriented FASTA

`write_fasta()` persists a group's contigs in the current
[`contig_order`](#rusty_dot.paf_io.CrossIndex.contig_order), reverse-complementing
any contig flagged by [`reversed_contigs()`](#rusty_dot.paf_io.CrossIndex.reversed_contigs).
Orientation is always expressed relative to the target group, which stays the
forward reference — so reorder the group you want to *change* as the **query**.

```python
cross = CrossIndex(k=15)
cross.load_fasta("assembly_a.fasta", group="a")
cross.load_fasta("assembly_b.fasta", group="b")
cross.compute_matches(query_group="b", target_group="a")

# Reorder + reorient B against a FIXED A (A is not moved or flipped):
cross.reorder_for_colinearity("b", "a", reorder_target=False)

cross.write_fasta("assembly_a.sorted.fasta", "a")  # unchanged reference
cross.write_fasta("assembly_b.sorted.fasta", "b")  # reordered + reoriented
```

Pass `reorder_target=True` (the default) to let both assemblies be reordered;
A remains the forward reference and B's contigs are the ones reverse-complemented.
`reorder_by_length(group)` sets a length-sorted order for a group before writing.
See the *Writing Reordered FASTA* tutorial for full worked examples.

## Custom group names

```python
cross = CrossIndex(k=15)
cross.load_fasta("genome_a.fasta", group="Group_A")
cross.load_fasta("genome_b.fasta", group="Group_B")

cross.compute_matches()  # auto-detects the two groups
q_sorted, t_sorted = cross.reorder_contigs()

# Plot with explicit group names
plotter = DotPlotter(cross)
plotter.plot(
    query_group="Group_A",
    target_group="Group_B",
    output_path="cross_plot.png",
)

# Or rename groups and use explicitly
cross.rename_group("Group_A", "query")
cross.rename_group("Group_B", "target")
cross.compute_matches(query_group="query", target_group="target")
q_sorted, t_sorted = cross.reorder_contigs(query_group="query", target_group="target")
```

## Class

::: rusty_dot.paf_io.CrossIndex
