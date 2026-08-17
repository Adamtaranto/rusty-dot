# Plot Style

This module provides opt-in Nature-journal-style formatting for rusty-dot
plots: Helvetica/Arial fonts at 5&ndash;7&thinsp;pt, 0.5&thinsp;pt line and
axes widths, no top/right spines, outward ticks, and 300&thinsp;dpi
tight-bbox figure export defaults.

Use the context manager to style a single plot:

```python
from rusty_dot import DotPlotter, SequenceIndex, nature_style

idx = SequenceIndex(k=10)
idx.add_sequence("seq1", "ACGTACGTACGTACGTACGT")
idx.add_sequence("seq2", "TACGTACGTACGTACGTACG")

plotter = DotPlotter(idx)
with nature_style():
    fig = plotter.plot(output_path="dotplot.png", dpi=300)
```

!!! note
    Explicit keyword arguments always override rcParams, and
    `DotPlotter.plot` passes its own `dpi` argument (default 150) to
    `savefig`.  Pass `dpi=300` explicitly, as above, to save at Nature's
    recommended resolution; the `savefig.dpi`/`savefig.bbox` defaults in
    the style apply to plain `plt.savefig()` calls that omit those
    arguments.

Or apply the style globally for the session:

```python
from rusty_dot import use_nature_style

use_nature_style()
```

## Constants

::: rusty_dot.style.NATURE_RC

## Functions

::: rusty_dot.style.nature_style

::: rusty_dot.style.use_nature_style
