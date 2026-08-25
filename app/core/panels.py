"""Panel-grid layout helpers for the interactive dotplot.

The app always passes *explicit* ``query_names`` / ``target_names`` lists to
the plot call so that (a) the (row, col) position of every panel in the
rendered grid is known — which is what maps a double-clicked panel in the
embedded HTML report back to a (query contig, target contig) pair — and
(b) contig reordering never mutates the cached alignment objects.

:func:`resolve_orders` computes those lists for every ordering mode in
:data:`core.state.ORDER_CHOICES`, and :func:`panel_pair` performs the
(row, col) → names lookup.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterable, Mapping

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from rusty_dot.paf_io import PafRecord


def panel_pair(
    query_names: list[str],
    target_names: list[str],
    row: int,
    col: int,
) -> tuple[str, str]:
    """Map a panel grid position to its (query, target) contig names.

    Panels are laid out row-major: row *i* is ``query_names[i]`` and column
    *j* is ``target_names[j]`` — the same order the name lists were passed
    to the plot call, and the same ``rd-panel-<row>-<col>`` ids used in the
    HTML report.

    Parameters
    ----------
    query_names : list[str]
        Plotted query contig names, in row order.
    target_names : list[str]
        Plotted target contig names, in column order.
    row : int
        Panel row index (0-based, from the top).
    col : int
        Panel column index (0-based, from the left).

    Returns
    -------
    tuple[str, str]
        ``(query_name, target_name)`` for the panel.

    Raises
    ------
    IndexError
        If *row* or *col* is negative or outside the grid.
    """
    if not (0 <= row < len(query_names)):
        raise IndexError(f'Panel row {row} outside grid of {len(query_names)} rows')
    if not (0 <= col < len(target_names)):
        raise IndexError(f'Panel col {col} outside grid of {len(target_names)} cols')
    return query_names[row], target_names[col]


def resolve_orders(
    mode: str,
    records: Iterable['PafRecord'],
    query_names: list[str],
    target_names: list[str],
    lengths: Mapping[str, int],
) -> tuple[list[str], list[str], set[str]]:
    """Compute explicit axis orders (and reversals) for an ordering mode.

    Pure function: the inputs are never mutated, so cached alignment objects
    can be re-laid-out freely as the user switches modes.

    Parameters
    ----------
    mode : str
        A key of :data:`core.state.ORDER_CHOICES`: ``'input'`` (keep both
        axes as given), ``'length'`` (both axes by descending length),
        ``'colinearity'`` (d-genies gravity ordering of both axes) or
        ``'colinearity_ref'`` (target axis kept fixed, query contigs
        gravity-ordered against it).
    records : iterable of rusty_dot.paf_io.PafRecord
        Alignment records used for the gravity ordering (ignored for
        ``'input'`` / ``'length'``).
    query_names : list[str]
        Query contig names in input order.
    target_names : list[str]
        Target contig names in input order.
    lengths : mapping of str to int
        Contig name → sequence length, used for ``'length'`` ordering.

    Returns
    -------
    tuple[list[str], list[str], set[str]]
        ``(query_order, target_order, reversed_queries)``.
        *reversed_queries* is non-empty only for the colinearity modes; it
        holds query contigs detected as reverse-oriented against their
        best-matching target.

    Raises
    ------
    ValueError
        If *mode* is not a recognised ordering mode.
    """
    if mode == 'input':
        return list(query_names), list(target_names), set()
    if mode == 'length':

        def by_len(name: str) -> tuple[int, str]:
            """Sort key: descending length, then name for stability.

            Parameters
            ----------
            name : str
                Contig name.

            Returns
            -------
            tuple[int, str]
                ``(-length, name)``.
            """
            return (-lengths.get(name, 0), name)

        return sorted(query_names, key=by_len), sorted(target_names, key=by_len), set()
    if mode in ('colinearity', 'colinearity_ref'):
        from rusty_dot.paf_io import compute_gravity_contigs  # noqa: PLC0415

        q_order, t_order, reversed_q = compute_gravity_contigs(
            records,
            list(query_names),
            list(target_names),
            sort_targets=(mode == 'colinearity'),
        )
        return q_order, t_order, reversed_q
    raise ValueError(f'Unknown contig-order mode {mode!r}')


def has_self_pair(query_names: Iterable[str], target_names: Iterable[str]) -> bool:
    """Whether the panel grid contains at least one self-comparison panel.

    Diagonal feature shading only draws on panels where a contig is
    compared against itself (``DotPlotter._draw_annotation_squares`` is
    reached only for such panels), so the control that toggles it is
    pointless — and misleading — on a plain cross-assembly comparison.

    Parameters
    ----------
    query_names, target_names : Iterable[str]
        The plotted row and column contig names, as produced by
        :func:`resolve_orders`.

    Returns
    -------
    bool
        ``True`` when any contig appears on both axes.  This is a
        superset of what actually gets shaded: the plotter additionally
        requires the two sequences to be the same length before treating
        a panel as a self-comparison, so a name that coincidentally
        collides between two different assemblies enables the control
        without guaranteeing squares.
    """
    return bool(set(query_names) & set(target_names))


def nav_tips(focused: bool, multi_panel: bool) -> list[tuple[str, str]]:
    """Build the navigation tips shown under the interactive report.

    Parameters
    ----------
    focused : bool
        Whether the focused single-pair view is active.
    multi_panel : bool
        Whether the grid has more than one panel.  Click-to-focus zooms one
        panel and dims the rest, so it is meaningless — and disabled in the
        report — when there is only one panel to choose from; advertising it
        anyway sends users clicking at something that will not respond.

    Returns
    -------
    list[tuple[str, str]]
        ``(action, effect)`` pairs in display order.
    """
    tips = [
        ('scroll', 'pan up/down'),
        ('Shift+scroll', 'pan left/right'),
        ('Cmd/Ctrl+scroll', 'zoom'),
        ('drag', 'zoom to region'),
    ]
    if not focused and multi_panel:
        tips.append(('click panel', 'focus'))
    tips.append(('click match', 'details'))
    tips.append(('Esc', 'reset'))
    if not focused:
        # Double-click drill-down works even on a single-panel overview.
        tips.append(('double-click panel', 'standalone view'))
    return tips


def filter_by_min_length(
    names: Iterable[str],
    lengths: Mapping[str, int],
    min_length: int,
) -> tuple[list[str], list[str]]:
    """Split contig names into those long enough to plot and those not.

    A grid panel per contig means a handful of large chromosomes can be
    buried under hundreds of short scaffolds, each drawn as a sliver.
    Dropping the short ones is a display choice only — callers are
    expected to keep the excluded names for anything that must stay
    complete, such as the reordered-FASTA export.

    Parameters
    ----------
    names : Iterable[str]
        Contig names for one axis, in input order.
    lengths : Mapping[str, int]
        ``name -> length`` in bases.  Names missing from the mapping are
        treated as length 0 and therefore excluded by any positive
        threshold.
    min_length : int
        Minimum length to keep, in bases.  Zero or negative keeps
        everything.

    Returns
    -------
    tuple[list[str], list[str]]
        ``(kept, excluded)``, each in the input order.
    """
    if min_length <= 0:
        return list(names), []
    kept: list[str] = []
    excluded: list[str] = []
    for name in names:
        (kept if lengths.get(name, 0) >= min_length else excluded).append(name)
    return kept, excluded
