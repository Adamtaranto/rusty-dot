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
