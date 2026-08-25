"""Shared, case-insensitive colour assignment for GFF feature types.

:class:`~rusty_dot.annotation.GffAnnotation` colours its own types by their
index into its own sorted type list, so the *same* type lands on different
palette entries in two different uploads -- ``gene`` green on the x track
and blue on the y track -- and ``CDS`` vs ``cds`` are treated as unrelated.
Fixing that in the library would recolour every existing user's plots, so
the app assigns colours itself and passes them down explicitly (the library
already supports ``GffAnnotation(colors=...)`` / ``set_colors``).

The rules here are:

* types are matched case-insensitively, with a small alias table so
  ``repeat_region`` and ``TE`` share the repeat colour;
* a handful of common biological types get reserved, conventional colours;
* everything else walks a 24-entry palette chosen to stay distinguishable
  from the reserved five, falling back to deterministic generated colours
  rather than silently repeating.

Assignment is over the *union* of types across all uploads, so query and
target always agree.  Everything is plain Python: no ``rusty_dot`` import,
so this stays unit-testable without a Shiny session or the compiled
extension.
"""

from __future__ import annotations

import colorsys
import hashlib

#: Conventional colours for the types users look for first.  Keys are
#: normalised (see :func:`normalise_type`).
RESERVED: dict[str, str] = {
    'gene': '#2ca02c',  # green
    'cds': '#e6c229',  # yellow
    'mrna': '#7b1e3a',  # maroon
    'exon': '#8c8c8c',  # grey
    'repeat': '#d62728',  # red
}

#: Display-only spelling variants folded onto a reserved key.  Deliberately
#: short and explicit: this affects colour only, never parsing, filtering or
#: what the plot reports as a feature's type.
_ALIASES: dict[str, str] = {
    'transcript': 'mrna',
    'messenger_rna': 'mrna',
    'coding_sequence': 'cds',
    'repeat_region': 'repeat',
    'repeat_unit': 'repeat',
    'dispersed_repeat': 'repeat',
    'tandem_repeat': 'repeat',
    'inverted_repeat': 'repeat',
    'ltr': 'repeat',
    'te': 'repeat',
    'transposable_element': 'repeat',
}

#: Non-repeating palette for everything else.  Ordered so neighbouring
#: types stay easy to tell apart, and screened against RESERVED.
PALETTE: tuple[str, ...] = (
    '#1f77b4',  # blue
    '#ff7f0e',  # orange
    '#9467bd',  # purple
    '#17becf',  # cyan
    '#8c564b',  # brown
    '#e377c2',  # pink
    '#bcbd22',  # olive
    '#393b79',  # navy
    '#5254a3',  # indigo
    '#637939',  # moss
    '#8c6d31',  # bronze
    '#843c39',  # brick
    '#7b4173',  # plum
    '#3182bd',  # steel
    '#31a354',  # emerald
    '#756bb1',  # violet
    '#636363',  # slate
    '#e6550d',  # tangerine
    '#fd8d3c',  # apricot
    '#74c476',  # sage
    '#9e9ac8',  # lilac
    '#969696',  # ash
    '#ad494a',  # rust
    '#a55194',  # orchid
)


def normalise_type(feature_type: str) -> str:
    """Fold a raw GFF type name to its colour key.

    Parameters
    ----------
    feature_type : str
        Raw feature type from GFF column 3 (e.g. ``'CDS'``, ``'mRNA'``,
        ``'repeat_region'``).

    Returns
    -------
    str
        Case-folded name, with alias variants mapped onto their reserved
        key.  ``'CDS'`` and ``'cds'`` both give ``'cds'``; ``'TE'`` and
        ``'repeat_region'`` both give ``'repeat'``.
    """
    key = feature_type.strip().casefold()
    return _ALIASES.get(key, key)


def _generated_color(key: str) -> str:
    """Deterministic fallback colour for a type past the palette's end.

    Derived from a hash of the (normalised) name, so it is stable across
    sessions and independent of how many other types are present -- and
    never silently repeats a palette or reserved colour by coincidence of
    ordering.
    """
    digest = hashlib.sha256(key.encode('utf-8')).digest()
    hue = digest[0] / 255.0
    sat = 0.45 + (digest[1] / 255.0) * 0.35
    val = 0.55 + (digest[2] / 255.0) * 0.30
    r, g, b = colorsys.hsv_to_rgb(hue, sat, val)
    return '#{:02x}{:02x}{:02x}'.format(int(r * 255), int(g * 255), int(b * 255))


def assign_shared_colors(types_by_role: dict[str, list[str]]) -> dict[str, str]:
    """Assign one colour per normalised type across every role.

    Parameters
    ----------
    types_by_role : dict[str, list[str]]
        ``role -> raw feature types`` for each upload (query, target, ...).

    Returns
    -------
    dict[str, str]
        ``normalised type -> hex colour``.  Reserved types keep their
        conventional colour; the rest take palette entries in sorted-name
        order, so the result is deterministic and independent of which
        role a type came from or the order roles are listed in.
    """
    keys = sorted({normalise_type(ft) for fts in types_by_role.values() for ft in fts})
    colors: dict[str, str] = {}
    palette_i = 0
    for key in keys:
        if key in RESERVED:
            colors[key] = RESERVED[key]
        elif palette_i < len(PALETTE):
            colors[key] = PALETTE[palette_i]
            palette_i += 1
        else:
            colors[key] = _generated_color(key)
    return colors


def color_map_for(feature_types: list[str], shared: dict[str, str]) -> dict[str, str]:
    """Expand shared colours back onto one annotation's raw type names.

    Parameters
    ----------
    feature_types : list[str]
        Raw type names as they appear in that upload.
    shared : dict[str, str]
        Output of :func:`assign_shared_colors`.

    Returns
    -------
    dict[str, str]
        ``raw feature type -> hex colour``, ready for
        :meth:`~rusty_dot.annotation.GffAnnotation.set_colors`.  Types
        absent from *shared* are omitted rather than defaulted, leaving the
        library's own colour in place.
    """
    out: dict[str, str] = {}
    for ft in feature_types:
        color = shared.get(normalise_type(ft))
        if color:
            out[ft] = color
    return out


def display_name(raw_names: list[str]) -> str:
    """Pick the label to show for a group of spellings of one type.

    Parameters
    ----------
    raw_names : list[str]
        Every raw spelling seen for one normalised type (e.g.
        ``['CDS', 'cds']``).

    Returns
    -------
    str
        The most common spelling, ties broken alphabetically so the label
        does not flicker between uploads.
    """
    counts: dict[str, int] = {}
    for name in raw_names:
        counts[name] = counts.get(name, 0) + 1
    return min(counts, key=lambda n: (-counts[n], n))
