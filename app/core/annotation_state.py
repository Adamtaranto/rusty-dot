"""Pure helpers for the app's GFF annotation controls.

The dynamic per-feature-type toggle/colour UI needs stable Shiny input ids
(alphanumeric + underscore) derived from arbitrary GFF type names, and a
way to apply the user's selections to a parsed
:class:`~rusty_dot.annotation.GffAnnotation`.  Everything here is plain
Python so it can be unit-tested without a Shiny session.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from rusty_dot.annotation import GffAnnotation

#: Roles a GFF upload can play (query = y axis, target = x axis).
ANNOTATION_ROLES: tuple[str, str] = ('query', 'target')


def slugify(feature_type: str) -> str:
    """Turn a GFF feature-type name into a Shiny-safe input-id fragment.

    Parameters
    ----------
    feature_type : str
        Raw feature type from GFF column 3 (e.g. ``'five_prime_UTR'``,
        ``'repeat region'``).

    Returns
    -------
    str
        Lowercased ``[a-z0-9_]`` slug; non-conforming runs collapse to a
        single underscore.  Distinct types can collide after slugging
        (e.g. ``'Gene'`` vs ``'gene'``) — use :func:`type_slug_map` to
        disambiguate.
    """
    slug = re.sub(r'[^a-z0-9]+', '_', feature_type.lower()).strip('_')
    return slug or 'type'


def type_slug_map(feature_types: list[str]) -> dict[str, str]:
    """Map each feature type to a unique input-id slug.

    Parameters
    ----------
    feature_types : list[str]
        Feature type names (any order, assumed unique).

    Returns
    -------
    dict[str, str]
        ``feature_type -> slug``; collisions get a numeric suffix in input
        order, so the mapping is deterministic for a given list.
    """
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for ft in feature_types:
        slug = slugify(ft)
        candidate = slug
        n = 2
        while candidate in used:
            candidate = f'{slug}_{n}'
            n += 1
        used.add(candidate)
        mapping[ft] = candidate
    return mapping


def merge_annotations(
    annotations: list['GffAnnotation | None'],
) -> 'GffAnnotation | None':
    """Combine several parsed annotations into one.

    Used where a single annotation object is required but the features may
    have come from more than one upload — diagonal shading on a
    self-comparison, where the query and target GFFs describe the same
    sequences, and (later) merging GenBank-derived features with an
    uploaded GFF for the same role.

    Parameters
    ----------
    annotations : list[GffAnnotation or None]
        Parsed annotations; ``None`` entries are skipped.

    Returns
    -------
    GffAnnotation or None
        A new annotation over the concatenated records, carrying the
        colours of the inputs (later entries win on conflict).  ``None``
        when there is nothing to merge.  A single non-``None`` input is
        returned as-is rather than copied.
    """
    from rusty_dot.annotation import GffAnnotation  # noqa: PLC0415

    present = [a for a in annotations if a is not None]
    if not present:
        return None
    if len(present) == 1:
        return present[0]
    colors: dict[str, str] = {}
    records: list = []
    for ann in present:
        records.extend(ann.records)
        colors.update({ft: ann.get_color(ft) for ft in ann.feature_types()})
    return GffAnnotation(records, colors=colors)


def apply_annotation_config(
    annotation: 'GffAnnotation',
    enabled: dict[str, bool],
    colors: dict[str, str],
) -> 'GffAnnotation | None':
    """Apply per-type visibility and colour selections to an annotation.

    Parameters
    ----------
    annotation : GffAnnotation
        The parsed upload (never mutated — a filtered copy is returned).
    enabled : dict[str, bool]
        ``feature_type -> visible``.  Types missing from the dict stay
        visible.
    colors : dict[str, str]
        ``feature_type -> colour`` overrides (e.g. from colour pickers).
        Empty/None values are ignored.

    Returns
    -------
    GffAnnotation or None
        A filtered, recoloured copy — or ``None`` when every type is
        toggled off (nothing to draw).
    """
    keep = [ft for ft in annotation.feature_types() if enabled.get(ft, True)]
    if not keep:
        return None
    filtered = annotation.keep_feature_types(keep)
    overrides = {ft: c for ft, c in colors.items() if c}
    if overrides:
        filtered.set_colors(overrides)
    return filtered
