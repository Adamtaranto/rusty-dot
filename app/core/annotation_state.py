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


def build_feature_rows(
    annotation: 'GffAnnotation | None',
    seqname: str,
    role: str,
) -> list[dict]:
    """Describe every feature on one sequence for the annotations table.

    Rows are positional over ``annotation.records``, so ``uid`` is stable
    for the lifetime of an upload set and needs no hashing.  Re-uploading
    a file renumbers them, which is why the app clears its per-feature
    overrides whenever a role's sources change.

    Parameters
    ----------
    annotation : GffAnnotation or None
        The merged annotation for *role*.
    seqname : str
        Only features on this sequence are listed.
    role : str
        ``'query'`` or ``'target'``; part of the uid and shown as a column.

    Returns
    -------
    list[dict]
        One row per feature, in file order.  ``start``/``end`` are
        **1-based inclusive** for display, matching the source file rather
        than the library's internal half-open coordinates.
    """
    if annotation is None:
        return []
    rows = []
    for i, rec in enumerate(annotation.records):
        if rec.seqname != seqname:
            continue
        attrs = rec.attr_dict()
        extra = {
            k: v for k, v in attrs.items() if k not in ('ID', 'Parent', 'Name') and v
        }
        rows.append(
            {
                'uid': f'{role}:{i}',
                'role': role,
                'source_file': getattr(rec, 'source_file', '') or '',
                'source': rec.source,
                'type': rec.feature_type,
                'seqname': rec.seqname,
                'start': rec.start + 1,
                'end': rec.end,
                'length': rec.end - rec.start,
                'strand': rec.strand,
                'score': rec.score,
                'id': rec.feature_id,
                'parent': rec.parent,
                'name': rec.name,
                'attributes': extra,
            }
        )
    return rows


def apply_feature_overrides(
    annotation: 'GffAnnotation | None',
    hidden_uids: frozenset[str],
    colors: dict[str, str],
    role: str,
) -> 'GffAnnotation | None':
    """Apply per-feature visibility and colour choices.

    Must be applied to the **unfiltered** merged annotation — the same one
    :func:`build_feature_rows` described — because uids are positional
    over ``records``, and filtering by type first would renumber them.
    Run :func:`apply_annotation_config` *after* this, so a type toggled
    off still wins over a feature toggled on; ``keep_feature_types``
    reuses the record objects, so the colours written here survive it.

    Colours are written onto the records themselves
    (``GffFeature.color``), which the drawing code prefers over the
    per-type colour.

    Parameters
    ----------
    annotation : GffAnnotation or None
        The merged, *unfiltered* annotation for *role*.
    hidden_uids : frozenset[str]
        ``role:index`` uids the user has switched off.
    colors : dict[str, str]
        ``uid -> hex colour`` overrides.
    role : str
        The role whose uids these are.

    Returns
    -------
    GffAnnotation or None
        A filtered copy, or ``None`` when nothing remains visible.  The
        input is returned unchanged when nothing is hidden (colours are
        applied in place either way).
    """
    from rusty_dot.annotation import GffAnnotation  # noqa: PLC0415

    if annotation is None:
        return None
    kept = []
    for i, rec in enumerate(annotation.records):
        uid = f'{role}:{i}'
        if uid in hidden_uids:
            continue
        override = colors.get(uid)
        if override:
            rec.color = override
        elif rec.color is not None:
            rec.color = None  # reset-to-type-colour
        kept.append(rec)
    if not kept:
        return None
    if len(kept) == len(annotation.records):
        return annotation
    return GffAnnotation(
        kept, colors={ft: annotation.get_color(ft) for ft in annotation.feature_types()}
    )


def replace_source(
    entries: tuple,
    kind: str,
    filename: str,
    annotation: 'GffAnnotation | None',
    key,
) -> tuple | None:
    """Add, replace or clear one annotation source, or report no change.

    A role can hold several sources at once (a GenBank file's own features
    plus an uploaded GFF), so only the entry of the matching *kind* is
    touched.

    Returning ``None`` for an unchanged source matters: the app's reactive
    values invalidate on *identity*, so re-setting an equal-but-new tuple
    would rebuild the feature-type controls and discard the user's toggles
    and colours on every run.

    Parameters
    ----------
    entries : tuple
        Current sources for the role, each
        ``{'kind', 'filename', 'annotation', 'key'}``.
    kind : str
        ``'gff'`` or ``'genbank'`` — which slot this call owns.
    filename : str
        Name of the uploaded file, recorded for display.
    annotation : GffAnnotation or None
        Parsed features, or ``None`` to clear this kind.
    key : hashable
        Identity of the upload, normally ``(content_digest, filename)``.
        The filename is part of it because re-uploading identical content
        under a new name must still refresh what the drill-down displays.

    Returns
    -------
    tuple or None
        The new entries tuple, or ``None`` when nothing changed and the
        caller should leave the reactive value alone.
    """
    current = next((e for e in entries if e['kind'] == kind), None)
    clearing = annotation is None or not len(annotation)
    if clearing:
        # Clearing a slot that is already empty is a no-op, not a reset.
        return (
            None if current is None else tuple(e for e in entries if e['kind'] != kind)
        )
    if current is not None and current.get('key') is not None and current['key'] == key:
        return None
    return tuple(
        [e for e in entries if e['kind'] != kind]
        + [
            {
                'kind': kind,
                'filename': filename,
                'annotation': annotation,
                'key': key,
            }
        ]
    )
