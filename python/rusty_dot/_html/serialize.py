"""Payload serialisation for the interactive HTML dotplot report.

The :class:`~rusty_dot.dotplot.DotPlotter` grid renderer captures, for every
panel it draws, the exact segments handed to matplotlib (already filtered,
mirrored and chained).  :func:`build_panel_payload` converts that raw capture
into the compact JSON-ready dictionary embedded in the HTML report.

Payload contract (consumed by ``report.js``)
--------------------------------------------
::

    {
      "panels": {
        "rd-panel-<row>-<col>": {
          "query":  <display query name>,
          "target": <display target name>,
          "qlen":   <query length>,
          "tlen":   <target length>,
          "segments": {
            "fwd":      [[qs, qe, ts, te], ...],
            "rev":      [[qs, qe, ts, te], ...],
            "identity": [[qs, qe, ts, te, identity, strand], ...]
          },
          "seqs": {"fwd": [...], "rev": [...], "identity": [...]}  # optional
        },
        ...
      },
      "has_sequences": true|false,
      "tracks": {                      # optional, 1x1 track layouts only
        "x": [{"gid": "rd-xtrack-<n>", "group": <int>, "type": ...,
               "seqname": ..., "start": ..., "end": ..., "strand": ...,
               "id": ..., "parent": ..., "name": ..., "source": ...,
               "source_file": ..., "color": ...}, ...],
        "y": [...]
      }
    }

``tracks`` entries are in draw order and carry their own SVG gid, so the
report addresses them directly rather than by index.  ``group`` is shared
by the parts of one multi-part feature (a spliced CDS), which is what lets
a modifier-click highlight the whole feature at once.  Side tracks only
exist on single-pair layouts, so the key sits at the top level rather than
under a panel.

The order of entries in each ``segments`` list matches the drawing order of
the corresponding :class:`~matplotlib.collections.LineCollection`, which in
turn matches the document order of the per-segment SVG elements inside the
``<g id="rd-matches-<row>-<col>-<layer>">`` group.  That ordering is the
contract that lets ``report.js`` map a clicked SVG element back to its
coordinates.
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Optional

_log = logging.getLogger(__name__)

# Soft cap on the total number of embedded residues (~2 MB of sequence text).
MAX_EMBEDDED_RESIDUES = 2_000_000

# Segment layers whose entries carry a query-side sequence slice.
_LAYERS = ('fwd', 'rev', 'identity')

# Complement table for reverse-complementing embedded sequence slices.
# Unmapped characters (ambiguity codes other than N) pass through unchanged.
_COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def _revcomp(seq: str) -> str:
    """Reverse-complement a DNA string (N and case preserved).

    Parameters
    ----------
    seq : str
        DNA sequence; characters outside ``ACGTNacgtn`` are reversed but not
        complemented.

    Returns
    -------
    str
        The reverse complement.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def _total_query_residues(panels: dict[str, dict[str, Any]]) -> int:
    """Sum the query-side span of every captured segment.

    Parameters
    ----------
    panels : dict
        Raw per-panel capture as produced by
        :meth:`~rusty_dot.dotplot.DotPlotter.plot` (``segments`` lists whose
        first two entries are the query start/end).

    Returns
    -------
    int
        Total number of query residues that embedding sequences would add
        to the payload.
    """
    total = 0
    for panel in panels.values():
        for layer in _LAYERS:
            for seg in panel['segments'][layer]:
                total += int(seg[1]) - int(seg[0])
    return total


def _panel_sequences(
    panel: dict[str, Any],
    query_seq: str,
) -> dict[str, list[str]]:
    """Slice the query sequence for every segment of one panel.

    When the panel was rendered with a reverse-oriented query
    (``reverse_query`` set by :meth:`~rusty_dot.dotplot.DotPlotter.plot`
    via *reverse_contigs*), the captured coordinates are mirrored
    (``q' = qlen - q``) relative to the stored forward-orientation sequence.
    The substring is then taken from the original (un-mirrored) region and
    reverse-complemented so it reads exactly as displayed on the plot.

    Parameters
    ----------
    panel : dict
        Raw capture entry for one panel (``segments`` lists of
        ``[qs, qe, ...]`` rows, plus ``qlen`` and ``reverse_query``).
    query_seq : str
        Full query sequence for this panel's row (forward orientation, as
        stored in the index).

    Returns
    -------
    dict
        ``{layer: [substring, ...]}`` with one query-side substring per
        segment, in segment order and display orientation.
    """
    reverse = bool(panel.get('reverse_query', False))
    qlen = int(panel['qlen'])
    seqs: dict[str, list[str]] = {}
    for layer in _LAYERS:
        slices = []
        for seg in panel['segments'][layer]:
            qs, qe = int(seg[0]), int(seg[1])
            if reverse:
                # Mirrored display coords -> original region, then revcomp
                # so the fragment matches the plotted orientation.
                slices.append(_revcomp(query_seq[qlen - qe : qlen - qs]))
            else:
                slices.append(query_seq[qs:qe])
        seqs[layer] = slices
    return seqs


def build_panel_payload(
    capture: dict[str, Any],
    get_sequence: Optional[Callable[[str], str]] = None,
    max_residues: int = MAX_EMBEDDED_RESIDUES,
) -> dict[str, Any]:
    """Build the JSON-ready payload embedded in the HTML report.

    Pure function: reads the raw capture assembled during
    :meth:`~rusty_dot.dotplot.DotPlotter.plot` and returns a new dictionary;
    the capture itself is not modified.

    Parameters
    ----------
    capture : dict
        Raw capture with a ``'panels'`` mapping of
        ``gid -> {query, target, query_id, qlen, tlen, segments}``.
        ``query_id`` is the internal (possibly group-prefixed) name used for
        sequence lookup and is not propagated to the payload.
    get_sequence : callable or None, optional
        ``name -> sequence`` accessor (duck-typed from the plotter's index).
        When ``None`` (e.g. the index is a PAF alignment with no stored
        sequences) no sequences are embedded.  Default is ``None``.
    max_residues : int, optional
        Upper bound on the total number of embedded residues.  When the
        capture's segments span more query residues than this, sequences are
        omitted and ``has_sequences`` is ``False``.  Default is
        :data:`MAX_EMBEDDED_RESIDUES` (~2 MB).

    Returns
    -------
    dict
        Payload dictionary with ``'panels'`` and ``'has_sequences'`` keys,
        ready for compact ``json.dumps``.
    """
    raw_panels: dict[str, dict[str, Any]] = capture['panels']

    embed = get_sequence is not None
    if embed:
        total = _total_query_residues(raw_panels)
        if total > max_residues:
            _log.info(
                'HTML report: %d match residues exceed the %d embedding cap; '
                'sequences omitted from the payload.',
                total,
                max_residues,
            )
            embed = False

    # Cache full query sequences per internal name: one row of panels shares
    # a query, so fetch each sequence at most once.
    seq_cache: dict[str, Optional[str]] = {}

    panels_out: dict[str, dict[str, Any]] = {}
    has_sequences = False
    for gid, panel in raw_panels.items():
        entry: dict[str, Any] = {
            'query': panel['query'],
            'target': panel['target'],
            'qlen': int(panel['qlen']),
            'tlen': int(panel['tlen']),
            'segments': {
                layer: [list(seg) for seg in panel['segments'][layer]]
                for layer in _LAYERS
            },
        }
        # Diagonal-panel GFF features: one dict per drawn patch, in draw
        # order — the report JS maps the SVG children of the panel's
        # 'rd-annot-<r>-<c>' group back to these entries by index.
        if panel.get('annotations'):
            entry['annotations'] = [dict(a) for a in panel['annotations']]
        if embed:
            qid = panel['query_id']
            if qid not in seq_cache:
                try:
                    seq_cache[qid] = get_sequence(qid)  # type: ignore[misc]
                except Exception:  # noqa: BLE001 - index without this seq
                    _log.warning(
                        'HTML report: could not fetch sequence for %r; '
                        'omitting sequences for its panels.',
                        qid,
                    )
                    seq_cache[qid] = None
            query_seq = seq_cache[qid]
            if query_seq is not None:
                entry['seqs'] = _panel_sequences(panel, query_seq)
                has_sequences = True
        panels_out[gid] = entry

    payload: dict[str, Any] = {
        'panels': panels_out,
        'has_sequences': has_sequences,
    }
    tracks = capture.get('tracks')
    if tracks and (tracks.get('x') or tracks.get('y')):
        payload['tracks'] = {
            axis: [_track_entry(axis, n, group, feat) for n, group, feat in entries]
            for axis, entries in tracks.items()
        }
    return payload


def _track_entry(axis: str, n: int, group: int, feat: Any) -> dict[str, Any]:
    """Describe one drawn side-track part for the report payload."""
    return {
        'gid': f'rd-{axis}track-{n}',
        'group': int(group),
        'type': feat.feature_type,
        'seqname': feat.seqname,
        'start': int(feat.start),
        'end': int(feat.end),
        'strand': feat.strand,
        'id': feat.feature_id,
        'parent': feat.parent,
        'name': feat.name,
        'source': feat.source,
        'source_file': getattr(feat, 'source_file', '') or '',
        'color': getattr(feat, 'color', None),
    }
