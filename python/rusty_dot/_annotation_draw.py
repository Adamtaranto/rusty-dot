"""Shared drawing helpers for GFF annotation tracks.

Used by :meth:`rusty_dot.dotplot.DotPlotter.plot` (focused single-pair
views) and :meth:`~rusty_dot.dotplot.DotPlotter.plot_single` to render
side annotation tracks: lane-packed feature shapes beside the dotplot
axes, with genome-browser conventions —

* unstranded features as rounded rectangles;
* stranded feature types (gene, mRNA, exon, CDS, ORF) as arrow polygons
  whose head points in the feature's direction (flipping when the axis is
  displayed reverse-complemented);
* multi-part groups (e.g. a multi-exon CDS sharing an ``ID``/``Parent``)
  on one lane, joined by a connector line through the part midpoints.

Lane assignment is a greedy interval scheduling over group extents, so
overlapping features stack instead of hiding each other.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Polygon, Rectangle

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .annotation import GffAnnotation, GffFeature

#: Feature types drawn as direction arrows (matched case-insensitively).
STRANDED_TYPES: frozenset[str] = frozenset({'gene', 'exon', 'mrna', 'cds', 'orf'})

#: Vertical padding inside a lane (fraction of lane height on each side).
_LANE_PAD = 0.12

#: Arrow-head length as a fraction of the feature length, capped at this
#: fraction of the sequence length so heads stay visually reasonable.
_HEAD_FRAC = 0.35
_HEAD_CAP_FRAC = 0.01


def assign_lanes(extents: list[tuple[int, int]], gap: int = 0) -> list[int]:
    """Assign non-overlapping display lanes to interval extents.

    Greedy interval scheduling: extents are processed in start order and
    placed in the first lane whose last occupant ends at least *gap* bases
    before the new start.  Input order is preserved in the returned list.

    Parameters
    ----------
    extents : list[tuple[int, int]]
        ``(start, end)`` half-open extents (e.g. whole-group spans).
    gap : int, optional
        Minimum bases required between neighbours in a lane.  Default ``0``.

    Returns
    -------
    list[int]
        Lane index (0-based) per input extent, in input order.
    """
    order = sorted(range(len(extents)), key=lambda i: extents[i][0])
    lane_ends: list[int] = []
    lanes = [0] * len(extents)
    for i in order:
        start, end = extents[i]
        for lane, lane_end in enumerate(lane_ends):
            if start >= lane_end + gap:
                lanes[i] = lane
                lane_ends[lane] = end
                break
        else:
            lanes[i] = len(lane_ends)
            lane_ends.append(end)
    return lanes


def _mirror(start: int, end: int, seq_len: int) -> tuple[int, int]:
    """Mirror a half-open interval for a reverse-displayed axis."""
    return seq_len - end, seq_len - start


def _is_stranded(feature: 'GffFeature') -> bool:
    """Whether a feature should be drawn as a direction arrow."""
    return feature.feature_type.lower() in STRANDED_TYPES and feature.strand in (
        '+',
        '-',
    )


def _xy(along: float, across: float, orientation: str) -> tuple[float, float]:
    """Map (position-along-sequence, position-across-lanes) to axes x/y."""
    if orientation == 'x':
        return along, across
    return across, along


def _arrow_points(
    start: float,
    end: float,
    lo: float,
    hi: float,
    forward: bool,
    head_len: float,
    orientation: str,
) -> list[tuple[float, float]]:
    """Build a 5-point arrow polygon for a stranded feature.

    Parameters
    ----------
    start, end : float
        Feature interval along the sequence axis (display coordinates).
    lo, hi : float
        Lane extent across the track.
    forward : bool
        Arrow head at *end* when ``True``, at *start* when ``False``.
    head_len : float
        Head length along the sequence axis.
    orientation : str
        ``'x'`` (track under the x axis) or ``'y'`` (track beside the y
        axis).

    Returns
    -------
    list[tuple[float, float]]
        Polygon vertices in axes coordinates.
    """
    mid = (lo + hi) / 2.0
    if forward:
        base = max(start, end - head_len)
        pts = [(start, lo), (base, lo), (end, mid), (base, hi), (start, hi)]
    else:
        base = min(end, start + head_len)
        pts = [(end, lo), (base, lo), (start, mid), (base, hi), (end, hi)]
    return [_xy(a, c, orientation) for a, c in pts]


def draw_track(
    ax: Axes,
    annotation: 'GffAnnotation',
    seq_name: str,
    seq_len: int,
    orientation: str,
    reverse: bool = False,
    max_lanes: int = 6,
    gap_bp: int = 0,
) -> int:
    """Draw a lane-packed annotation track on a side axis.

    The axis along the sequence must already share limits with the main
    dotplot panel (``sharex``/``sharey``); this function only sets the
    across-track limits.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The track axes.
    annotation : GffAnnotation
        Features and colours.
    seq_name : str
        Sequence whose features are drawn.
    seq_len : int
        Sequence length in bases.
    orientation : str
        ``'x'`` for the track below the x axis (features run left-right),
        ``'y'`` for the track beside the y axis (features run up-down).
    reverse : bool, optional
        Mirror coordinates and flip strands for an axis displayed
        reverse-complemented.  Default ``False``.
    max_lanes : int, optional
        Lanes beyond this are clamped onto the last lane to bound track
        height.  Default ``6``.
    gap_bp : int, optional
        Minimum bases between lane neighbours.  Default ``0``.

    Returns
    -------
    int
        Number of lanes used (0 when the sequence has no features).
    """
    groups = annotation.iter_groups(seq_name)
    if not groups:
        ax.set_axis_off()
        return 0

    extents = []
    for _, parts in groups:
        lo = min(p.start for p in parts)
        hi = max(p.end for p in parts)
        if reverse:
            lo, hi = _mirror(lo, hi, seq_len)
        extents.append((lo, hi))
    lanes = [min(lane, max_lanes - 1) for lane in assign_lanes(extents, gap_bp)]
    n_lanes = max(lanes) + 1

    head_cap = seq_len * _HEAD_CAP_FRAC
    for ((_gid, _parent), parts), lane in zip(groups, lanes):
        lane_lo = lane + _LANE_PAD
        lane_hi = lane + 1 - _LANE_PAD
        color = annotation.get_color(parts[0].feature_type)
        mids: list[float] = []
        for feat in parts:
            start, end = float(feat.start), float(feat.end)
            strand = feat.strand
            if reverse:
                start, end = _mirror(feat.start, feat.end, seq_len)
                strand = {'+': '-', '-': '+'}.get(strand, strand)
            mids.append((start + end) / 2.0)
            if _is_stranded(feat):
                head_len = min(_HEAD_FRAC * (end - start), head_cap)
                ax.add_patch(
                    Polygon(
                        _arrow_points(
                            start,
                            end,
                            lane_lo,
                            lane_hi,
                            forward=(strand == '+'),
                            head_len=head_len,
                            orientation=orientation,
                        ),
                        closed=True,
                        facecolor=color,
                        edgecolor='none',
                        zorder=2,
                    )
                )
            else:
                length = end - start
                rounding = min(0.25 * length, 0.005 * seq_len)
                xy = _xy(start, lane_lo, orientation)
                if orientation == 'x':
                    w, h = length, lane_hi - lane_lo
                else:
                    w, h = lane_hi - lane_lo, length
                    xy = _xy(start, lane_lo, orientation)
                if rounding <= 0:
                    ax.add_patch(
                        Rectangle(xy, w, h, facecolor=color, edgecolor='none', zorder=2)
                    )
                else:
                    # mutation_aspect rescales the across-track rounding so the
                    # corners look round despite bp-vs-lane axis anisotropy.
                    aspect = (0.25 * (1 - 2 * _LANE_PAD)) / rounding
                    if orientation == 'y':
                        aspect = 1.0 / aspect if aspect else 1.0
                    ax.add_patch(
                        FancyBboxPatch(
                            xy,
                            w,
                            h,
                            boxstyle=f'round,pad=0,rounding_size={rounding}',
                            mutation_aspect=aspect,
                            facecolor=color,
                            edgecolor='none',
                            zorder=2,
                        )
                    )
        if len(parts) > 1:
            # Connector through the part midpoints at lane centre.
            centre = lane + 0.5
            xs, ys = zip(*(_xy(m, centre, orientation) for m in sorted(mids)))
            ax.add_line(Line2D(xs, ys, color=color, linewidth=0.8, zorder=1))

    if orientation == 'x':
        ax.set_ylim(0, n_lanes)
        ax.invert_yaxis()  # lane 0 nearest the main panel
    else:
        # Lane 0 nearest the main panel: the track sits left of the y axis,
        # so lanes grow leftwards.
        ax.set_xlim(n_lanes, 0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(
        left=False,
        right=False,
        top=False,
        bottom=False,
        labelleft=False,
        labelright=False,
        labeltop=False,
        labelbottom=False,
    )
    return n_lanes


def annotation_legend_handles(
    annotation: 'GffAnnotation',
    feature_types: Optional[list[str]] = None,
) -> list[Patch]:
    """Build proxy legend handles for an annotation's feature types.

    Parameters
    ----------
    annotation : GffAnnotation
        Colour source.
    feature_types : list[str], optional
        Subset/order of types to include.  Defaults to every type present.

    Returns
    -------
    list[matplotlib.patches.Patch]
        One coloured patch per feature type, labelled with the type name.
    """
    types = feature_types if feature_types is not None else annotation.feature_types()
    return [
        Patch(facecolor=annotation.get_color(ft), edgecolor='none', label=ft)
        for ft in types
    ]
