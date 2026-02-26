"""GFF feature annotation support for rusty-dot dotplots.

Provides the :class:`GffAnnotation` class for loading GFF features and
mapping them onto dotplot axes.

Examples
--------
>>> from rusty_dot.annotation import GffAnnotation
>>> ann = GffAnnotation.from_file("features.gff")
>>> ann.feature_types()
['CDS', 'gene', 'repeat_region']
>>> coding = ann.keep_feature_types(['CDS', 'gene'])
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Generator, Optional, Union

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

_log = logging.getLogger(__name__)

# Supported colour palette names (matplotlib qualitative colormaps).
SUPPORTED_PALETTES: list[str] = [
    'tab10',
    'tab20',
    'Set1',
    'Set2',
    'Set3',
    'Paired',
    'Dark2',
    'Accent',
]


@dataclass
class GffFeature:
    """A single feature record parsed from a GFF file.

    Coordinates are stored in 0-based half-open ``[start, end)`` format,
    converted from the 1-based inclusive GFF convention on load.

    Parameters
    ----------
    seqname : str
        Sequence (chromosome / contig) name.
    source : str
        Source field from column 2 of the GFF.
    feature_type : str
        Feature type (column 3), e.g. ``'gene'``, ``'CDS'``.
    start : int
        0-based start coordinate (inclusive).
    end : int
        0-based end coordinate (exclusive).
    score : float or None
        Alignment score, or ``None`` if the GFF field is ``'.'``.
    strand : str
        ``'+'``, ``'-'``, or ``'.'``.
    frame : int or None
        Reading frame (0, 1, or 2), or ``None`` if the GFF field is ``'.'``.
    attributes : str
        Raw attributes string from column 9 of the GFF.
    """

    seqname: str
    source: str
    feature_type: str
    start: int
    end: int
    score: Optional[float]
    strand: str
    frame: Optional[int]
    attributes: str


def _parse_gff(
    gff_path: Union[str, Path],
) -> Generator[GffFeature, None, None]:
    """Parse a GFF3 file and yield :class:`GffFeature` records.

    Comment lines (starting with ``#``) and lines with fewer than 8
    tab-separated fields are silently skipped.  Coordinates are converted
    from the 1-based inclusive GFF convention to 0-based half-open
    ``[start, end)`` format.

    Parameters
    ----------
    gff_path : str or Path
        Path to the GFF file.

    Yields
    ------
    GffFeature
        One record per valid GFF data line.

    Raises
    ------
    FileNotFoundError
        If *gff_path* does not exist.
    """
    with open(gff_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            seqname = parts[0]
            source = parts[1]
            feature_type = parts[2]
            start = int(parts[3]) - 1  # GFF is 1-based inclusive → 0-based
            end = int(parts[4])  # GFF end is 1-based inclusive → keep as exclusive
            score_str = parts[5]
            score = None if score_str == '.' else float(score_str)
            strand = parts[6]
            frame_str = parts[7]
            frame = None if frame_str == '.' else int(frame_str)
            attributes = parts[8] if len(parts) > 8 else ''
            yield GffFeature(
                seqname=seqname,
                source=source,
                feature_type=feature_type,
                start=start,
                end=end,
                score=score,
                strand=strand,
                frame=frame,
                attributes=attributes,
            )


class GffAnnotation:
    """Feature annotations loaded from a GFF file for dotplot overlays.

    Each unique feature type is automatically assigned a colour from a
    qualitative Matplotlib colour palette.  Colours can be overridden per
    feature type via the *colors* argument or the :meth:`set_colors` method.

    Parameters
    ----------
    records : list[GffFeature]
        Pre-parsed feature records.
    colors : dict[str, str], optional
        Mapping of ``feature_type → colour`` (any Matplotlib colour
        specification).  Feature types not listed here receive automatic
        palette colours.
    palette : str, optional
        Name of a Matplotlib qualitative palette used for automatic colour
        assignment.  Must be one of :data:`SUPPORTED_PALETTES`.
        Default is ``'tab10'``.

    Examples
    --------
    >>> from rusty_dot.annotation import GffAnnotation
    >>> ann = GffAnnotation.from_file("features.gff", palette="Set2")
    >>> ann.keep_feature_types(["gene"]).sequence_names()
    ['chr1', 'chr2']
    """

    def __init__(
        self,
        records: list[GffFeature],
        colors: Optional[dict[str, str]] = None,
        palette: str = 'tab10',
    ) -> None:
        """Initialise a GffAnnotation.

        Parameters
        ----------
        records : list[GffFeature]
            Pre-parsed feature records.
        colors : dict[str, str], optional
            Mapping of ``feature_type → colour``.
        palette : str, optional
            Palette name for automatic colour assignment. Default ``'tab10'``.
        """
        if palette not in SUPPORTED_PALETTES:
            raise ValueError(
                f"Unknown palette {palette!r}. "
                f"Choose from: {SUPPORTED_PALETTES}"
            )
        self._records: list[GffFeature] = list(records)
        self._palette = palette
        self._colors: dict[str, str] = {}
        if colors:
            self._colors.update(colors)
        self._assign_colors()

    # ------------------------------------------------------------------
    # Class-method constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_file(
        cls,
        gff_path: Union[str, Path],
        colors: Optional[dict[str, str]] = None,
        palette: str = 'tab10',
    ) -> 'GffAnnotation':
        """Load a GFF file and return a :class:`GffAnnotation`.

        Parameters
        ----------
        gff_path : str or Path
            Path to the GFF3 annotation file.
        colors : dict[str, str], optional
            Override colours for specific feature types.
        palette : str, optional
            Palette name for automatic colour assignment.
            Default is ``'tab10'``.

        Returns
        -------
        GffAnnotation
            Populated annotation object.
        """
        records = list(_parse_gff(gff_path))
        return cls(records=records, colors=colors, palette=palette)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _assign_colors(self) -> None:
        """Assign palette colours to feature types not already coloured."""
        types = self.feature_types()
        cmap = plt.get_cmap(self._palette)
        n = getattr(cmap, 'N', 256)
        for i, ft in enumerate(types):
            if ft not in self._colors:
                self._colors[ft] = mcolors.to_hex(cmap(i % n))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_colors(self, colors: dict[str, str]) -> None:
        """Override colours for specific feature types.

        Parameters
        ----------
        colors : dict[str, str]
            Mapping of ``feature_type → colour`` (any Matplotlib colour
            specification, e.g. ``'red'``, ``'#ff0000'``).
        """
        self._colors.update(colors)

    def feature_types(self) -> list[str]:
        """Return sorted list of unique feature types in the annotation.

        Returns
        -------
        list[str]
            Sorted unique feature type names.
        """
        return sorted({r.feature_type for r in self._records})

    def sequence_names(self) -> list[str]:
        """Return sorted list of unique sequence names in the annotation.

        Returns
        -------
        list[str]
            Sorted unique sequence names.
        """
        return sorted({r.seqname for r in self._records})

    def keep_feature_types(self, feature_types: list[str]) -> 'GffAnnotation':
        """Return a new annotation containing only the specified feature types.

        Parameters
        ----------
        feature_types : list[str]
            Feature type names to retain.

        Returns
        -------
        GffAnnotation
            Filtered annotation with the same colours and palette.
        """
        keep = set(feature_types)
        new_records = [r for r in self._records if r.feature_type in keep]
        return GffAnnotation(
            new_records,
            colors=dict(self._colors),
            palette=self._palette,
        )

    def drop_feature_types(self, feature_types: list[str]) -> 'GffAnnotation':
        """Return a new annotation with the specified feature types removed.

        Parameters
        ----------
        feature_types : list[str]
            Feature type names to remove.

        Returns
        -------
        GffAnnotation
            Filtered annotation with the same colours and palette.
        """
        drop = set(feature_types)
        new_records = [r for r in self._records if r.feature_type not in drop]
        return GffAnnotation(
            new_records,
            colors=dict(self._colors),
            palette=self._palette,
        )

    def filter_by_sequence(self, sequence_names: list[str]) -> 'GffAnnotation':
        """Return a new annotation containing only features from specified sequences.

        Parameters
        ----------
        sequence_names : list[str]
            Sequence names to retain.

        Returns
        -------
        GffAnnotation
            Filtered annotation with the same colours and palette.
        """
        keep = set(sequence_names)
        new_records = [r for r in self._records if r.seqname in keep]
        return GffAnnotation(
            new_records,
            colors=dict(self._colors),
            palette=self._palette,
        )

    def get_features_for_sequence(self, seq_name: str) -> list[GffFeature]:
        """Return all features belonging to a specific sequence.

        Parameters
        ----------
        seq_name : str
            Sequence name to look up.

        Returns
        -------
        list[GffFeature]
            Features for *seq_name* in the order they appear in the file.
        """
        return [r for r in self._records if r.seqname == seq_name]

    def get_color(self, feature_type: str) -> str:
        """Return the colour assigned to a feature type.

        Parameters
        ----------
        feature_type : str
            Feature type name.

        Returns
        -------
        str
            Hex colour string, or ``'#888888'`` if the type is unknown.
        """
        return self._colors.get(feature_type, '#888888')

    @property
    def records(self) -> list[GffFeature]:
        """All feature records held by this annotation.

        Returns
        -------
        list[GffFeature]
            Immutable view of the internal record list.
        """
        return list(self._records)

    def __len__(self) -> int:
        """Return the number of feature records."""
        return len(self._records)

    def __repr__(self) -> str:
        """Return a string representation."""
        return (
            f'GffAnnotation('
            f'{len(self._records)} records, '
            f'{len(self.feature_types())} types, '
            f'palette={self._palette!r})'
        )
