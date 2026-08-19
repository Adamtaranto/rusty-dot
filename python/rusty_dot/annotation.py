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

from dataclasses import dataclass, field
import gzip
import logging
from pathlib import Path
from typing import Generator, Iterable, Optional, Union
from urllib.parse import unquote

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

_log = logging.getLogger(__name__)

#: Magic bytes identifying gzip-compressed input.
_GZIP_MAGIC = b'\x1f\x8b'


def parse_attributes(attributes: str) -> dict[str, str]:
    """Parse a GFF3 column-9 attributes string into a dict.

    Splits on ``;`` then on the first ``=`` per field, percent-decoding
    values per the GFF3 specification.  Empty fields and trailing
    separators are tolerated; a field without ``=`` (GFF2-style
    ``key "value"``) falls back to splitting on the first space.

    Parameters
    ----------
    attributes : str
        Raw column-9 text (e.g. ``'ID=gene1;Name=abc%3B1'``).

    Returns
    -------
    dict[str, str]
        Attribute key → decoded value.  Later duplicates of a key win.
    """
    result: dict[str, str] = {}
    for chunk in attributes.split(';'):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '=' in chunk:
            key, _, value = chunk.partition('=')
        elif ' ' in chunk:
            key, _, value = chunk.partition(' ')
            value = value.strip('"')
        else:
            key, value = chunk, ''
        result[key.strip()] = unquote(value.strip())
    return result


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
    #: Lazily populated cache for :meth:`attr_dict` (not part of equality).
    _attr_cache: Optional[dict[str, str]] = field(
        default=None, repr=False, compare=False
    )

    def attr_dict(self) -> dict[str, str]:
        """Return the parsed column-9 attributes (cached).

        Returns
        -------
        dict[str, str]
            Attribute key → percent-decoded value.
        """
        if self._attr_cache is None:
            self._attr_cache = parse_attributes(self.attributes)
        return self._attr_cache

    @property
    def feature_id(self) -> Optional[str]:
        """The GFF3 ``ID`` attribute, or ``None`` when absent.

        Returns
        -------
        str or None
            Feature identifier.
        """
        return self.attr_dict().get('ID')

    @property
    def parent(self) -> Optional[str]:
        """The GFF3 ``Parent`` attribute (first value of a comma list).

        Returns
        -------
        str or None
            Parent feature identifier, or ``None`` when absent.
        """
        parent = self.attr_dict().get('Parent')
        if parent is None:
            return None
        return parent.split(',', 1)[0]

    @property
    def name(self) -> Optional[str]:
        """The GFF3 ``Name`` attribute, falling back to ``ID``.

        Returns
        -------
        str or None
            Display name for the feature.
        """
        attrs = self.attr_dict()
        return attrs.get('Name') or attrs.get('ID')


def _parse_gff_lines(lines: Iterable[str]) -> Generator[GffFeature, None, None]:
    """Parse GFF3 lines and yield :class:`GffFeature` records.

    Comment/directive lines (starting with ``#``) and lines with fewer than
    8 tab-separated fields are silently skipped; an embedded ``##FASTA``
    section terminates parsing.  Coordinates are converted from the 1-based
    inclusive GFF convention to 0-based half-open ``[start, end)`` format.

    Parameters
    ----------
    lines : Iterable[str]
        GFF text lines (with or without trailing newlines).

    Yields
    ------
    GffFeature
        One record per valid GFF data line.
    """
    for line in lines:
        line = line.rstrip('\n').rstrip('\r')
        if line.startswith('##FASTA'):
            break
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


def _parse_gff(
    gff_path: Union[str, Path],
) -> Generator[GffFeature, None, None]:
    """Parse a (optionally gzip-compressed) GFF3 file.

    Parameters
    ----------
    gff_path : str or Path
        Path to the GFF file (``.gff``/``.gff3``, optionally ``.gz``).

    Yields
    ------
    GffFeature
        One record per valid GFF data line.

    Raises
    ------
    FileNotFoundError
        If *gff_path* does not exist.
    """
    path = Path(gff_path)
    with open(path, 'rb') as probe:
        is_gzip = probe.read(2) == _GZIP_MAGIC
    opener = gzip.open if is_gzip else open
    with opener(path, 'rt') as fh:
        yield from _parse_gff_lines(fh)


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
                f'Unknown palette {palette!r}. Choose from: {SUPPORTED_PALETTES}'
            )
        self._records: list[GffFeature] = list(records)
        self._palette = palette
        self._colors: dict[str, str] = {}
        # Lazy per-sequence index (records are immutable after construction);
        # avoids a full scan per panel when plotting many-contig grids.
        self._by_seq: Optional[dict[str, list[GffFeature]]] = None
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

    @classmethod
    def from_text(
        cls,
        text: str,
        colors: Optional[dict[str, str]] = None,
        palette: str = 'tab10',
    ) -> 'GffAnnotation':
        """Build a :class:`GffAnnotation` from GFF3 text.

        Parameters
        ----------
        text : str
            GFF3 content.
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
        records = list(_parse_gff_lines(text.splitlines()))
        return cls(records=records, colors=colors, palette=palette)

    @classmethod
    def from_bytes(
        cls,
        data: bytes,
        colors: Optional[dict[str, str]] = None,
        palette: str = 'tab10',
    ) -> 'GffAnnotation':
        """Build a :class:`GffAnnotation` from raw (optionally gzipped) bytes.

        Intended for in-memory uploads (e.g. the browser app) where no file
        path exists.  Gzip input is detected by magic bytes.

        Parameters
        ----------
        data : bytes
            Raw ``.gff``/``.gff3`` file content, optionally gzip-compressed.
        colors : dict[str, str], optional
            Override colours for specific feature types.
        palette : str, optional
            Palette name for automatic colour assignment.
            Default is ``'tab10'``.

        Returns
        -------
        GffAnnotation
            Populated annotation object.

        Raises
        ------
        ValueError
            If gzip decompression fails or the content is not UTF-8 text.
        """
        if data[:2] == _GZIP_MAGIC:
            try:
                data = gzip.decompress(data)
            except (OSError, EOFError) as exc:
                raise ValueError(f'Could not decompress gzip input: {exc}') from exc
        try:
            text = data.decode('utf-8')
        except UnicodeDecodeError as exc:
            raise ValueError('Input is not UTF-8 text — is this a GFF3 file?') from exc
        return cls.from_text(text, colors=colors, palette=palette)

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
        if self._by_seq is None:
            by_seq: dict[str, list[GffFeature]] = {}
            for r in self._records:
                by_seq.setdefault(r.seqname, []).append(r)
            self._by_seq = by_seq
        return list(self._by_seq.get(seq_name, []))

    def iter_groups(
        self, seq_name: str
    ) -> list[tuple[tuple[Optional[str], Optional[str]], list[GffFeature]]]:
        """Group a sequence's features into multi-part units.

        Features sharing both an ``ID`` and a ``Parent`` attribute form one
        group (e.g. the parts of a multi-exon CDS, which GFF3 records as
        several lines with the same ``ID``); features without an ``ID`` are
        singleton groups.  Parts within a group are ordered by start
        coordinate.

        Parameters
        ----------
        seq_name : str
            Sequence name to look up.

        Returns
        -------
        list[tuple[tuple[str or None, str or None], list[GffFeature]]]
            ``((feature_id, parent), parts)`` pairs in first-appearance
            order.
        """
        groups: dict[object, list[GffFeature]] = {}
        for feat in self.get_features_for_sequence(seq_name):
            if feat.feature_id is not None and feat.parent is not None:
                key: object = (feat.feature_type, feat.feature_id, feat.parent)
            else:
                key = id(feat)  # singleton group
            groups.setdefault(key, []).append(feat)
        result = []
        for parts in groups.values():
            parts = sorted(parts, key=lambda f: f.start)
            head = parts[0]
            result.append(((head.feature_id, head.parent), parts))
        return result

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
