"""Tests for the dotplot visualization module."""

import os
import pytest

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.dotplot import DotPlotter


@pytest.fixture
def dotplot_index():
    """Create an index for dotplot testing."""
    idx = SequenceIndex(k=4)
    idx.add_sequence("seq1", "ACGTACGTACGTACGTACGT")
    idx.add_sequence("seq2", "TACGTACGTACGTACGTACG")
    idx.add_sequence("seq3", "GCGCGCGCGCGCGCGCGCGC")
    return idx


def test_dotplotter_creation(dotplot_index):
    """Test DotPlotter can be created."""
    plotter = DotPlotter(dotplot_index)
    assert plotter.index is dotplot_index


def test_plot_all_vs_all(dotplot_index, tmp_path):
    """Test all-vs-all dotplot generation."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / "dotplot.png")
    plotter.plot(output_path=output)
    assert os.path.exists(output)
    assert os.path.getsize(output) > 0


def test_plot_single(dotplot_index, tmp_path):
    """Test single pair dotplot generation."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / "single.png")
    plotter.plot_single("seq1", "seq2", output_path=output)
    assert os.path.exists(output)


def test_plot_subset(dotplot_index, tmp_path):
    """Test plotting a subset of sequences."""
    plotter = DotPlotter(dotplot_index)
    output = str(tmp_path / "subset.png")
    plotter.plot(
        query_names=["seq1", "seq2"],
        target_names=["seq2", "seq3"],
        output_path=output,
    )
    assert os.path.exists(output)


def test_plot_empty_index_raises():
    """Test that plotting with empty index raises ValueError."""
    idx = SequenceIndex(k=4)
    plotter = DotPlotter(idx)
    with pytest.raises(ValueError):
        plotter.plot()
