"""Tests for the Nature-journal-style plot formatting module."""

import matplotlib
import pytest

from rusty_dot import (
    NATURE_RC,
    DotPlotter,
    SequenceIndex,
    nature_style,
    use_nature_style,
)


@pytest.fixture
def restore_rcparams():
    """Snapshot matplotlib rcParams and restore them after the test."""
    saved = matplotlib.rcParams.copy()
    yield
    matplotlib.rcParams.update(saved)


def _assert_nature_rcparams_active() -> None:
    """Assert that key NATURE_RC values are active in matplotlib.rcParams."""
    assert matplotlib.rcParams['font.family'] == ['sans-serif']
    assert matplotlib.rcParams['font.sans-serif'][:2] == ['Helvetica', 'Arial']
    assert matplotlib.rcParams['font.size'] == 7.0
    assert matplotlib.rcParams['xtick.labelsize'] == 6.0
    assert matplotlib.rcParams['axes.spines.top'] is False
    assert matplotlib.rcParams['axes.spines.right'] is False
    assert matplotlib.rcParams['axes.linewidth'] == 0.5
    assert matplotlib.rcParams['lines.linewidth'] == 0.5
    assert matplotlib.rcParams['xtick.major.width'] == 0.5
    assert matplotlib.rcParams['xtick.direction'] == 'out'
    assert matplotlib.rcParams['savefig.dpi'] == 300
    assert matplotlib.rcParams['savefig.bbox'] == 'tight'


class TestNatureRc:
    """Tests for the NATURE_RC constant."""

    def test_keys_are_valid_rcparams(self):
        """Every key in NATURE_RC must be a valid matplotlib rcParam."""
        for key in NATURE_RC:
            assert key in matplotlib.rcParams, f'invalid rcParam key: {key}'

    def test_font_sizes_in_nature_range(self):
        """All font sizes fall in the 5-7 pt Nature range."""
        size_keys = [k for k in NATURE_RC if k.endswith('size')]
        for key in size_keys:
            assert 5.0 <= NATURE_RC[key] <= 7.0, key


class TestNatureStyleContextManager:
    """Tests for the nature_style context manager."""

    def test_applies_inside_context(self):
        """Key rcParams take Nature values inside the with block."""
        with nature_style():
            _assert_nature_rcparams_active()

    def test_restores_rcparams_on_exit(self):
        """rcParams are restored to their prior values after exit."""
        before = dict(matplotlib.rcParams)
        with nature_style():
            pass
        after = dict(matplotlib.rcParams)
        for key in NATURE_RC:
            assert after[key] == before[key], key

    def test_restores_rcparams_on_exception(self):
        """rcParams are restored even when the block raises."""
        before_top = matplotlib.rcParams['axes.spines.top']
        before_size = matplotlib.rcParams['font.size']
        with pytest.raises(RuntimeError):
            with nature_style():
                raise RuntimeError('boom')
        assert matplotlib.rcParams['axes.spines.top'] == before_top
        assert matplotlib.rcParams['font.size'] == before_size


class TestUseNatureStyle:
    """Tests for the global use_nature_style function."""

    def test_applies_globally(self, restore_rcparams):
        """use_nature_style() mutates rcParams outside any context."""
        use_nature_style()
        _assert_nature_rcparams_active()


class TestStyledDotPlotterFigure:
    """Render a small DotPlotter figure under the style."""

    @pytest.fixture
    def index(self):
        """Build a small SequenceIndex with two synthetic contigs."""
        idx = SequenceIndex(k=10)
        idx.add_sequence('seq1', 'ACGTACGTACGTACGTACGT')
        idx.add_sequence('seq2', 'TACGTACGTACGTACGTACG')
        return idx

    def test_figure_reflects_style(self, index):
        """Figure produced inside nature_style() has styled axes."""
        import matplotlib.pyplot as plt

        with nature_style():
            fig = DotPlotter(index).plot(query_names=['seq1'], target_names=['seq2'])
        try:
            axes = fig.get_axes()
            assert axes, 'expected at least one axes in the figure'
            for ax in axes:
                assert not ax.spines['top'].get_visible()
                assert not ax.spines['right'].get_visible()
                assert ax.spines['left'].get_linewidth() == 0.5
                for tick in ax.get_xticklabels() + ax.get_yticklabels():
                    assert tick.get_fontsize() == 6.0
        finally:
            plt.close(fig)

    def test_default_figure_differs(self, index):
        """Without the style, default figures keep top/right spines."""
        import matplotlib.pyplot as plt

        fig = DotPlotter(index).plot(query_names=['seq1'], target_names=['seq2'])
        try:
            ax = fig.get_axes()[0]
            assert ax.spines['top'].get_visible()
            assert ax.spines['right'].get_visible()
        finally:
            plt.close(fig)
