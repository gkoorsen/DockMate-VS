"""Results rendering selects the dedicated page without changing the workflow."""

from types import SimpleNamespace

import pytest

from dockmate_vs.gui.app import DockMateVSApp


@pytest.mark.parametrize("protocol", [False, True])
def test_rendering_selects_results_and_summary(monkeypatch, tmp_path, protocol):
    app = object.__new__(DockMateVSApp)
    app.results_tab = "results"
    app.results_summary_tab = "summary"
    app.results_charts_tab = "charts"
    selected = []
    app.workflow_notebook = SimpleNamespace(select=lambda tab: selected.append(("main", tab)))
    app.results_notebook = SimpleNamespace(select=lambda tab: selected.append(("results", tab)))
    app.after_idle = lambda callback: None
    for method in ("_clear_frame", "_populate_summary_tab", "_populate_charts_tab",
                   "_populate_protocol_report", "_populate_protocol_charts"):
        monkeypatch.setattr(app, method, lambda *args: None)
    if protocol:
        app._render_protocol_results(tmp_path / "results.csv", tmp_path / "summary.md")
    else:
        app._render_results({}, [])
    assert selected == [("results", "summary"), ("main", "results")]


@pytest.mark.parametrize("selected, expected_mode", [
    ("results", None), ("poses", None), ("screening", "screening"),
    ("protocol", "protocol_development"),
])
def test_view_navigation_preserves_workflow_mode(selected, expected_mode):
    app = object.__new__(DockMateVSApp)
    modes = []
    scroll_positions = []
    app._main_canvas = SimpleNamespace(yview_moveto=scroll_positions.append)
    app.workflow_notebook = SimpleNamespace(select=lambda: selected)
    app.results_tab, app.pose_viewer_tab = "results", "poses"
    app.screening_tab, app.protocol_tab = "screening", "protocol"
    app.mode_var = SimpleNamespace(set=modes.append)
    app._update_mode = lambda: None
    app._update_pair_count = lambda: None
    app.after_idle = lambda callback: None
    app.last_results_path = None
    app._on_workflow_changed()
    assert modes == ([] if expected_mode is None else [expected_mode])
    assert scroll_positions == [0]
