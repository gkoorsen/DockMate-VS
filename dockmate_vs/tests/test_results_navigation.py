"""Results rendering selects the dedicated page without changing the workflow."""

from types import SimpleNamespace
from pathlib import Path

import pandas as pd
import pytest

import dockmate_vs.gui.app as app_module
from dockmate_vs.gui.app import DockMateVSApp


class FakeResultsTable:
    def __init__(self):
        self.columns = ("compound", "score", "note")
        self.headings = {
            "compound": "Compound",
            "score": "Docking score",
            "note": "Note",
        }
        self.rows = {
            "row1": ("Ligand A", "-7.5", "primary"),
            "row2": ("Ligand B", "-6.8", "contains\ttab"),
        }
        self.clipboard = ""

    def __getitem__(self, key):
        assert key == "columns"
        return self.columns

    def heading(self, column, option):
        assert option == "text"
        return self.headings[column]

    def get_children(self, parent=""):
        assert parent == ""
        return tuple(self.rows)

    def item(self, item_id, option):
        assert option == "values"
        return self.rows[item_id]

    def clipboard_clear(self):
        self.clipboard = ""

    def clipboard_append(self, value):
        self.clipboard += value

    def winfo_toplevel(self):
        return None


def test_results_table_copy_includes_displayed_headings_and_rows():
    table = FakeResultsTable()
    statuses = []

    DockMateVSApp._copy_treeview_table(table, SimpleNamespace(set=statuses.append))

    assert table.clipboard == (
        "Compound\tDocking score\tNote\n"
        "Ligand A\t-7.5\tprimary\n"
        "Ligand B\t-6.8\t\"contains\ttab\"\n"
    )
    assert statuses == ["Copied 2 row(s)"]


def test_results_table_exports_displayed_data_to_excel(monkeypatch, tmp_path):
    table = FakeResultsTable()
    output_path = tmp_path / "pose_recovery.xlsx"
    statuses = []
    messages = []
    monkeypatch.setattr(
        app_module.filedialog,
        "asksaveasfilename",
        lambda **_kwargs: str(output_path),
    )
    monkeypatch.setattr(
        app_module.messagebox,
        "showinfo",
        lambda *args, **kwargs: messages.append((args, kwargs)),
    )

    exported = DockMateVSApp._export_treeview_table(
        table,
        "Protocol pose recovery",
        SimpleNamespace(set=statuses.append),
    )

    frame = pd.read_excel(output_path, dtype=str)
    assert exported == output_path
    assert list(frame.columns) == ["Compound", "Docking score", "Note"]
    assert frame.to_dict("records") == [
        {"Compound": "Ligand A", "Docking score": "-7.5", "Note": "primary"},
        {
            "Compound": "Ligand B",
            "Docking score": "-6.8",
            "Note": "contains\ttab",
        },
    ]
    assert statuses == ["Exported 2 row(s)"]
    assert len(messages) == 1


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


def test_results_scroll_target_uses_innermost_registered_canvas():
    app = object.__new__(DockMateVSApp)
    outer_canvas = SimpleNamespace(master=None)
    outer_content = SimpleNamespace(master=outer_canvas)
    inner_canvas = SimpleNamespace(master=outer_content)
    inner_content = SimpleNamespace(master=inner_canvas)
    leaf = SimpleNamespace(master=inner_content)
    app._results_scroll_canvases = [outer_canvas, inner_canvas]

    assert app._results_scroll_target(leaf) is inner_canvas
    assert app._results_scroll_target(outer_content) is outer_canvas
    assert app._results_scroll_target(SimpleNamespace(master=None)) is None


def test_pose_output_file_resolves_after_run_folder_is_copied(tmp_path):
    run_dir = tmp_path / "anathi_template_updated_list_30_matched_decoys"
    pose_file = (
        run_dir
        / "Mpro_ligand_001"
        / "variants"
        / "ligand_001_v1"
        / "docked.pdbqt"
    )
    pose_file.parent.mkdir(parents=True)
    pose_file.write_text("MODEL 1\nENDMDL\n")

    stale_path = (
        "/old/output/anathi_template_updated_list_30_matched_decoys/"
        "Mpro_ligand_001/variants/ligand_001_v1/docked.pdbqt"
    )

    resolved = DockMateVSApp._resolve_pose_output_file(
        stale_path,
        run_dir / "redock_results.csv",
    )

    assert resolved == pose_file


def test_pose_output_file_resolves_from_protocol_results_folder(tmp_path):
    protocol_dir = (
        tmp_path / "anathi_template_updated_list_30_matched_decoys" / "protocol_development"
    )
    pose_file = (
        protocol_dir
        / "Mpro_protocol_case"
        / "variants"
        / "ligand_001_v1"
        / "docked.pdbqt"
    )
    pose_file.parent.mkdir(parents=True)
    pose_file.write_text("MODEL 1\nENDMDL\n")

    stale_path = (
        "/old/output/anathi_template_updated_list_30_matched_decoys/protocol_development/"
        "Mpro_protocol_case/variants/ligand_001_v1/docked.pdbqt"
    )

    resolved = DockMateVSApp._resolve_pose_output_file(
        stale_path,
        protocol_dir / "protocol_development_results.csv",
    )

    assert resolved == pose_file


def test_pose_viewer_case_label_distinguishes_complex_and_variant():
    label = DockMateVSApp._pose_case_label(
        {
            "pdb_id": "5REE",
            "ligand": "T1M",
            "display_name": "Gancaonin P",
            "output_file": Path(
                "5REE_T1M_Gancaonin_P/variants/Gancaonin P_v1/docked.pdbqt"
            ),
        },
        6,
    )

    assert label == "0007 | 5REE/T1M | Gancaonin P | variant Gancaonin P_v1"
    assert DockMateVSApp._pose_label_key(label) == (
        "0007_|_5REE/T1M_|_GANCAONIN_P_|_VARIANT_GANCAONIN_P_V1"
    )
