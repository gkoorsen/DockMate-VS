from pathlib import Path

import dockmate_vs.gui.app as app_module
from dockmate_vs.gui.app import DockMateVSApp
from dockmate_vs.gui.folder_picker import FolderPickerDialog, subdirectories


def test_subdirectories_returns_only_folders_in_name_order(tmp_path: Path):
    (tmp_path / "zeta").mkdir()
    (tmp_path / "Alpha").mkdir()
    (tmp_path / "notes.txt").write_text("not a folder")

    assert [path.name for path in subdirectories(tmp_path)] == ["Alpha", "zeta"]


def test_folder_picker_uses_nearest_existing_parent(tmp_path: Path):
    missing = tmp_path / "not-created" / "nested"

    assert FolderPickerDialog._usable_initial_directory(missing) == tmp_path


def test_folder_picker_uses_file_parent(tmp_path: Path):
    results = tmp_path / "redock_results.json"
    results.write_text('{"results": []}')

    assert FolderPickerDialog._usable_initial_directory(results) == tmp_path


class _PathVar:
    def __init__(self, value: Path) -> None:
        self.value = value

    def get(self) -> str:
        return str(self.value)


def test_results_folder_button_uses_navigable_picker(monkeypatch, tmp_path: Path):
    selected = tmp_path / "screening_run"
    selected.mkdir()
    app = object.__new__(DockMateVSApp)
    app.output_var = _PathVar(tmp_path)
    loaded = []
    app._load_results_selection = loaded.append
    picker_calls = []

    def fake_picker(parent, **options):
        picker_calls.append((parent, options))
        return selected

    monkeypatch.setattr(app_module, "choose_directory", fake_picker)

    app._browse_results_folder()

    assert loaded == [selected]
    assert picker_calls == [
        (
            app,
            {
                "title": "Select one completed docking run folder",
                "initial_dir": tmp_path,
            },
        )
    ]


def test_pose_folder_button_uses_navigable_picker(monkeypatch, tmp_path: Path):
    selected = tmp_path / "screening_run"
    selected.mkdir()
    app = object.__new__(DockMateVSApp)
    app.output_var = _PathVar(tmp_path)
    loaded = []
    app._load_pose_results_selection = loaded.append

    monkeypatch.setattr(
        app_module,
        "choose_directory",
        lambda parent, **options: selected,
    )

    app._browse_pose_results_folder()

    assert loaded == [selected]
