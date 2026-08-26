import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from dockmate_vs.cli import main
from dockmate_vs.headless import (
    HeadlessDockMateRunner,
    load_campaign_config,
    regenerate_report,
    run_campaign,
)
from dockmate_vs.gui.app import DockMateVSApp


def _write_screening_workbook(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "PDB_ID": "1ABC",
                "Ligand": "LIG",
                "Target_Ligand": "active",
                "SMILES": "CCO",
                "decoy compound": "decoy",
                "decoy SMILES": "CCC",
            },
            {
                "PDB_ID": "1ABC",
                "Ligand": "LIG",
                "Target_Ligand": "sample",
                "SMILES": "CCN",
            },
        ]
    ).to_excel(path, index=False)


def test_campaign_config_resolves_paths_and_forces_screening_score_selection(tmp_path):
    workbook = tmp_path / "campaign.xlsx"
    _write_screening_workbook(workbook)
    config_path = tmp_path / "campaign.yml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "input_file": workbook.name,
                "output_dir": "results",
                "single": {
                    "variant_select_by": "rmsd",
                    "size_override": [20, 22, 24],
                },
                "protocol_sweep": {"box_definitions": ["margin:3", "20x22x24"]},
            }
        )
    )

    config = load_campaign_config(config_path, "screen")

    assert config["mode"] == "screening"
    assert config["input_file"] == str(workbook)
    assert config["output_dir"] == tmp_path / "results"
    assert config["single"]["variant_select_by"] == "score"
    assert config["single"]["size_override"].tolist() == [20.0, 22.0, 24.0]
    assert [box["label"] for box in config["protocol_sweep"]["box_definitions"]] == [
        "margin:3",
        "20x22x24",
    ]


def test_headless_runner_does_not_initialise_tk():
    runner = HeadlessDockMateRunner()

    assert runner.progress_dialog is None
    assert runner._network_phase_complete is False
    assert "tk" not in runner.__dict__


def test_screen_campaign_uses_existing_pipeline_without_tk(monkeypatch, tmp_path):
    workbook = tmp_path / "campaign.xlsx"
    _write_screening_workbook(workbook)
    output_dir = tmp_path / "results"
    config_path = tmp_path / "campaign.yml"
    config_path.write_text(
        yaml.safe_dump({"input_file": str(workbook), "output_dir": str(output_dir)})
    )
    observed = {}

    def fake_worker(self, pairs, config):
        observed["pairs"] = pairs
        observed["config"] = config
        output_dir.mkdir(parents=True)
        (output_dir / "redock_results.json").write_text(
            '{"results": [{"docking_completed": true}]}'
        )
        self._queue.put(("done", output_dir / "redock_results.json"))

    monkeypatch.setattr(HeadlessDockMateRunner, "_run_worker", fake_worker)

    result = run_campaign(config_path, "screen")

    assert result == output_dir / "redock_results.json"
    assert len(observed["pairs"]) == 3
    assert observed["config"]["planned_cases"] == {
        "total": 3,
        "actives": 1,
        "decoys": 1,
        "samples": 1,
    }


def test_report_command_regenerates_screening_outputs(tmp_path):
    result_path = tmp_path / "redock_results.json"
    row = {
        "pdb_id": "1ABC",
        "ligand_resname": "LIG",
        "ligand_chain": "A",
        "mode": "screening",
        "engine": "smina",
        "protocol": "single",
        "best_rmsd": 999.9,
        "success": False,
        "runtime_sec": 1.0,
        "best_score": -7.5,
        "dock_name": "sample",
        "docking_completed": True,
    }
    result_path.write_text(json.dumps({"results": [row]}))

    report = regenerate_report(tmp_path)

    assert report == tmp_path / "redock_summary.md"
    assert "Docking Analysis Summary" in report.read_text()
    assert json.loads((tmp_path / "redock_summary.json").read_text())["total_cases"] == 1


def test_cli_returns_nonzero_for_missing_campaign(tmp_path, capsys):
    status = main(["screen", "--config", str(tmp_path / "missing.yml")])

    assert status == 1
    assert "Campaign config not found" in capsys.readouterr().err


def test_config_rejects_sampling_without_size(tmp_path):
    workbook = tmp_path / "campaign.xlsx"
    _write_screening_workbook(workbook)
    config_path = tmp_path / "campaign.yml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "input_file": str(workbook),
                "output_dir": str(tmp_path / "results"),
                "sampling": {"enabled": True},
            }
        )
    )

    with pytest.raises(ValueError, match="sampling.size"):
        load_campaign_config(config_path, "screen")


def test_gui_docker_command_preserves_host_paths_for_pose_loading(tmp_path):
    input_file = tmp_path / "inputs" / "campaign.xlsx"
    input_file.parent.mkdir()
    input_file.write_text("placeholder")
    output_dir = tmp_path / "results"
    campaign_file = output_dir / "dockmate_campaign.screen.yml"

    command = DockMateVSApp._docker_run_command(
        image="dockmate-vs:test",
        platform_name="linux/amd64",
        input_file=input_file,
        output_dir=output_dir,
        campaign_file=campaign_file,
        command_name="screen",
    )

    assert f"{input_file.parent}:{input_file.parent}:ro" in command
    assert f"{output_dir}:{output_dir}" in command
    assert command[command.index("-w") + 1] == str(output_dir)
    assert str(campaign_file) in command
    assert command[-3:] == ["screen", "--config", str(campaign_file)]


def test_gui_docker_command_avoids_duplicate_mount_for_shared_input_output(tmp_path):
    input_file = tmp_path / "campaign.xlsx"
    campaign_file = tmp_path / "dockmate_campaign.screen.yml"

    command = DockMateVSApp._docker_run_command(
        "dockmate-vs:test", "linux/amd64", input_file, tmp_path, campaign_file, "screen"
    )

    mounts = [command[index + 1] for index, item in enumerate(command) if item == "-v"]
    assert mounts == [f"{tmp_path}:{tmp_path}"]
