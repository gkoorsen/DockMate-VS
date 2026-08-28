"""Headless campaign execution for DockMate-VS.

The validated docking implementation currently lives on ``DockMateVSApp``.
This adapter deliberately skips Tk initialisation while reusing those workflow
methods, which keeps GUI and command-line campaigns behaviourally identical.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from dataclasses import fields
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

import numpy as np
import pandas as pd
import yaml
from loguru import logger

from dockmate_vs.gui.app import DockMateVSApp, RedockResult


DEFAULT_CONFIG: Dict[str, Any] = {
    "threshold": 2.0,
    "filters": {
        "exclude_additives": False,
        "exclude_cofactors": False,
        "smiles_handling": "automatic",
    },
    "sampling": {
        "enabled": False,
        "size": None,
        "seed": None,
        "include_all_controls": True,
        "strategy": "stratified_by_structure",
    },
    "single": {
        "engine": "smina",
        "box_margin": 4.0,
        "apo_site_mode": "auto",
        "site_definition_mode": "auto",
        "site_residues": "",
        "size_override": None,
        "water_handling": "remove_all",
        "exhaustiveness": 16,
        "num_modes": 20,
        "energy_range": 3.0,
        "cpu": 4,
        "seed": 42,
        "timeout": 1200,
        "scoring": "vina",
        "smina_binary": "smina",
        "vina_binary": "vina",
        "rdock_root": None,
        "rdock_runs": 20,
        "rdock_seed": 42,
        "rdock_radius": None,
        "ligand_variant_mode": "adaptive",
        "variant_select_by": "score",
        "max_tautomers": 8,
        "max_conformers": 10,
        "n_cpus": 4,
    },
    "adaptive": {},
    "rescore": {
        "enable": False,
        "scoring": "vina",
        "smina_binary": "smina",
    },
    "protocol_sweep": {
        "water_modes": ["remove_all", "retain_all", "selective"],
        "engines": ["smina"],
        "box_definitions": [
            {"label": "margin:4", "box_margin": 4.0, "size_override": None}
        ],
        "rescore_methods": ["none"],
        "exhaustiveness": [8, 16, 32],
        "seeds": [42],
    },
}


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    merged: Dict[str, Any] = {}
    for key, value in base.items():
        merged[key] = _deep_merge(value, {}) if isinstance(value, dict) else value
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _resolve_path(value: Any, base_dir: Path) -> Path:
    path = Path(str(value)).expanduser()
    return path.resolve() if path.is_absolute() else (base_dir / path).resolve()


def _normalise_box_definitions(values: Any) -> list[dict]:
    if isinstance(values, str):
        return DockMateVSApp._parse_protocol_box_definitions(values)
    if not isinstance(values, list) or not values:
        raise ValueError("protocol_sweep.box_definitions must be a non-empty list")

    definitions = []
    for value in values:
        if isinstance(value, str):
            definitions.extend(DockMateVSApp._parse_protocol_box_definitions(value))
            continue
        if not isinstance(value, dict):
            raise ValueError("Each box definition must be a mapping or string")
        box = dict(value)
        size = box.get("size_override")
        if size is not None:
            if not isinstance(size, (list, tuple)) or len(size) != 3:
                raise ValueError("A fixed size_override must contain three dimensions")
            size = [float(item) for item in size]
            if any(item <= 0 for item in size):
                raise ValueError("Box dimensions must be greater than zero")
            box["size_override"] = size
            box.setdefault("box_margin", None)
            box.setdefault("label", "x".join(f"{item:g}" for item in size))
        else:
            margin = float(box.get("box_margin", 4.0))
            if margin <= 0:
                raise ValueError("Box margins must be greater than zero")
            box["box_margin"] = margin
            box["size_override"] = None
            box.setdefault("label", f"margin:{margin:g}")
        definitions.append(box)
    return definitions


def load_campaign_config(path: Path, mode: str) -> dict:
    """Load and validate a YAML/JSON campaign file."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.is_file():
        raise ValueError(f"Campaign config not found: {config_path}")
    try:
        payload = yaml.safe_load(config_path.read_text())
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid campaign YAML: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("Campaign config must contain a top-level mapping")

    config = _deep_merge(DEFAULT_CONFIG, payload)
    config["mode"] = "protocol_development" if mode == "protocol" else "screening"
    if not config.get("input_file"):
        raise ValueError("Campaign config requires input_file")
    if not config.get("output_dir"):
        raise ValueError("Campaign config requires output_dir")

    base_dir = config_path.parent
    config["input_file"] = str(_resolve_path(config["input_file"], base_dir))
    config["output_dir"] = _resolve_path(config["output_dir"], base_dir)
    if not Path(config["input_file"]).is_file():
        raise ValueError(f"Input workbook not found: {config['input_file']}")

    config["threshold"] = float(config["threshold"])
    if config["threshold"] <= 0:
        raise ValueError("threshold must be greater than zero")

    single = config["single"]
    if single.get("engine") not in {"vina", "smina", "rdock"}:
        raise ValueError("single.engine must be vina, smina, or rdock")
    if single.get("water_handling") not in {"remove_all", "retain_all", "selective"}:
        raise ValueError(
            "single.water_handling must be remove_all, retain_all, or selective"
        )
    single["box_margin"] = float(single["box_margin"])
    single["exhaustiveness"] = int(single["exhaustiveness"])
    single["num_modes"] = int(single["num_modes"])
    single["energy_range"] = float(single["energy_range"])
    single["cpu"] = int(single["cpu"])
    single["seed"] = int(single["seed"])
    single["timeout"] = int(single["timeout"])
    single["n_cpus"] = int(single.get("n_cpus") or single["cpu"])
    single["max_tautomers"] = int(single["max_tautomers"])
    single["max_conformers"] = int(single["max_conformers"])
    if single["box_margin"] <= 0:
        raise ValueError("single.box_margin must be greater than zero")
    if min(single["exhaustiveness"], single["num_modes"], single["cpu"], single["n_cpus"]) <= 0:
        raise ValueError(
            "single exhaustiveness, num_modes, cpu, and n_cpus must be greater than zero"
        )
    if not 1 <= single["max_tautomers"] <= 20:
        raise ValueError("single.max_tautomers must be between 1 and 20")
    if not 1 <= single["max_conformers"] <= 30:
        raise ValueError("single.max_conformers must be between 1 and 30")
    if single.get("size_override") is not None:
        size = single["size_override"]
        if not isinstance(size, (list, tuple)) or len(size) != 3:
            raise ValueError("single.size_override must contain three dimensions")
        single["size_override"] = np.asarray(size, dtype=float)
    if not single.get("rdock_root"):
        single["rdock_root"] = os.environ.get("RBT_ROOT") or sys.prefix
    single["rdock_root"] = str(Path(single["rdock_root"]).expanduser())

    sampling = config["sampling"]
    sampling["enabled"] = bool(sampling.get("enabled"))
    if sampling["enabled"]:
        if sampling.get("size") is None:
            raise ValueError("sampling.size is required when sampling is enabled")
        sampling["size"] = int(sampling["size"])
        if sampling["size"] <= 0:
            raise ValueError("sampling.size must be greater than zero")
        if sampling.get("seed") is not None:
            sampling["seed"] = int(sampling["seed"])

    sweep = config["protocol_sweep"]
    sweep["box_definitions"] = _normalise_box_definitions(sweep["box_definitions"])
    for key in ("water_modes", "engines", "rescore_methods", "exhaustiveness", "seeds"):
        if not isinstance(sweep.get(key), list) or not sweep[key]:
            raise ValueError(f"protocol_sweep.{key} must be a non-empty list")
    sweep["exhaustiveness"] = [int(item) for item in sweep["exhaustiveness"]]
    sweep["seeds"] = [int(item) for item in sweep["seeds"]]
    if any(item <= 0 for item in sweep["exhaustiveness"]):
        raise ValueError("protocol_sweep.exhaustiveness values must be greater than zero")
    unsupported_engines = set(sweep["engines"]) - {"vina", "smina", "rdock"}
    if unsupported_engines:
        raise ValueError(
            "Unsupported protocol engine(s): " + ", ".join(sorted(unsupported_engines))
        )
    unsupported_waters = set(sweep["water_modes"]) - {
        "remove_all", "retain_all", "selective"
    }
    if unsupported_waters:
        raise ValueError(
            "Unsupported water mode(s): " + ", ".join(sorted(unsupported_waters))
        )
    DockMateVSApp._parse_protocol_rescore_methods(
        ",".join(str(item) for item in sweep["rescore_methods"])
    )

    # Screening must never select a ligand variant using native-pose RMSD.
    if mode == "screen":
        single["variant_select_by"] = "score"
    return config


class _EventSink:
    def __init__(self) -> None:
        self.terminal_event: Optional[tuple] = None
        self.error: Optional[str] = None
        self._last_progress = (-1, -1)

    def put(self, event: tuple) -> None:
        kind = event[0]
        if kind == "log":
            print(event[1], flush=True)
        elif kind == "progress":
            current, total, label = event[1:4]
            marker = (current, total)
            if marker != self._last_progress:
                print(f"[{current}/{total}] {label}", flush=True)
                self._last_progress = marker
        elif kind in {"preflight_failed", "protocol_incompatible"}:
            self.error = str(event[1])
            if kind == "protocol_incompatible":
                self.error = f"{event[2]} ({event[1]})"
            self.terminal_event = event
        elif kind in {"done", "protocol_done", "cancelled"}:
            self.terminal_event = event


class HeadlessDockMateRunner(DockMateVSApp):
    """Reuse DockMateVSApp workflow methods without creating a Tk root."""

    def __init__(self) -> None:
        # Deliberately do not call DockMateVSApp.__init__ (which creates Tk).
        self.progress_dialog = None
        self._queue = _EventSink()
        self._network_phase_complete = False


def _load_pairs(runner: HeadlessDockMateRunner, config: dict) -> list[dict]:
    filters = config.get("filters", {})
    pairs, columns = runner._load_pairs_from_excel(
        Path(config["input_file"]),
        exclude_additives=bool(filters.get("exclude_additives")),
        exclude_cofactors=bool(filters.get("exclude_cofactors")),
    )
    if not pairs:
        raise ValueError("No valid PDB/ligand pairs were found in the workbook")
    config["columns"] = columns
    return pairs


def run_campaign(config_path: Path, mode: str) -> Path:
    """Execute a protocol-development or screening campaign synchronously."""
    if mode not in {"protocol", "screen"}:
        raise ValueError(f"Unsupported campaign mode: {mode}")
    config = load_campaign_config(config_path, mode)
    runner = HeadlessDockMateRunner()
    pairs = _load_pairs(runner, config)

    if mode == "protocol":
        actives = runner._protocol_active_pairs(pairs)
        if not actives:
            raise ValueError(
                "Protocol development requires at least one active/decoy control pair"
            )
        runner._run_protocol_worker(actives, config)
        expected = config["output_dir"] / "protocol_development" / "protocol_development_results.csv"
    else:
        sampling = config["sampling"]
        if sampling.get("enabled"):
            pairs = runner._apply_random_sample(
                pairs, int(sampling["size"]), sampling.get("seed")
            )
        config["planned_cases"] = {
            "total": len(pairs),
            "actives": sum(pair.get("control_label") == 1 for pair in pairs),
            "decoys": sum(pair.get("control_label") == 0 for pair in pairs),
            "samples": sum(pair.get("control_label") is None for pair in pairs),
        }
        runner._run_worker(pairs, config)
        expected = config["output_dir"] / "redock_results.json"

    if runner._queue.error:
        raise RuntimeError(runner._queue.error)
    if not expected.is_file():
        raise RuntimeError(f"Campaign ended without producing {expected}")
    if mode == "protocol":
        frame = DockMateVSApp._read_results_csv(expected)
        if frame.empty or "status" not in frame or not (frame["status"] == "complete").any():
            raise RuntimeError(
                "Protocol campaign produced no successfully completed conditions; "
                f"review {expected}"
            )
    else:
        payload = json.loads(expected.read_text())
        if not any(item.get("docking_completed") is True for item in payload.get("results", [])):
            raise RuntimeError(
                "Screening campaign produced no successfully completed dockings; "
                f"review {expected}"
            )
    return expected


def _none_if_missing(value: Any) -> Any:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    return value


def _results_from_rows(rows: Iterable[dict]) -> list[RedockResult]:
    allowed = {field.name for field in fields(RedockResult)}
    results = []
    for row in rows:
        values = {key: _none_if_missing(value) for key, value in row.items() if key in allowed}
        if values.get("control_label") is not None:
            values["control_label"] = int(values["control_label"])
        for key in ("success", "docking_completed"):
            if isinstance(values.get(key), str):
                values[key] = values[key].strip().lower() in {"1", "true", "yes"}
        results.append(RedockResult(**values))
    return results


def regenerate_report(run_path: Path, threshold: Optional[float] = None) -> Path:
    """Regenerate a Markdown/JSON report from an existing run."""
    selected = Path(run_path).expanduser().resolve()
    result_path = DockMateVSApp._result_file_for_selection(selected)
    if result_path is None:
        raise ValueError(f"No supported DockMate-VS results found under {selected}")

    runner = HeadlessDockMateRunner()
    if result_path.name == "protocol_development_results.csv":
        rows = DockMateVSApp._read_results_csv(result_path).to_dict("records")
        return runner._write_protocol_report(rows, result_path.parent, threshold or 2.0)

    if result_path.suffix.lower() == ".json":
        payload = json.loads(result_path.read_text())
        rows = payload.get("results", [])
    else:
        rows = DockMateVSApp._read_results_csv(result_path).to_dict("records")
    results = _results_from_rows(rows)

    summary_path = result_path.with_name("redock_summary.json")
    if threshold is None and summary_path.is_file():
        try:
            threshold = float(json.loads(summary_path.read_text()).get("threshold", 2.0))
        except (ValueError, TypeError, json.JSONDecodeError):
            threshold = 2.0
    threshold = threshold or 2.0
    summary = runner._build_summary(results, threshold)
    canonical_results_path = result_path.with_name("redock_results.json")
    runner._write_summary_files(canonical_results_path, summary)
    return canonical_results_path.with_name("redock_summary.md")


def doctor(strict: bool = False) -> int:
    """Print dependency diagnostics; return a process-style status code."""
    required_imports = {
        "numpy": "numpy",
        "pandas": "pandas",
        "openpyxl": "openpyxl",
        "rdkit": "rdkit",
        "Bio": "biopython",
        "yaml": "pyyaml",
    }
    missing_required = []
    print("Python dependencies:")
    for module, package in required_imports.items():
        try:
            __import__(module)
            print(f"  OK      {package}")
        except ImportError:
            print(f"  MISSING {package}")
            missing_required.append(package)

    print("Docking tools:")
    required_tools = ("vina", "smina", "obabel")
    optional_tools = ("rbdock", "rbcavity", "fpocket")
    missing_tools = []
    for name in required_tools + optional_tools:
        location = shutil.which(name)
        status = "OK" if location else ("OPTIONAL" if name in optional_tools else "MISSING")
        print(f"  {status:<9}{name}: {location or 'not found'}")
        if not location and name in required_tools:
            missing_tools.append(name)

    print("Host visualization integrations (not included in the core container):")
    for name in ("pymol", "ligplot"):
        print(f"  {'OK' if shutil.which(name) else 'OPTIONAL':<9}{name}")

    if missing_required or (strict and missing_tools):
        logger.error(
            "Dependency check failed: imports={}, tools={}", missing_required, missing_tools
        )
        return 1
    return 0
