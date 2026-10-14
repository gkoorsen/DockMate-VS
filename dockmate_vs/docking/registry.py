"""Explicit engine registration shared by campaigns, the GUI, and diagnostics."""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional, Tuple, Type
import re
import shutil
from collections.abc import Mapping as MappingABC

from .base import DockingEngine


@dataclass(frozen=True)
class EngineCapabilities:
    supports_exhaustiveness: bool = True
    requires_cavity_preparation: bool = False
    site_definition_label: str = "not_applicable"
    supports_score_only_rescoring: bool = True
    binary_names: Tuple[str, ...] = ()
    receptor_format: str = "pdbqt"
    ligand_format: str = "pdbqt"
    output_suffix: str = ".pdbqt"
    scoring_functions: Tuple[str, ...] = ()
    score_units: str = ""
    docking_score_label: str = ""
    color: str = "#60758A"


@dataclass(frozen=True)
class EngineRegistration:
    name: str
    adapter: Type[DockingEngine]
    capabilities: EngineCapabilities
    binary_config_keys: Tuple[str, ...] = ()
    root_config_key: Optional[str] = None
    version_args: Tuple[str, ...] = ("--version",)
    required_by_doctor: bool = False

    def binary_paths(self, config: Optional[Mapping] = None) -> dict:
        config = config or {}
        options = (config.get("engine_options") or {}).get(self.name, {})
        overrides = options.get("binaries", {})
        root = config.get(self.root_config_key) if self.root_config_key else None
        paths = {}
        for index, binary in enumerate(self.capabilities.binary_names):
            key = self.binary_config_keys[index] if index < len(self.binary_config_keys) else None
            candidate = overrides.get(binary) or (config.get(key) if key else None)
            if not candidate and root:
                candidate = str(Path(root).expanduser() / "bin" / binary)
            candidate = str(Path(candidate).expanduser()) if candidate else binary
            paths[binary] = shutil.which(candidate) or candidate
        return paths


_ENGINES = {}


def register_engine(name, adapter, capabilities=None, **metadata):
    """Register an adapter explicitly; accidental replacement is an error."""
    if not isinstance(name, str) or not re.fullmatch(r"[a-z][a-z0-9_-]*", name):
        raise ValueError("Engine names must be lowercase identifiers")
    if name in _ENGINES:
        raise ValueError(f"Engine already registered: {name}")
    if not isinstance(adapter, type) or not issubclass(adapter, DockingEngine):
        raise TypeError("An engine adapter must implement DockingEngine")
    capabilities = capabilities or EngineCapabilities()
    if capabilities.receptor_format not in {"pdb", "pdbqt"}:
        raise ValueError("Supported receptor formats are pdb and pdbqt")
    if capabilities.ligand_format not in {"sdf", "pdbqt"}:
        raise ValueError("Supported ligand formats are sdf and pdbqt")
    if capabilities.output_suffix not in {".pdbqt", ".sdf", ".sd"}:
        raise ValueError("Engine output must be PDBQT or SDF/SD")
    registration = EngineRegistration(name, adapter, capabilities, **metadata)
    _ENGINES[name] = registration
    return registration


def registered_engines():
    return tuple(_ENGINES)


def get_engine(name):
    key = str(getattr(name, "value", name)).lower()
    try:
        return _ENGINES[key]
    except KeyError:
        raise ValueError(f"Unknown docking engine {name!r}; registered: {', '.join(_ENGINES)}") from None


def capabilities_for(name):
    """Unknown historical engine names remain displayable, but cannot execute."""
    registration = _ENGINES.get(str(name).lower())
    return registration.capabilities if registration else EngineCapabilities()


def create_engine(name, config):
    registration = get_engine(name)
    validate_engine_options(config.get("engine_options"))
    settings = dict(config)
    settings["engine_options"] = dict((config.get("engine_options") or {}).get(registration.name, {}))
    settings["engine_options"]["binaries"] = registration.binary_paths(config)
    return registration.adapter.from_config(settings)


def validate_engine_options(options):
    """Validate the common options envelope; adapters validate their own values."""
    if options is None:
        return
    if not isinstance(options, MappingABC):
        raise ValueError("engine_options must be a mapping keyed by registered engine name")
    for name, values in options.items():
        registration = get_engine(name)
        if not isinstance(values, MappingABC):
            raise ValueError(f"engine_options.{name} must be a mapping")
        binaries = values.get("binaries", {})
        if not isinstance(binaries, MappingABC):
            raise ValueError(f"engine_options.{name}.binaries must be a mapping")
        for binary, path in binaries.items():
            if binary not in registration.capabilities.binary_names:
                raise ValueError(f"Unknown binary {binary!r} for engine {name}")
            if not isinstance(path, str) or not path.strip():
                raise ValueError(f"Binary override {name}.{binary} requires a path")


def docking_score_label(name, scoring=None):
    capability = capabilities_for(name)
    label = capability.docking_score_label
    if len(capability.scoring_functions) > 1:
        label = str(scoring or "").strip()
        if not label:
            return f"{str(name).capitalize()} docking (method unknown)"
        label = {"vina": "Vina", "vinardo": "Vinardo"}.get(label.lower(), label)
    return f"{label or str(name).capitalize() or 'Unknown'} (docking)"
