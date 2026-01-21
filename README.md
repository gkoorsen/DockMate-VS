# Redock Analysis GUI

Standalone GUI for redocking analysis with single or adaptive docking, ligand filtering, controls/decoys, and pose inspection.

## Quick start

```
python scripts/launch_redock_analysis_gui.py
```

Or install and run:

```
pip install -e .
redock-gui
```

## External tools

The GUI shells out to compiled tools that pip does not install automatically. Install with your package manager before running:

- **Smina or Vina** (required for docking): point the GUI to the binary location.
- **rDock** (optional): set the rDock root folder in the GUI if you want the rDock protocols.
- **Open Babel** (`obabel` CLI) is needed for ligand conversions and rescoring.
- **OpenMM**, **PDBFixer**, **Reduce**, and **PROPKA** are required/used by the protein preparation pipeline. They are not installed via `pip install -e .`; install them separately, for example with conda:

  ```
  conda install -c conda-forge openmm pdbfixer reduce propka
  ```

  Some installations also require PyMOL/Reduce; refer to each project for macOS/ARM64 binaries.

## Python dependencies

Core dependencies are in `pyproject.toml`. You may prefer conda for RDKit and OpenMM.

## Project layout

- `docking_platform_gui/gui/redock_analysis.py`: main GUI.
- `docking_platform_gui/adaptive_docking.py`: adaptive protocol logic.
- `docking_platform_gui/preparation/`: protein/ligand preparation.
- `docking_platform_gui/docking/`: docking engine wrappers.
- `scripts/launch_redock_analysis_gui.py`: launch script.

## Notes

- Rescoring uses `smina --score_only`.
- Adaptive docking uses the 5-protocol flow implemented in `adaptive_docking.py`.
