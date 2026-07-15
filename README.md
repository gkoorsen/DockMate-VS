# Redock Analysis GUI

Standalone GUI for redocking analysis with single or adaptive docking, ligand filtering, controls/decoys, and pose inspection.

## Install

> **On Windows?** Follow the step-by-step student guide:
> [docs/SETUP_WSL.md](docs/SETUP_WSL.md) — installs WSL/Ubuntu and the GUI from
> scratch.

### Recommended: conda (Linux/Ubuntu, macOS)

`environment.yml` pins the Python stack **and** the compiled tools (docking
engines, Open Babel, fpocket, OpenMM, ...) that pip cannot install:

```
conda env create -f environment.yml
conda activate docking_platform_gui
redock-gui
```

### Alternative: pip

```
pip install -e .                 # core dependencies
pip install -e '.[prep,pdbqt]'   # + protein-prep pipeline & PDBQT helpers
```

On Debian/Ubuntu the GUI also needs the Tk toolkit, which is a separate system
package:

```
sudo apt install python3-tk
```

The external binaries below still have to be installed separately when using pip
(see "External tools").

You can also launch without installing:

```
python scripts/launch_redock_analysis_gui.py
```

## External tools

The GUI shells out to compiled tools that pip does not install automatically.
`environment.yml` provides most of them; otherwise install with your package
manager before running:

- **Smina or Vina** (required for docking): point the GUI to the binary location.
- **rDock** (optional): set the rDock root folder in the GUI if you want the rDock protocols. Not packaged in `environment.yml`.
- **Open Babel** (`obabel` CLI) is needed for ligand conversions and rescoring.
- **OpenMM**, **PDBFixer**, **Reduce**, and **PROPKA** are used by the protein preparation pipeline. Install via the `prep` extra (`pip install -e '.[prep]'`) or conda:

  ```
  conda install -c conda-forge openmm pdbfixer reduce propka
  ```

  Without them the GUI still launches and docks against already-prepared receptors; the preparation pipeline raises a clear error when invoked.
- **LigPlot / LigPlus** (optional): interaction diagrams. Put the `ligplot` binary on `PATH` or set the `LIGPLOT_BIN` environment variable (the `HET_GROUP_DICTIONARY` variable points to its `components.cif`).

## Python dependencies

Core dependencies are declared in `pyproject.toml`; optional extras (`prep`,
`pdbqt`) cover the heavier preparation stack. You may prefer conda for RDKit,
OpenMM, and the compiled tools.

## Project layout

- `docking_platform_gui/gui/redock_analysis.py`: main GUI.
- `docking_platform_gui/adaptive_docking.py`: adaptive protocol logic.
- `docking_platform_gui/preparation/`: protein/ligand preparation.
- `docking_platform_gui/docking/`: docking engine wrappers.
- `scripts/launch_redock_analysis_gui.py`: launch script.

## Notes

- Rescoring uses `smina --score_only`.
- Adaptive docking uses the 5-protocol flow implemented in `adaptive_docking.py`.
