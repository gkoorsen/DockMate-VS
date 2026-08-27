# Third-party notices

The optional DockMate-VS core container installs third-party scientific tools
from Conda packages. These tools remain governed by their own licences:

- AutoDock Vina: Apache License 2.0, https://github.com/ccsb-scripps/AutoDock-Vina
- Smina: GNU General Public License v2.0, https://github.com/mwojcikowski/smina
- rDock: GNU Lesser General Public License v3.0, https://github.com/CBDD/rDock
- fpocket: MIT License, https://github.com/Discngine/fpocket
- Open Babel: GNU General Public License v2.0, https://github.com/openbabel/openbabel
- OpenMM: MIT and LGPL components, https://github.com/openmm/openmm
- PDBFixer: MIT License, https://github.com/openmm/pdbfixer

The container is built from the `mambaorg/micromamba` image. Micromamba and the
base-image components retain their own licences and notices:
https://github.com/mamba-org/micromamba-docker

The image does **not** contain PyMOL or LigPlot+. DockMate-VS launches those as
optional host applications. LigPlot+ must be obtained separately under the
licence offered by its authors; it must not be copied into a redistributed
DockMate-VS image without permission.

Python and Conda dependencies listed in `pyproject.toml`, `environment.yml`, and
`docker/environment.yml` are not relicensed by DockMate-VS. Package metadata and
licence files installed in the environment provide the exact notices for the
versions resolved during a build.

The two ACES example workbooks contain a compact, deterministically selected
subset of active and matched-decoy structures from DUD-E. DUD-E describes the
database as free to use; the source data remain attributable to the Irwin and
Shoichet Laboratories and are not relicensed under DockMate-VS's MIT licence.
See https://dude.docking.org/targets/aces and cite Mysinger MM, Carchia M, Irwin
JJ, Shoichet BK, *J. Med. Chem.* 2012, https://doi.org/10.1021/jm300687e. The
workbook Metadata sheets and `scripts/prepare_dude_aces_examples.py` record the
source checksums and deterministic subset procedure.
