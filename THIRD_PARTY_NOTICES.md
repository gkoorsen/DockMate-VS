# Third-party notices

The optional DockMate-VS core container installs third-party scientific tools
from Conda packages. These tools remain governed by their own licences:

- AutoDock Vina: Apache License 2.0, https://github.com/ccsb-scripps/AutoDock-Vina
- Smina: GNU General Public License v2.0, https://github.com/mwojcikowski/smina
- rDock: GNU Lesser General Public License v3.0, https://github.com/CBDD/rDock
- fpocket: MIT License, https://github.com/Discngine/fpocket
- Open Babel: GNU General Public License v2.0, https://github.com/openbabel/openbabel
- OpenMM: MIT and LGPL components, https://github.com/openmm/openmm

The image does **not** contain PyMOL or LigPlot+. DockMate-VS launches those as
optional host applications. LigPlot+ must be obtained separately under the
licence offered by its authors; it must not be copied into a redistributed
DockMate-VS image without permission.

Conda package metadata and licence files inside the image provide the exact
notices for the versions installed during a build.
