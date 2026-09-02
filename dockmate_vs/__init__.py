"""DockMate-VS: docking protocol development and virtual screening."""

from importlib.metadata import PackageNotFoundError, version


try:
    __version__ = version("dockmate-vs")
except PackageNotFoundError:
    __version__ = "0.1.0"


__all__ = ["__version__"]
