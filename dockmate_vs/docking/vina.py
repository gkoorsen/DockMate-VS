"""Vina uses the shared PDBQT adapter with its fixed baseline scoring mode."""

from .smina import SminaDockingEngine


class VinaDockingEngine(SminaDockingEngine):
    binary_name = "vina"
    fixed_scoring = "vina"
