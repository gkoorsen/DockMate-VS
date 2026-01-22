#!/usr/bin/env python3
"""
Launch script for redock analysis GUI.

Usage:
    python scripts/launch_redock_analysis_gui.py
"""

import argparse
import sys
from pathlib import Path
import multiprocessing as mp  # ← ADDED: For parallel processing fix

sys.path.insert(0, str(Path(__file__).parent.parent))


def check_dependencies():
    missing = []

    try:
        import tkinter  # noqa: F401
    except ImportError:
        missing.append("tkinter")

    try:
        import pandas  # noqa: F401
    except ImportError:
        missing.append("pandas")

    try:
        from rdkit import Chem  # noqa: F401
    except ImportError:
        missing.append("rdkit")

    try:
        from loguru import logger  # noqa: F401
    except ImportError:
        missing.append("loguru")

    try:
        import openpyxl  # noqa: F401
    except ImportError:
        missing.append("openpyxl")

    if missing:
        print("ERROR: Missing required dependencies:")
        for dep in missing:
            print(f"  - {dep}")
        print("\nInstall missing dependencies:")
        print("  pip install pandas openpyxl rdkit loguru")
        sys.exit(1)


def main():
    check_dependencies()

    parser = argparse.ArgumentParser(description="Launch redock analysis GUI")
    parser.add_argument(
        "--vina",
        action="store_true",
        help="Use Vina (instead of Smina) for adaptive docking"
    )
    parser.add_argument(
        "--vina-binary",
        type=str,
        default=None,
        help="Path to Vina binary"
    )
    parser.add_argument(
        "--smina-binary",
        type=str,
        default=None,
        help="Path to Smina binary"
    )
    parser.add_argument(
        "--serial",
        action="store_true",
        help="Force serial processing (disable parallel RMSD calculation)"
    )
    args = parser.parse_args()

    from loguru import logger
    from docking_platform_gui.gui.redock_analysis import RedockAnalysisApp
    import os

    # Set serial mode if requested
    if args.serial:
        os.environ['LIGAND_PREP_FORCE_SERIAL'] = '1'
        logger.info("Serial processing mode enabled via --serial flag")

    logger.add(
        "redock_gui.log",
        rotation="10 MB",
        level="INFO",
        format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}"
    )

    logger.info("Starting redock analysis GUI")

    try:
        app = RedockAnalysisApp(
            use_vina_default=args.vina,
            vina_binary_default=args.vina_binary,
            smina_binary_default=args.smina_binary
        )
        app.mainloop()
    except Exception as exc:
        logger.error(f"Fatal error: {exc}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    logger.info("Redock analysis GUI closed")


if __name__ == "__main__":
    # ← ADDED: Fix multiprocessing for GUI compatibility
    # Use 'spawn' method instead of 'fork' to avoid deadlocks with Tkinter
    # This enables parallel RMSD processing in ligand preparation
    try:
        mp.set_start_method('spawn', force=True)
    except RuntimeError:
        # Start method already set (e.g., in tests or when imported)
        pass
    
    main()