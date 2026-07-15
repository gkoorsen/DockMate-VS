"""CLI entry point for the redock analysis GUI."""

import argparse
import sys
from pathlib import Path


def check_dependencies() -> None:
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
        print("  pip install -e .            # core Python dependencies")
        if "tkinter" in missing:
            print("  sudo apt install python3-tk # Tk toolkit (Debian/Ubuntu)")
        print("  # or create the full environment: conda env create -f environment.yml")
        sys.exit(1)


def main() -> None:
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
    args = parser.parse_args()

    from loguru import logger
    from docking_platform_gui.gui.redock_analysis import RedockAnalysisApp

    logger.add(
        str(Path("redock_gui.log")),
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
    main()
