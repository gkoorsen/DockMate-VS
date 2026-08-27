"""Command-line entry point for DockMate-VS."""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path
from typing import Optional


def check_gui_dependencies() -> None:
    """Fail early with an actionable message when GUI dependencies are absent."""
    unavailable = []
    checks = (
        ("tkinter", "tkinter"),
        ("pandas", "pandas"),
        ("rdkit", "rdkit"),
        ("loguru", "loguru"),
        ("openpyxl", "openpyxl"),
    )
    for module, package in checks:
        try:
            __import__(module)
        except ImportError as exc:
            unavailable.append((package, str(exc)))

    if unavailable:
        print("ERROR: Required GUI dependencies could not be imported:", file=sys.stderr)
        for dependency, error in unavailable:
            print(f"  - {dependency}: {error}", file=sys.stderr)
        print("\nInstall the DockMate-VS environment from environment.yml.", file=sys.stderr)
        raise SystemExit(1)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dockmate-vs",
        description="Docking protocol development and virtual screening",
    )
    parser.add_argument(
        "--vina", action="store_true", help="Use Vina instead of Smina for adaptive docking"
    )
    parser.add_argument("--vina-binary", default=None, help="Path to the Vina binary")
    parser.add_argument("--smina-binary", default=None, help="Path to the Smina binary")

    commands = parser.add_subparsers(dest="command")
    for name, help_text in (
        ("protocol", "Run a headless protocol-development campaign"),
        ("screen", "Run a headless virtual-screening campaign"),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument(
            "--config", required=True, type=Path, help="Campaign YAML or JSON file"
        )

    report = commands.add_parser("report", help="Regenerate reports for a completed run")
    report.add_argument("--run", required=True, type=Path, help="Run folder or results file")
    report.add_argument(
        "--threshold", type=float, default=None, help="Override the RMSD success threshold"
    )

    doctor = commands.add_parser("doctor", help="Check Python and docking-tool dependencies")
    doctor.add_argument(
        "--strict",
        action="store_true",
        help="Return an error when a core docking executable is unavailable",
    )
    return parser


def _configure_logging() -> None:
    from loguru import logger

    try:
        logger.add(
            str(Path("dockmate-vs.log")),
            rotation="10 MB",
            level="INFO",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {message}",
        )
    except OSError as exc:
        # Read-only container mounts and managed environments should not stop
        # an otherwise valid campaign; Loguru's stderr sink remains available.
        logger.warning("Could not create dockmate-vs.log: {}", exc)


def _run_headless(args: argparse.Namespace) -> int:
    from dockmate_vs.headless import doctor, regenerate_report, run_campaign

    if args.command in {"protocol", "screen"}:
        output = run_campaign(args.config, args.command)
        print(f"Completed: {output}")
        return 0
    if args.command == "report":
        output = regenerate_report(args.run, args.threshold)
        print(f"Report written: {output}")
        return 0
    if args.command == "doctor":
        return doctor(args.strict)
    raise ValueError(f"Unknown command: {args.command}")


def _launch_gui(args: argparse.Namespace) -> int:
    check_gui_dependencies()
    from loguru import logger
    from dockmate_vs.gui.app import DockMateVSApp

    logger.info("Starting DockMate-VS")
    app = DockMateVSApp(
        use_vina_default=args.vina,
        vina_binary_default=args.vina_binary,
        smina_binary_default=args.smina_binary,
    )
    app.mainloop()
    logger.info("DockMate-VS closed")
    return 0


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging()

    try:
        return _run_headless(args) if args.command else _launch_gui(args)
    except KeyboardInterrupt:
        print("Cancelled.", file=sys.stderr)
        return 130
    except Exception as exc:
        try:
            from loguru import logger

            logger.exception("Fatal error: {}", exc)
        except ImportError:
            traceback.print_exc()
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
