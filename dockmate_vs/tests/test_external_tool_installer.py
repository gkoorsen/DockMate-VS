import os
import shutil
import subprocess
from pathlib import Path

import pytest


BASH = shutil.which("bash")
INSTALLER = Path(__file__).resolve().parents[2] / "scripts/install_external_tools.sh"


@pytest.mark.skipif(BASH is None, reason="external-tool installer requires bash")
def test_external_tool_installer_has_valid_syntax_and_help():
    subprocess.run([BASH, "-n", str(INSTALLER)], check=True)
    result = subprocess.run(
        [BASH, str(INSTALLER), "--help"],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "--ligplus-root PATH" in result.stdout
    assert "does not edit shell startup files" in result.stdout


@pytest.mark.skipif(BASH is None, reason="external-tool installer requires bash")
def test_external_tool_installer_dry_run_does_not_modify_environment(tmp_path):
    prefix = tmp_path / "dockmate-vs"
    bin_dir = prefix / "bin"
    bin_dir.mkdir(parents=True)

    python = bin_dir / "python"
    python.write_text("#!/usr/bin/env bash\nexit 1\n")
    python.chmod(0o755)
    conda = tmp_path / "conda"
    conda.write_text("#!/usr/bin/env bash\nexit 0\n")
    conda.chmod(0o755)

    environment = os.environ.copy()
    environment.update(
        {
            "CONDA_PREFIX": str(prefix),
            "CONDA_DEFAULT_ENV": "dockmate-vs",
            "CONDA_EXE": str(conda),
            "PATH": f"{bin_dir}:/usr/bin:/bin",
        }
    )
    result = subprocess.run(
        [BASH, str(INSTALLER), "--dry-run"],
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert "DRY RUN:" in result.stdout
    assert "rdkit\\<2026" in result.stdout
    assert "no files or environments were changed" in result.stdout
    assert not (prefix / "etc").exists()
