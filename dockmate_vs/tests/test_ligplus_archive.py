"""Archive installation uses synthetic files, not licensed LigPlot+ assets."""

import importlib.util
import io
import os
from pathlib import Path
import shlex
import stat
import subprocess
import sys
import tarfile
import zipfile

import pytest

from dockmate_vs.tests.test_external_tool_installer import BASH, INSTALLER


spec = importlib.util.spec_from_file_location("install_ligplus_archive", INSTALLER.with_name("install_ligplus_archive.py"))
archive_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(archive_module)


def make_archive(path, wrapped=True, extra=None, missing_binary=False, version=b"test jar"):
    root = "LigPlus/" if wrapped else ""
    binary = str(archive_module.executable_paths(sys.platform)[0])
    files = {root + "LigPlus.jar": version, root + "lib/data/components.cif": b"test dictionary"}
    if not missing_binary:
        files[root + binary] = b"#!/bin/sh\nexit 0\n"
    files.update(extra or {})
    if path.suffix == ".zip":
        with zipfile.ZipFile(path, "w") as archive:
            for name, data in files.items():
                archive.writestr(name, data)
    else:
        modes = {".tar": "w", ".gz": "w:gz", ".bz2": "w:bz2", ".xz": "w:xz"}
        with tarfile.open(path, modes[path.suffix]) as archive:
            for name, data in files.items():
                member = tarfile.TarInfo(name)
                member.size, member.mode = len(data), 0o644
                archive.addfile(member, io.BytesIO(data))
    return path


@pytest.mark.parametrize("extension", [".zip", ".tar", ".tar.gz", ".tar.bz2", ".tar.xz"])
@pytest.mark.parametrize("wrapped", [False, True])
def test_archive_formats_layouts_permissions_and_reuse(tmp_path, extension, wrapped):
    archive = make_archive(tmp_path / ("LigPlot download" + extension), wrapped)
    prefix = tmp_path / "environment with spaces"
    prefix.mkdir()
    planned = archive_module.install_archive(archive, prefix, check=True)
    assert not (prefix / "opt").exists()
    installed = archive_module.install_archive(archive, prefix)
    assert installed == planned
    assert (installed / "LigPlus.jar").read_bytes() == b"test jar"
    binary = installed / archive_module.executable_paths(sys.platform)[0]
    assert os.access(binary, os.X_OK)
    timestamp = binary.stat().st_mtime_ns
    assert archive_module.install_archive(archive, prefix) == installed
    assert binary.stat().st_mtime_ns == timestamp
    assert not list(installed.parent.glob(".staging-*"))


def test_new_archive_keeps_previous_installation(tmp_path):
    first = make_archive(tmp_path / "first.zip")
    old = archive_module.install_archive(first, tmp_path)
    second = make_archive(tmp_path / "second.zip", version=b"new jar")
    new = archive_module.install_archive(second, tmp_path)
    assert old != new
    assert (old / "LigPlus.jar").read_bytes() == b"test jar"
    assert (new / "LigPlus.jar").read_bytes() == b"new jar"


@pytest.mark.parametrize("name", ["../escaped", "/tmp/escaped", "C:/escaped", "..\\escaped"])
@pytest.mark.parametrize("extension", [".tar", ".zip"])
def test_unsafe_paths_rejected_without_writes(tmp_path, name, extension):
    archive = make_archive(tmp_path / ("bad" + extension), extra={name: b"bad"})
    with pytest.raises(ValueError, match="Unsafe archive path"):
        archive_module.install_archive(archive, tmp_path)
    assert not (tmp_path / "opt").exists()


@pytest.mark.parametrize("kind", ["tar_symlink", "tar_hardlink", "tar_fifo", "zip_symlink"])
def test_links_and_special_files_rejected(tmp_path, kind):
    path = make_archive(tmp_path / ("links.zip" if kind.startswith("zip") else "links.tar"))
    if kind.startswith("zip"):
        with zipfile.ZipFile(path, "a") as archive:
            member = zipfile.ZipInfo("link")
            member.create_system = 3
            member.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive.writestr(member, "../outside")
    else:
        with tarfile.open(path, "a") as archive:
            member = tarfile.TarInfo("link")
            member.type = {"tar_symlink": tarfile.SYMTYPE, "tar_hardlink": tarfile.LNKTYPE,
                           "tar_fifo": tarfile.FIFOTYPE}[kind]
            member.linkname = "../outside"
            archive.addfile(member)
    with pytest.raises(ValueError, match="links and special files"):
        archive_module.install_archive(path, tmp_path)
    assert not (tmp_path / "opt").exists()


def test_incomplete_or_ambiguous_distribution_rejected(tmp_path):
    path = make_archive(tmp_path / "incomplete.zip", missing_binary=True)
    with pytest.raises(ValueError, match="executable for this platform"):
        archive_module.install_archive(path, tmp_path)
    path = make_archive(tmp_path / "ambiguous.zip", extra={"Other/LigPlus.jar": b"other"})
    with pytest.raises(ValueError, match="exactly one LigPlus.jar"):
        archive_module.install_archive(path, tmp_path)


def test_failed_extraction_does_not_publish_or_remove_previous_install(tmp_path, monkeypatch):
    first = make_archive(tmp_path / "first.zip")
    old = archive_module.install_archive(first, tmp_path)
    second = make_archive(tmp_path / "second.zip", version=b"new")
    planned = archive_module.install_archive(second, tmp_path, check=True)
    def fail(*args):
        raise OSError("simulated disk error")
    monkeypatch.setattr(archive_module.shutil, "copyfileobj", fail)
    with pytest.raises(OSError, match="simulated disk error"):
        archive_module.install_archive(second, tmp_path)
    assert (old / "LigPlus.jar").read_bytes() == b"test jar"
    assert not planned.exists()
    assert not list(old.parent.glob(".staging-*"))


@pytest.fixture
def installer_env(tmp_path):
    prefix = tmp_path / "dockmate env"
    binaries = prefix / "bin"
    binaries.mkdir(parents=True)
    for name in ("vina", "smina", "rbdock", "rbcavity", "obabel", "reduce", "mk_prepare_receptor.py",
                 "fpocket", "pymol", "java", "gs", "mamba"):
        path = binaries / name
        path.write_text("#!/bin/sh\nexit 0\n")
        path.chmod(0o755)
    python = binaries / "python"
    python.write_text('#!/bin/sh\nif [ "$1" = "-c" ]; then exit 0; fi\n'
                      f'exec {shlex.quote(sys.executable)} "$@"\n')
    python.chmod(0o755)
    env = dict(os.environ, CONDA_PREFIX=str(prefix), CONDA_DEFAULT_ENV="dockmate-vs",
               PATH=f"{binaries}:/usr/bin:/bin")
    return prefix, env


@pytest.mark.skipif(BASH is None, reason="requires bash")
@pytest.mark.parametrize("dry_run", [False, True])
def test_shell_archive_registration_and_dry_run(tmp_path, installer_env, dry_run):
    prefix, env = installer_env
    archive = make_archive(tmp_path / "LigPlot download.zip")
    args = [BASH, str(INSTALLER), "--package-manager", str(prefix / "bin/mamba"),
            "--ligplus-archive", str(archive)]
    if dry_run:
        args.append("--dry-run")
    result = subprocess.run(args, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    if dry_run:
        assert "install validated archive" in result.stdout
        assert not (prefix / "opt").exists()
        assert not (prefix / "etc").exists()
    else:
        target = (prefix / "bin/ligplot").resolve(strict=True)
        assert target.is_relative_to(prefix / "opt/ligplus")
        hooks = prefix / "etc/conda/activate.d"
        assert "LIGPLUS_ROOT=" in (hooks / "dockmate-ligplus-root.sh").read_text()
        assert "HET_GROUP_DICTIONARY=" in (hooks / "dockmate-ligplus-dictionary.sh").read_text()


@pytest.mark.skipif(BASH is None, reason="requires bash")
@pytest.mark.parametrize("case", ["missing", "conflict", "invalid", "argument_missing"])
def test_shell_archive_validation_before_environment_changes(tmp_path, installer_env, case):
    prefix, env = installer_env
    archive = tmp_path / "input.zip"
    if case == "invalid":
        archive.write_bytes(b"not an archive")
    elif case == "conflict":
        make_archive(archive)
    args = [BASH, str(INSTALLER), "--ligplus-archive"]
    if case != "argument_missing":
        args.append(str(archive))
    if case == "conflict":
        args.extend(["--ligplus-root", str(prefix)])
    result = subprocess.run(args, env=env, capture_output=True, text=True)
    assert result.returncode != 0
    assert not (prefix / "etc").exists()
    assert not (prefix / "opt").exists()
