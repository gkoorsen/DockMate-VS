"""Validate and stage a user-supplied LigPlot+ archive without downloading it."""

import argparse
import hashlib
from pathlib import Path, PurePosixPath
import shutil
import stat
import sys
import tarfile
import tempfile
import zipfile


def member_path(name):
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError(f"Unsafe archive path: {name}")
    return path


def executable_paths(system):
    platforms = ("exe_mac64", "exe_mac") if system == "darwin" else ("exe_linux64", "exe_linux")
    return [PurePosixPath("lib", platform, "ligplot") for platform in platforms]


def archive_entries(archive):
    entries = []
    seen = set()
    members = archive.infolist() if isinstance(archive, zipfile.ZipFile) else archive.getmembers()
    for member in members:
        if isinstance(member, zipfile.ZipInfo):
            name, directory = member.filename, member.is_dir()
            mode = member.external_attr >> 16
            if stat.S_IFMT(mode) not in (0, stat.S_IFREG, stat.S_IFDIR):
                raise ValueError(f"Archive links and special files are not supported: {name}")
        else:
            name, directory, mode = member.name, member.isdir(), member.mode
            if not (member.isfile() or directory):
                raise ValueError(f"Archive links and special files are not supported: {name}")
        path = member_path(name)
        if path == PurePosixPath(".") and directory:
            continue
        if path in seen:
            raise ValueError(f"Duplicate archive path: {name}")
        seen.add(path)
        entries.append((member, path, directory, mode))
    return entries


def install_archive(archive_path, prefix, check=False):
    archive_path = archive_path.resolve(strict=True)
    prefix = prefix.resolve(strict=True)
    if archive_path.name.endswith(".zip"):
        archive = zipfile.ZipFile(archive_path)
    elif archive_path.name.endswith((".tar.gz", ".tgz", ".tar.bz2", ".tbz2", ".tar.xz", ".txz", ".tar")):
        archive = tarfile.open(archive_path, "r:*")
    else:
        raise ValueError("Unsupported LigPlot+ archive; use .zip, .tar, .tar.gz, .tar.bz2, or .tar.xz")

    with archive:
        entries = archive_entries(archive)
        files = {path for _, path, directory, _ in entries if not directory}
        roots = [path.parent for path in files if path.name == "LigPlus.jar"]
        if len(roots) != 1:
            raise ValueError("Expected exactly one LigPlus.jar in the archive")
        root = roots[0]
        binaries = executable_paths(sys.platform)
        if not any(root / binary in files for binary in binaries):
            raise ValueError("LigPlot executable for this platform was not found in the archive")

        digest = hashlib.sha256()
        with archive_path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                digest.update(block)
        destination = prefix / "opt" / "ligplus" / digest.hexdigest()
        if destination.exists() or destination.is_symlink():
            if destination.is_symlink() or not (destination / "LigPlus.jar").is_file():
                raise ValueError(f"Existing archive installation is invalid: {destination}")
            if not any((destination / binary).is_file() and
                       (destination / binary).stat().st_mode & 0o111 for binary in binaries):
                raise ValueError(f"Existing installation has no executable LigPlot binary: {destination}")
            return destination
        if check:
            return destination

        # Stage on the destination filesystem; publish only after extraction succeeds.
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".staging-", dir=destination.parent) as temporary:
            staging = Path(temporary)
            for member, path, directory, mode in entries:
                target = staging / path
                if directory:
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                source = archive.open(member) if isinstance(archive, zipfile.ZipFile) else archive.extractfile(member)
                with source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
                target.chmod(0o644 | (mode & 0o111))
            extracted = staging / root
            # ZIP distributions may not retain executable permission bits.
            for binary in extracted.glob("lib/exe_*/*"):
                if binary.is_file():
                    binary.chmod(binary.stat().st_mode | 0o111)
            extracted.rename(destination)
        return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("prefix", type=Path)
    parser.add_argument("--check", action="store_true", help="Validate and print the destination without writing files")
    args = parser.parse_args()
    try:
        print(install_archive(args.archive, args.prefix, args.check))
    except (OSError, ValueError, tarfile.TarError, zipfile.BadZipFile, RuntimeError, EOFError) as exc:
        parser.exit(1, f"LigPlot+ archive error: {exc}\n")


if __name__ == "__main__":
    main()
