"""Focused tests for reconstruction input validation; run with unittest."""

import io
import json
import tarfile
import tempfile
import unittest
from pathlib import Path

from reconstruct import MEMBERS, read_sources, reconstruct


class SourceValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.archive = self.root / "source.tar.gz"
        self.extracted = self.root / "extracted"
        self.extracted.mkdir()

    def make_archive(self, special=None, duplicate=False, missing=False, extra=False):
        with tarfile.open(self.archive, "w:gz") as handle:
            for index, name in enumerate(MEMBERS.values()):
                if missing and index == 1:
                    continue
                member = tarfile.TarInfo(name)
                data = b"example source record\n"
                member.size = len(data)
                if special and index == 0:
                    member.type = special
                    member.linkname = "../../outside"
                    member.size = 0
                handle.addfile(member, io.BytesIO(data) if member.isfile() else None)
                if duplicate and index == 0:
                    handle.addfile(member, io.BytesIO(data))
            if extra:
                member = tarfile.TarInfo("../../outside")
                member.size = 1
                handle.addfile(member, io.BytesIO(b"x"))

    def test_reads_only_exact_members(self):
        self.make_archive(extra=True)
        read_sources(self.archive, self.extracted)
        self.assertEqual({p.name for p in self.extracted.iterdir()}, set(MEMBERS))
        self.assertEqual((self.extracted / "actives.smi").read_bytes(), b"example source record\n")

    def test_rejects_symlink(self):
        self.make_archive(special=tarfile.SYMTYPE)
        with self.assertRaisesRegex(ValueError, "Invalid source member"):
            read_sources(self.archive, self.extracted)

    def test_rejects_hardlink(self):
        self.make_archive(special=tarfile.LNKTYPE)
        with self.assertRaisesRegex(ValueError, "Invalid source member"):
            read_sources(self.archive, self.extracted)

    def test_rejects_duplicate_member(self):
        self.make_archive(duplicate=True)
        with self.assertRaisesRegex(ValueError, "Expected one archive member"):
            read_sources(self.archive, self.extracted)

    def test_rejects_missing_member(self):
        self.make_archive(missing=True)
        with self.assertRaisesRegex(ValueError, "Expected one archive member"):
            read_sources(self.archive, self.extracted)

    def test_rejects_wrong_archive_before_creating_output(self):
        self.make_archive()
        (self.root / "source.json").write_text(json.dumps({"archive_sha256": "0" * 64}))
        output = self.root / "new-output"
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            reconstruct(self.archive, output, self.root)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
