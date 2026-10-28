"""Regression checks for the verifier's success/failure exit status.

Synthetic regeneration replaces the expensive bootstrap subprocess in these
tests. A separate real regeneration is used for release/package verification.
"""
from contextlib import redirect_stdout
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('mw_verifier', Path(__file__).with_name('verify_mw_controlled.py'))
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


class VerifierExitTests(unittest.TestCase):
    def check_exit(self, got_change=None, packaged_change=None):
        baseline = pd.read_csv(ROOT / 'publication/tables/esr1_mw_controlled.csv')
        got, packaged = baseline.copy(), baseline.copy()
        if got_change:
            got_change(got)
        if packaged_change:
            packaged_change(packaged)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'verification').mkdir()
            (root / 'publication/tables').mkdir(parents=True)
            packaged.to_csv(root / 'publication/tables/esr1_mw_controlled.csv', index=False)
            def regenerate(*args, **kwargs):
                got.to_csv(root / 'verification/esr1_mw_controlled.csv', index=False)
            output = io.StringIO()
            with patch.object(verifier, 'ROOT', root), patch.object(verifier.subprocess, 'run', side_effect=regenerate), redirect_stdout(output):
                status = verifier.main()
        return status, output.getvalue()

    def test_matching_table_and_expected_missing_intervals_pass(self):
        self.assertEqual(self.check_exit()[0], 0)

    def test_packaged_only_mismatch_fails(self):
        status, output = self.check_exit(packaged_change=lambda d: d.__setitem__('n', d.n + 1))
        self.assertEqual(status, 1)
        self.assertIn('packaged table differs from regenerated table', output)

    def test_nonfinite_expected_numbers_fail_even_if_tables_agree(self):
        for value in [float('nan'), float('inf'), -float('inf')]:
            with self.subTest(value=value):
                def change(d): d.loc[0, 'auc'] = value
                status, output = self.check_exit(change, change)
                self.assertEqual(status, 1)
                self.assertIn('expected a finite number', output)

    def test_unexpected_values_for_missing_intervals_fail(self):
        for value in [0.0, float('inf')]:
            with self.subTest(value=value):
                def change(d): d.loc[8, 'ci_low'] = value
                status, output = self.check_exit(change, change)
                self.assertEqual(status, 1)
                self.assertIn('expected a missing value', output)

    def test_finite_numerical_mismatch_fails_even_if_tables_agree(self):
        def change(d): d.loc[0, 'auc'] = .9
        self.assertEqual(self.check_exit(change, change)[0], 1)

    def test_extra_rows_fail_even_if_tables_agree(self):
        def change(d):
            d.loc[len(d)] = d.iloc[0]
            d.loc[len(d)-1, 'analysis'] = 'Unexpected analysis'
        self.assertEqual(self.check_exit(change, change)[0], 1)

    def test_duplicate_rows_fail(self):
        def change(d): d.loc[1, 'analysis'] = d.loc[0, 'analysis']
        self.assertEqual(self.check_exit(change, change)[0], 1)


if __name__ == '__main__':
    unittest.main()
