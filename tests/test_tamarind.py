"""Tamarind adapter contract tests; scientific parity is checked in the container."""
import copy
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("tamarind_adapter", ROOT / "integrations/tamarind/predict.py")
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)

SITE = {"probability": 0.91, "coordinates_angstrom": [1.1, 2.2, 3.3],
        "coordinating_residues": ["A_3", "A_5", "A_9"],
        "residue_combination": "EHH", "template_id": "example", "template_rmsd_angstrom": 0.2}


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.app = Path(self.temporary.name)
        (self.app / "inputs").mkdir()
        self.structure = self.app / "inputs/sample.cif"
        self.structure.write_text("data_test\n_atom_site.id 1\n")

    def result(self, sites=None, total=None):
        sites = [copy.deepcopy(SITE)] if sites is None else sites
        total = len(sites) if total is None else total
        return {"sites": sites, "site_count": total, "truncated": total > len(sites)}

    def test_uploaded_file_and_scalar_parameters_reach_engine(self):
        predictor = Mock(return_value=self.result())
        result = adapter.run(self.app, {"structure": str(self.structure), "histidine_rotamers": "false",
                                       "minimum_probability": "0.4", "max_sites": "12.0"}, predictor)
        predictor.assert_called_once_with(structure_text=self.structure.read_text(), file_format="cif",
                                          database_backend="pglite", histidine_rotamers=False,
                                          minimum_probability=0.4, max_sites=12)
        saved = json.loads((self.app / "out/prediction.json").read_text())
        self.assertEqual(result, saved)
        self.assertEqual(saved["software_revision"], adapter.ENGINE_REVISION)
        self.assertEqual(len(saved["input_sha256"]), 64)
        with (self.app / "out/sites.csv").open() as stream:
            row = next(csv.DictReader(stream))
        self.assertEqual(json.loads(row["coordinating_residues"]), SITE["coordinating_residues"])
        self.assertEqual(float(row["z_angstrom"]), 3.3)

    def test_staged_file_fallback_and_relative_paths(self):
        for value in ("", "sample.cif", "inputs/sample.cif", str(self.structure)):
            with self.subTest(value=value):
                self.assertEqual(adapter.find_structure(self.app, {"structure": value}), self.structure)
        self.assertEqual(adapter.read_parameters({}),
                         {"histidine_rotamers": True, "minimum_probability": 0.0001, "max_sites": 100})

    def test_ambiguous_missing_or_external_file_is_rejected(self):
        second = self.app / "inputs/second.pdb"
        second.write_text("ATOM  ")
        with self.assertRaises(ValueError):
            adapter.find_structure(self.app, {})
        outside = self.app / "outside.pdb"
        outside.write_text("ATOM  ")
        for path in (outside, self.app / "inputs/missing.pdb"):
            with self.assertRaises(ValueError):
                adapter.find_structure(self.app, {"structure": str(path)})
        link = self.app / "inputs/link.pdb"
        link.symlink_to(outside)
        with self.assertRaises(ValueError):
            adapter.find_structure(self.app, {"structure": str(link)})

    def test_nonfinite_out_of_bounds_and_fractional_settings_are_rejected(self):
        for key, values in {"histidine_rotamers": ["maybe", ""],
                            "minimum_probability": ["nan", "inf", "-0.1", "1.1"],
                            "max_sites": ["nan", "inf", "0", "101", "2.5"]}.items():
            for value in values:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    adapter.read_parameters({key: value})

    def test_empty_and_binary_input_fail_before_prediction(self):
        for data in (b"", b"\xff", b"x" * (adapter.MAX_BYTES + 1)):
            self.structure.write_bytes(data)
            predictor = Mock()
            with self.assertRaises(ValueError):
                adapter.run(self.app, {}, predictor)
            predictor.assert_not_called()

    def test_no_sites_has_valid_json_and_header_only_csv(self):
        result = adapter.run(self.app, {}, Mock(return_value=self.result([])))
        self.assertEqual(result["returned_site_count"], 0)
        self.assertFalse(result["truncated"])
        with (self.app / "out/sites.csv").open() as stream:
            self.assertEqual(list(csv.DictReader(stream)), [])

    def test_truncation_preserves_total_and_returned_counts(self):
        result = adapter.run(self.app, {"max_sites": "1"}, Mock(return_value=self.result(total=10)))
        self.assertEqual((result["site_count"], result["returned_site_count"], result["truncated"]), (10, 1, True))

    def test_engine_failure_does_not_produce_success_artifacts(self):
        with self.assertRaises(RuntimeError):
            adapter.run(self.app, {}, Mock(side_effect=RuntimeError("engine failed")))
        self.assertFalse((self.app / "out/prediction.json").exists())

    def test_failed_worker_logs_are_bounded_and_cleanup_still_runs(self):
        log = io.StringIO()
        with contextlib.redirect_stderr(log), self.assertRaisesRegex(RuntimeError, 'worker failed'):
            with adapter.diagnostic_workdir() as directory:
                work = Path(directory)
                (work / 'worker.log').write_text('discard this prefix\n' + 'x' * 40000 + '\nworker detail')
                (work / 'database.log').write_text('database detail')
                raise RuntimeError('worker failed')
        self.assertFalse(work.exists())
        self.assertIn('worker detail', log.getvalue())
        self.assertIn('database detail', log.getvalue())
        self.assertNotIn('discard this prefix', log.getvalue())
        self.assertLess(len(log.getvalue()), 33000)

    def test_successful_worker_does_not_emit_diagnostics(self):
        log = io.StringIO()
        with contextlib.redirect_stderr(log), adapter.diagnostic_workdir() as directory:
            work = Path(directory)
            (work / 'worker.log').write_text('successful worker')
        self.assertFalse(work.exists())
        self.assertEqual(log.getvalue(), '')


if __name__ == "__main__":
    unittest.main()

