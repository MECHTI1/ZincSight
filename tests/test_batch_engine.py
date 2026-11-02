"""Batch identity, normalization, filtering and subprocess lifecycle contracts."""
from contextlib import contextmanager, redirect_stderr
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, ROOT / "integrations/tamarind" / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


engine = load("batch_engine", "batch_engine.py")
worker = load("batch_worker", "batch_worker.py")
PDB = """MODEL        1
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C
ATOM      2  CA  ALA A   2       2.000   3.000   4.000  1.00 20.00           C
ATOM      3  CA  ALA A   3       3.000   4.000   5.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      4  CA  ALA A   1      91.000  92.000  93.000  1.00 20.00           C
ENDMDL
END
"""
FIELDS = ["structure_id", "prob", "predicted_ion_pos", "chain_resi_1", "chain_resi_2",
          "chain_resi_3", "chain_resi_4", "resi_comb_ordered", "template_id", "template_rmsd"]


def row(key, probability, template):
    return dict(zip(FIELDS, [key, probability, "[1.1, 2.2, 3.3]", "A_1", "A_2", "A_3", "",
                             "{C,H,H}", template, 0.25]))


class BatchEngineTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.work = Path(directory.name)
        self.item = {"key": "a", "source": "one.pdb", "structure_text": PDB, "file_format": "pdb"}

    def validate(self, structures=None, **changes):
        options = dict(histidine_rotamers=True, minimum_probability=0.0001,
                       max_sites=None, cpu_workers=4, timeout=1740)
        options.update(changes)
        return engine._validate_inputs([self.item] if structures is None else structures, **options)

    def fake_engine(self):
        root = self.work / "engine"
        (root / "src").mkdir(parents=True)
        support = root / "integrations/mcp"
        support.mkdir(parents=True)
        for filename in ("worker.py", "embedded_database.py"):
            (support / filename).write_text("# fake worker support\n")
        return root

    def test_invalid_duplicate_keys_content_and_parameters_fail_before_starting_engine(self):
        for key in ("../escape", "x/y", "", "a.b", "x" * 65):
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.validate([{**self.item, "key": key}])
        for structures in ([], [self.item, self.item]):
            with self.assertRaises(ValueError):
                self.validate(structures)
        for option, value in (("histidine_rotamers", "false"), ("max_sites", True),
                              ("max_sites", 0), ("cpu_workers", 9), ("cpu_workers", 1.5),
                              ("minimum_probability", float("nan")), ("timeout", -1)):
            with self.subTest(option=option), self.assertRaises(ValueError):
                self.validate(**{option: value})
        for content in ("", "ATOM  \x00", "not a PDB"):
            with self.assertRaises(ValueError):
                self.validate([{**self.item, "structure_text": content}])

    def test_native_batches_are_not_limited_by_preview_counts_or_byte_ceilings(self):
        structures = [{**self.item, "key": f"s{number}"} for number in range(137)]
        self.assertEqual(len(self.validate(structures)), 137)
        # Reuse one string object to exercise >5 MB each and >100 MB combined
        # without making a second full batch-sized allocation in this test.
        content = "REMARK 123 padding\n" * 300_000 + PDB
        structures = [{**self.item, "key": f"s{number}", "structure_text": content} for number in range(21)]
        self.assertGreater(len(content.encode()), 5_000_000)
        self.assertGreater(len(content.encode()) * len(structures), 100_000_000)
        self.assertEqual(len(self.validate(structures, max_sites=1000)), 21)

    def test_disk_staged_structures_validate_and_copy_without_batch_text_materialization(self):
        original = self.work / 'private-staged.pdb'
        original.write_text(PDB)
        item = {key: value for key, value in self.item.items() if key != 'structure_text'}
        item['structure_path'] = original
        validated, = self.validate([item])
        self.assertEqual(validated['structure_path'], original.resolve())
        self.assertNotIn('structure_text', validated)
        root = self.fake_engine()
        pinned = SimpleNamespace(LIMITATIONS='limitations', parse_sites=Mock())

        def launch(*args, **kwargs):
            work = Path(kwargs['cwd'])
            self.assertEqual((work / 'inputs/a.pdb').read_bytes(), original.read_bytes())
            self.assertEqual(json.loads((work / 'batch-inputs.json').read_text())[0]['path'], 'inputs/a.pdb')
            (work / 'status.json').write_text(json.dumps({'has_sites': False, 'protein_residue_counts': {'a': 3}}))
            return SimpleNamespace(pid=12345, wait=Mock(return_value=0))

        with patch.object(engine, '_load_runner', return_value=pinned), \
                patch.object(engine.subprocess, 'Popen', side_effect=launch), patch.object(engine.os, 'killpg'):
            result, = engine.run_predictions([item], histidine_rotamers=False,
                                             minimum_probability=0, cpu_workers=8, engine_root=root)
        self.assertEqual(result['source'], 'one.pdb')

    def test_disk_staged_missing_ambiguous_link_and_binary_files_are_rejected(self):
        item = {key: value for key, value in self.item.items() if key != 'structure_text'}
        original = self.work / 'input.pdb'
        original.write_text(PDB)
        link = self.work / 'linked.pdb'
        link.symlink_to(original)
        for path in (self.work, self.work / 'missing.pdb', link, ''):
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.validate([{**item, 'structure_path': path}])
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            self.validate([{**self.item, 'structure_path': original}])
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            self.validate([item])
        for contents in (b'', b'ATOM  \x00', b'ATOM  \xff', b'not a structure'):
            original.write_bytes(contents)
            with self.subTest(contents=contents), self.assertRaises(ValueError):
                self.validate([{**item, 'structure_path': original}])

    def test_first_model_count_excludes_nucleic_acids_and_hetero_residues(self):
        original = self.work / "input.pdb"
        original.write_text(PDB)
        query = self.work / "query"
        query.mkdir()
        writer = Mock()
        item = {**self.item, "path": "input.pdb"}
        path, count = worker.normalize_structure(item, self.work, query, writer)
        self.assertEqual((path, count), (query / "a.cif", 3))
        model, input_path, file_format, output_path = writer.call_args.args
        self.assertEqual(len(list(model.get_residues())), 3)
        self.assertEqual((input_path, file_format, output_path), (original, "pdb", path))
        hetero = "HETATM    5  CA  ALA H   4       5.000   6.000   7.000  1.00 20.00           C\n"
        original.write_text(PDB.replace("ENDMDL", hetero + "ENDMDL", 1))
        self.assertEqual(worker.normalize_structure(item, self.work, query, writer)[1], 3)
        original.write_text(PDB.replace("ALA", " DA"))
        with self.assertRaisesRegex(ValueError, "Structure a.*found 0"):
            worker.normalize_structure(item, self.work, query, writer)
        original.write_text("\n".join(PDB.splitlines()[1:3]) + "\nEND\n")
        with self.assertRaisesRegex(ValueError, "at least 3.*found 2"):
            worker.normalize_structure(item, self.work, query, writer)

    def test_first_model_has_no_artificial_maximum_protein_residue_count(self):
        original = self.work / "large.pdb"
        original.write_text("\n".join(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C"
            for i in range(1, 5002)) + "\nEND\n")
        query = self.work / "query"
        query.mkdir()
        writer = Mock()
        item = {**self.item, "path": "large.pdb"}
        self.assertEqual(worker.normalize_structure(item, self.work, query, writer)[1], 5001)
        writer.assert_called_once()

    def test_one_database_and_one_native_call_receive_all_inputs_and_requested_cpu_count(self):
        manifest = [{"key": key, "source": key + ".pdb", "path": key + ".pdb", "file_format": "pdb"}
                    for key in ("a", "b")]
        (self.work / "batch-inputs.json").write_text(json.dumps(manifest))
        entered = Mock()

        @contextmanager
        def embedded_database(root, launcher):
            entered(root, launcher)
            yield {}

        def predict(paths, rotamers, output, cpus, create_pse):
            self.assertEqual(paths, [str(self.work / "query" / (key + ".cif")) for key in ("a", "b")])
            self.assertEqual((rotamers, cpus, create_pse), (True, 4, False))
            target = Path(output) / "table"
            target.mkdir(parents=True)
            (target / "table_predicted_zn_binding_sites.csv").write_text("structure_id,prob\n")
            return "results.zip"

        native = Mock(side_effect=predict)
        modules = {"embedded_database": SimpleNamespace(embedded_database=embedded_database),
                   "pinned_worker": SimpleNamespace(write_first_model=Mock()),
                   "src": ModuleType("src"), "src.primary_script": SimpleNamespace(main=native)}
        with patch.dict(sys.modules, modules), patch.object(worker, "__file__", str(self.work / "batch_worker.py")), \
                patch.object(sys, "argv", ["batch_worker.py", "1", "4"]), \
                patch.dict(os.environ, {"ZINCSIGHT_PGLITE_LAUNCHER": "/pinned/server.mjs"}), \
                patch.object(worker, "normalize_structure", side_effect=lambda item, root, query, *args: (query / (item["key"] + ".cif"), 3)):
            worker.main()
        entered.assert_called_once_with(self.work, Path("/pinned/server.mjs"))
        native.assert_called_once()
        self.assertEqual(json.loads((self.work / "status.json").read_text()),
                         {"has_sites": True, "protein_residue_counts": {"a": 3, "b": 3}})

    def test_structure_grouping_uses_pinned_parser_and_keeps_zero_site_inputs(self):
        roots = [Path("/opt/zincsight"), ROOT.parent / "socket-engine"]
        pinned = next((path for path in roots if (path / "integrations/mcp/runner.py").is_file()), None)
        if pinned is None:
            self.skipTest("Requires the pinned scientific engine installed by the Dockerfile.")
        parse_sites = engine._load_runner(pinned).parse_sites
        table = self.work / "predictions.csv"
        with table.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows([row("b", 0.8, "B"), row("a", 0.4, "low"), row("a", 0.9, "A")])
        grouped = engine._parse_batch_sites(table, ["a", "b", "c"], 0.5, parse_sites)
        self.assertEqual(list(grouped), ["a", "b", "c"])
        self.assertEqual([[site["template_id"] for site in grouped[key]] for key in grouped], [["A"], ["B"], []])
        with table.open("a", newline="") as stream:
            csv.DictWriter(stream, fieldnames=FIELDS).writerow(row("unknown", 0.9, "unexpected"))
        with self.assertRaisesRegex(ValueError, "unknown structure_id"):
            engine._parse_batch_sites(table, ["a", "b", "c"], 0.5, parse_sites)

    def test_no_sites_returns_original_order_and_independent_provenance(self):
        root = self.fake_engine()
        pinned = SimpleNamespace(LIMITATIONS="limitations", parse_sites=Mock())

        def launch(*args, **kwargs):
            work = Path(kwargs["cwd"])
            (work / "status.json").write_text(json.dumps({"has_sites": False,
                                                        "protein_residue_counts": {"b": 31, "a": 3}}))
            return SimpleNamespace(pid=12345, wait=Mock(return_value=0))

        with patch.object(engine, "_load_runner", return_value=pinned), \
                patch.object(engine.subprocess, "Popen", side_effect=launch) as popen, \
                patch.object(engine.os, "killpg") as kill:
            result = engine.run_predictions([self.item, {**self.item, "key": "b", "source": "two.pdb"}],
                                            histidine_rotamers=True, minimum_probability=0.5,
                                            max_sites=1, cpu_workers=4, engine_root=root)
        self.assertEqual([item["source"] for item in result], ["one.pdb", "two.pdb"])
        self.assertEqual([item["protein_residue_count"] for item in result], [3, 31])
        self.assertTrue(all(item["sites"] == [] and item["site_count"] == 0 and not item["truncated"] for item in result))
        pinned.parse_sites.assert_not_called()
        popen.assert_called_once()
        self.assertTrue(popen.call_args.kwargs["start_new_session"])
        self.assertTrue(popen.call_args.kwargs["env"]["ZINCSIGHT_PGLITE_LAUNCHER"].endswith("integrations/mcp/pglite/server.mjs"))
        kill.assert_called_once()

    def test_threshold_counts_and_truncation_are_independent_per_structure(self):
        root = self.fake_engine()
        sites = [{"probability": 0.9}, {"probability": 0.8}]
        pinned = SimpleNamespace(LIMITATIONS="limitations", parse_sites=Mock())

        def launch(*args, **kwargs):
            work = Path(kwargs["cwd"])
            (work / "status.json").write_text(json.dumps({"has_sites": True,
                                                        "protein_residue_counts": {"a": 3, "b": 4}}))
            return SimpleNamespace(pid=12345, wait=Mock(return_value=0))

        with patch.object(engine, "_load_runner", return_value=pinned), \
                patch.object(engine.subprocess, "Popen", side_effect=launch), \
                patch.object(engine.os, "killpg"), \
                patch.object(engine, "_parse_batch_sites", return_value={"a": sites, "b": []}) as parse:
            result = engine.run_predictions([self.item, {**self.item, "key": "b"}],
                                            histidine_rotamers=False, minimum_probability=0.5,
                                            max_sites=1, cpu_workers=4, engine_root=root)
        self.assertEqual([(item["site_count"], len(item["sites"]), item["truncated"]) for item in result],
                         [(2, 1, True), (0, 0, False)])
        self.assertEqual(parse.call_args.args[1:], (["a", "b"], 0.5, pinned.parse_sites))

    def test_default_returns_all_sites_above_one_hundred_without_truncation(self):
        root = self.fake_engine()
        sites = [{"probability": 0.9, "template_id": str(i)} for i in range(211)]
        pinned = SimpleNamespace(LIMITATIONS="limitations", parse_sites=Mock())

        def launch(*args, **kwargs):
            work = Path(kwargs["cwd"])
            (work / "status.json").write_text(json.dumps({"has_sites": True,
                                                        "protein_residue_counts": {"a": 5001}}))
            return SimpleNamespace(pid=12345, wait=Mock(return_value=0))

        with patch.object(engine, "_load_runner", return_value=pinned), \
                patch.object(engine.subprocess, "Popen", side_effect=launch), \
                patch.object(engine.os, "killpg"), \
                patch.object(engine, "_parse_batch_sites", return_value={"a": sites}):
            result, = engine.run_predictions([self.item], histidine_rotamers=False,
                                             minimum_probability=0.5, cpu_workers=4, engine_root=root)
        self.assertEqual(result["sites"], sites)
        self.assertEqual(result["site_count"], 211)
        self.assertFalse(result["truncated"])
        self.assertEqual(result["protein_residue_count"], 5001)

    def test_timeout_kills_worker_session_logs_failure_and_cleans_temporary_files(self):
        root = self.fake_engine()
        process = SimpleNamespace(pid=12345, wait=Mock(side_effect=[subprocess.TimeoutExpired("worker", 1), -9]))
        seen = []

        def launch(*args, **kwargs):
            work = Path(kwargs["cwd"])
            seen.append(work)
            (work / "worker.log").write_text("discard prefix\n" + "x" * 40000 + "\ninput a failed")
            (work / "database.log").write_text("database detail")
            return process

        output = io.StringIO()
        with patch.object(engine, "_load_runner", return_value=Mock()), \
                patch.object(engine.subprocess, "Popen", side_effect=launch), \
                patch.object(engine.os, "killpg") as kill, redirect_stderr(output), \
                self.assertRaisesRegex(RuntimeError, "exceeded the 1-second limit"):
            engine.run_predictions([self.item], histidine_rotamers=True, minimum_probability=0,
                                   max_sites=100, cpu_workers=2, engine_root=root, timeout=1)
        kill.assert_called_once_with(12345, engine.signal.SIGKILL)
        self.assertFalse(seen[0].exists())
        self.assertIn("input a failed", output.getvalue())
        self.assertIn("database detail", output.getvalue())
        self.assertNotIn("discard prefix", output.getvalue())
        self.assertLess(len(output.getvalue()), 33500)


if __name__ == "__main__":
    unittest.main()
