"""Public result layout and lossless packaging, using only fake artifact bytes."""
import copy
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

MODULE = Path(__file__).resolve().parents[1] / "integrations/tamarind/results_packager.py"
SPEC = importlib.util.spec_from_file_location("results_packager", MODULE)
packager = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(packager)


class ResultsPackagerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name) / "out"
        self.output.mkdir()

    def fixture(self, specs):
        results = []
        for number, spec in enumerate(specs, 1):
            key = f"s{number:04d}_input"
            sites = spec.get("sites", [{"probability": 0.91526484,
                                      "coordinating_residues": ["A_48", "B-2_36"],
                                      "coordinates_angstrom": [1.123456789, -2.0, 3.0],
                                      "template_id": "32", "template_rmsd_angstrom": 0.7248545171476303}])
            artifacts = {"mmcif": f"structures/{key}.cif", "pymol_session": f"structures/{key}.pse",
                         "site_annotations": f"structures/{key}_sites.json"}
            if spec.get("pdb", True):
                artifacts["pdb"] = f"structures/{key}.pdb"
            result = {"key": key, "source": spec.get("source", "0_1KLS.pdb"), "artifacts": artifacts,
                      "sites": sites, "site_count": spec.get("site_count", len(sites)),
                      "truncated": spec.get("truncated", False),
                      "input_provenance": spec.get("provenance", {}),
                      "prediction_json": f"predictions/{key}.json"}
            for kind, name in artifacts.items():
                path = self.output / name
                path.parent.mkdir(exist_ok=True)
                path.write_bytes(f"{kind} details for {key}\n".encode())
            prediction = self.output / result["prediction_json"]
            prediction.parent.mkdir(exist_ok=True)
            prediction.write_text(json.dumps(result, ensure_ascii=False))
            site_csv = self.output / "sites" / f"{key}.csv"
            site_csv.parent.mkdir(exist_ok=True)
            site_csv.write_text(f"full scientific site details for {key}\n")
            results.append(result)
        for name, text in {"CITATION.txt": "ZincSight DOI:10.1002/pro.70350\n",
                           "CITATION.bib": "@article{full-bibtex}\n", "citation.json": "{\"doi\":\"full-details\"}",
                           "structures.csv": "original structures CSV\n", "sites.csv": "original full sites CSV\n",
                           "RESULTS_README.txt": "full instructions\n",
                           "prediction.json": json.dumps({"structures": results}, ensure_ascii=False)}.items():
            (self.output / name).write_text(text)
        (self.output / "zincsight_results.zip").write_bytes(b"outdated previous bundle")
        return results

    def snapshot(self):
        return {p.relative_to(self.output).as_posix(): p.read_bytes()
                for p in self.output.rglob("*") if p.is_file()}

    def read_csv(self, name):
        with (self.output / name).open(newline="") as stream:
            return list(csv.DictReader(stream))

    def assert_unchanged(self, before):
        self.assertEqual(self.snapshot(), before)
        self.assertFalse(list(self.output.parent.glob(".out.packaging-*")))

    def test_compact_public_layout_and_complete_original_bundle(self):
        results = self.fixture([{}, {"source": "nested/2A0S-assembly1.cif", "pdb": False,
                                     "truncated": True, "site_count": 36}])
        originals = self.snapshot()
        input_copy = copy.deepcopy(results)
        metadata = packager.finalize_outputs(self.output, structures=results)
        self.assertEqual(results, input_copy)
        mappings = metadata["structures"]
        self.assertEqual(metadata["bundle_path"], "ZincSight_results.zip")
        self.assertEqual(mappings[results[0]["key"]]["viewer_file"], "01_1KLS_ZincSight.pdb")
        self.assertEqual(mappings[results[1]["key"]]["viewer_file"], "02_2A0S-assembly1_ZincSight.cif")
        expected_public = {"00_Summary.csv", "01_Predicted_sites.csv", "CITATION.txt", "prediction.json",
                           "ZincSight_results.zip"}
        for result in results:
            mapping = mappings[result["key"]]
            expected_public.update((mapping["viewer_file"], mapping["pymol_session"]))
            viewer_kind = "pdb" if "pdb" in result["artifacts"] else "mmcif"
            self.assertEqual((self.output / mapping["viewer_file"]).read_bytes(),
                             originals[result["artifacts"][viewer_kind]])
            self.assertEqual((self.output / mapping["pymol_session"]).read_bytes(),
                             originals[result["artifacts"]["pymol_session"]])
            self.assertEqual(mapping["bundle_members"], {**result["artifacts"],
                             "prediction_json": result["prediction_json"], "sites_csv": f"sites/{result['key']}.csv"})
        self.assertEqual(set(self.snapshot()), expected_public)
        self.assertEqual((self.output / "prediction.json").read_bytes(), originals["prediction.json"])
        with zipfile.ZipFile(self.output / metadata["bundle_path"]) as archive:
            names = archive.namelist()
            self.assertEqual(len(names), len(set(names)))
            self.assertNotIn("zincsight_results.zip", names)
            self.assertNotIn("ZincSight_results.zip", names)
            for name, contents in originals.items():
                if name != "zincsight_results.zip":
                    self.assertEqual(archive.read(name), contents, name)
            self.assertEqual(json.loads(archive.read("public_manifest.json")), metadata)
            old_batch = json.loads(archive.read("prediction.json"))
            for result in old_batch["structures"]:
                self.assertEqual(json.loads(archive.read(result["prediction_json"])), result)
                self.assertTrue(all(path in names for path in result["artifacts"].values()))

    def test_small_tables_preserve_numeric_values_and_show_truncation(self):
        results = self.fixture([{}, {"sites": [], "source": "2_empty.pdb"},
                                {"source": "3_full.cif", "truncated": True, "site_count": 36}])
        packager.finalize_outputs(self.output, structures=results)
        summary = self.read_csv("00_Summary.csv")
        self.assertEqual(tuple(summary[0]), packager.SUMMARY_FIELDS)
        self.assertEqual([row["status"] for row in summary], ["Complete", "No sites predicted", "Truncated"])
        self.assertEqual([row["sites"] for row in summary], ["1", "0", "36"])
        self.assertEqual(summary[1]["top_probability"], "")
        self.assertEqual(float(summary[0]["top_probability"]), results[0]["sites"][0]["probability"])
        sites = self.read_csv("01_Predicted_sites.csv")
        self.assertEqual(tuple(sites[0]), packager.SITE_FIELDS)
        self.assertEqual(len(sites), 2)
        self.assertEqual(sites[0]["coordinating_residues"], "A_48, B-2_36")
        self.assertEqual(float(sites[0]["probability"]), results[0]["sites"][0]["probability"])
        self.assertEqual(float(sites[0]["template_rmsd_angstrom"]), results[0]["sites"][0]["template_rmsd_angstrom"])
        self.assertEqual(sites[0]["template_id"], "32")

    def test_zero_predictions_keep_viewer_session_and_header_only_site_table(self):
        results = self.fixture([{"sites": []}])
        metadata = packager.finalize_outputs(self.output, structures=results)
        mapping = metadata["structures"][results[0]["key"]]
        self.assertTrue((self.output / mapping["viewer_file"]).is_file())
        self.assertTrue((self.output / mapping["pymol_session"]).is_file())
        self.assertEqual(self.read_csv("01_Predicted_sites.csv"), [])
        with (self.output / "01_Predicted_sites.csv").open() as stream:
            self.assertEqual(next(csv.reader(stream)), list(packager.SITE_FIELDS))

    def test_labels_preserve_archive_names_strip_only_one_upload_prefix_and_avoid_collisions(self):
        results = self.fixture([{"source": "01_02_same.pdb"}, {"source": "01_02_same.pdb"},
                                {"source": "nested/01_02_same.cif", "provenance": {"archive": "batch.zip"}},
                                {"source": "folder\\0_My weird 🧬 name.PDB.GZ"}])
        metadata = packager.finalize_outputs(self.output, structures=results)
        maps = list(metadata["structures"].values())
        self.assertEqual([m["display_name"] for m in maps], ["02_same", "02_same", "01_02_same", "My weird 🧬 name"])
        self.assertEqual(maps[0]["viewer_file"], "01_02_same_ZincSight.pdb")
        self.assertEqual(maps[1]["viewer_file"], "02_02_same_ZincSight.pdb")
        self.assertEqual(maps[3]["viewer_file"], "04_My_weird_name_ZincSight.pdb")
        self.assertEqual(len({m["viewer_file"] for m in maps}), 4)
        self.assertEqual(json.loads((self.output / "prediction.json").read_text())["structures"], results)

    def test_display_formula_escape_does_not_change_original_json(self):
        results = self.fixture([{"source": "0_=SUM(1).pdb"}])
        metadata = packager.finalize_outputs(self.output, structures=results)
        self.assertEqual(self.read_csv("00_Summary.csv")[0]["structure"], "'=SUM(1)")
        self.assertEqual(self.read_csv("01_Predicted_sites.csv")[0]["structure"], "'=SUM(1)")
        self.assertEqual(metadata["structures"][results[0]["key"]]["display_name"], "=SUM(1)")
        self.assertEqual(json.loads((self.output / "prediction.json").read_text())["structures"][0]["source"], "0_=SUM(1).pdb")

    def test_unexpected_detail_files_and_other_zip_are_retained_in_bundle(self):
        results = self.fixture([{}])
        extra = self.output / "details/nested/engine.log"
        extra.parent.mkdir(parents=True)
        extra.write_bytes(b"diagnostic details\x00\xff")
        (self.output / "other.zip").write_bytes(b"another useful artifact")
        packager.finalize_outputs(self.output, structures=results)
        self.assertFalse((self.output / "details").exists())
        self.assertFalse((self.output / "other.zip").exists())
        with zipfile.ZipFile(self.output / "ZincSight_results.zip") as archive:
            self.assertEqual(archive.read("details/nested/engine.log"), b"diagnostic details\x00\xff")
            self.assertEqual(archive.read("other.zip"), b"another useful artifact")

    def test_invalid_or_missing_artifacts_leave_original_outputs_unchanged(self):
        cases = [lambda r: r[0]["artifacts"].pop("mmcif"),
                 lambda r: r[0]["artifacts"].update(pdb="../outside.pdb"),
                 lambda r: r[0]["artifacts"].update(pdb="/tmp/outside.pdb"),
                 lambda r: r[0]["artifacts"].update(pdb="structures/missing.pdb"),
                 lambda r: r[0].update(key="../bad"),
                 lambda r: r[0].update(sites=[{"probability": float("nan"), "coordinating_residues": []}])]
        results = self.fixture([{}])
        original = self.snapshot()
        for change in cases:
            with self.subTest(change=change):
                invalid = copy.deepcopy(results)
                change(invalid)
                with self.assertRaises(ValueError):
                    packager.finalize_outputs(self.output, structures=invalid)
                self.assert_unchanged(original)
        with self.assertRaises(ValueError):
            packager.finalize_outputs(self.output, structures=[])
        self.assert_unchanged(original)

    def test_duplicate_keys_and_reserved_names_fail_without_changes(self):
        results = self.fixture([{}])
        original = self.snapshot()
        with self.assertRaises(ValueError):
            packager.finalize_outputs(self.output, structures=results + results)
        self.assert_unchanged(original)
        (self.output / "00_Summary.csv").write_text("already exists")
        original = self.snapshot()
        with self.assertRaises(ValueError):
            packager.finalize_outputs(self.output, structures=results)
        self.assert_unchanged(original)

    def test_symlink_output_or_member_cannot_copy_external_files(self):
        results = self.fixture([{}])
        outside = self.output.parent / "private.txt"
        outside.write_text("private")
        (self.output / "linked").symlink_to(outside)
        with self.assertRaises(ValueError):
            packager.finalize_outputs(self.output, structures=results)
        self.assertTrue((self.output / "linked").is_symlink())
        alias = self.output.parent / "alias"
        alias.symlink_to(self.output, target_is_directory=True)
        with self.assertRaises(ValueError):
            packager.finalize_outputs(alias, structures=results)
        self.assertEqual(outside.read_text(), "private")

    def test_staging_copy_and_zip_failures_leave_original_outputs_unchanged(self):
        results = self.fixture([{}])
        original = self.snapshot()
        with patch.object(packager.shutil, "copy2", side_effect=OSError("copy failure")):
            with self.assertRaisesRegex(OSError, "copy failure"):
                packager.finalize_outputs(self.output, structures=results)
        self.assert_unchanged(original)
        with patch.object(packager.zipfile.ZipFile, "write", side_effect=OSError("ZIP failure")):
            with self.assertRaisesRegex(OSError, "ZIP failure"):
                packager.finalize_outputs(self.output, structures=results)
        self.assert_unchanged(original)

    def test_publication_failure_restores_original_directory(self):
        results = self.fixture([{}])
        original = self.snapshot()
        replace = os.replace
        def fail_publication(source, destination):
            if Path(source).name == "public":
                raise OSError("publication failure")
            return replace(source, destination)
        with patch.object(packager.os, "replace", side_effect=fail_publication):
            with self.assertRaisesRegex(OSError, "publication failure"):
                packager.finalize_outputs(self.output, structures=results)
        self.assert_unchanged(original)


if __name__ == "__main__":
    unittest.main()
