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
DETAILED_FIELDS = ("structure_id", "source", "site_rank", "probability", "x_angstrom",
                   "y_angstrom", "z_angstrom", "coordinating_residues", "residue_combination",
                   "template_id", "template_rmsd_angstrom")


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
        with (self.output / "sites.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=DETAILED_FIELDS)
            writer.writeheader()
            for result in results:
                for rank, site in enumerate(result["sites"], 1):
                    x, y, z = site["coordinates_angstrom"]
                    source = result["source"]
                    writer.writerow({"structure_id": result["key"],
                                     "source": "'" + source if source.startswith(("=", "+", "-", "@")) else source,
                                     "site_rank": rank, "probability": site["probability"],
                                     "x_angstrom": x, "y_angstrom": y, "z_angstrom": z,
                                     "coordinating_residues": json.dumps(site["coordinating_residues"]),
                                     "residue_combination": site.get("residue_combination", "{H,H,H,E}"),
                                     "template_id": site["template_id"],
                                     "template_rmsd_angstrom": site["template_rmsd_angstrom"]})
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
        self.assertFalse(list(self.output.glob(".packaging-*")))

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
        expected_public = {"00_Summary.csv", "sites.csv", "CITATION.txt", "prediction.json",
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
        self.assertEqual((self.output / "sites.csv").read_bytes(), originals["sites.csv"])
        self.assertEqual(metadata["predicted_sites_file"], "sites.csv")
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
        sites = self.read_csv("sites.csv")
        self.assertEqual(tuple(sites[0]), DETAILED_FIELDS)
        self.assertEqual(len(sites), 2)
        self.assertEqual(json.loads(sites[0]["coordinating_residues"]), ["A_48", "B-2_36"])
        self.assertEqual(sites[0]["structure_id"], results[0]["key"])
        self.assertEqual(sites[0]["site_rank"], "1")
        self.assertEqual([float(sites[0][axis + "_angstrom"]) for axis in "xyz"],
                         results[0]["sites"][0]["coordinates_angstrom"])
        self.assertEqual(float(sites[0]["probability"]), results[0]["sites"][0]["probability"])
        self.assertEqual(float(sites[0]["template_rmsd_angstrom"]), results[0]["sites"][0]["template_rmsd_angstrom"])
        self.assertEqual(sites[0]["template_id"], "32")

    def test_zero_predictions_keep_viewer_session_and_header_only_site_table(self):
        results = self.fixture([{"sites": []}])
        metadata = packager.finalize_outputs(self.output, structures=results)
        mapping = metadata["structures"][results[0]["key"]]
        self.assertTrue((self.output / mapping["viewer_file"]).is_file())
        self.assertTrue((self.output / mapping["pymol_session"]).is_file())
        self.assertEqual(self.read_csv("sites.csv"), [])
        with (self.output / "sites.csv").open() as stream:
            self.assertEqual(next(csv.reader(stream)), list(DETAILED_FIELDS))

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
        self.assertEqual(self.read_csv("sites.csv")[0]["source"], "0_=SUM(1).pdb")
        self.assertEqual(metadata["structures"][results[0]["key"]]["display_name"], "=SUM(1)")
        self.assertEqual(json.loads((self.output / "prediction.json").read_text())["structures"][0]["source"], "0_=SUM(1).pdb")

    def test_unknown_root_files_directories_and_archives_are_untouched_and_not_archived(self):
        results = self.fixture([{}])
        extra = self.output / "details/nested/engine.log"
        extra.parent.mkdir(parents=True)
        extra.write_bytes(b"diagnostic details\x00\xff")
        (self.output / "other.zip").write_bytes(b"another useful artifact")
        unknown = {path: (path.stat().st_ino, path.read_bytes() if path.is_file() else None)
                   for path in (extra, extra.parent, extra.parent.parent, self.output / "other.zip")}
        packager.finalize_outputs(self.output, structures=results)
        for path, (inode, contents) in unknown.items():
            self.assertEqual(path.stat().st_ino, inode)
            if contents is not None:
                self.assertEqual(path.read_bytes(), contents)
        with zipfile.ZipFile(self.output / "ZincSight_results.zip") as archive:
            self.assertNotIn("details/nested/engine.log", archive.namelist())
            self.assertNotIn("other.zip", archive.namelist())

    def platform_files(self):
        values = {"output.log": b"before packaging\n", ".exit_code": b"running\n",
                  "settings.json": b'{"platform":"settings"}', "inputs/upload.pdb": b"platform upload",
                  ".platform_state/state.json": b'{"platform":"state"}'}
        for name, value in values.items():
            path = self.output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(value)
        paths = [self.output / name for name in values]
        paths.extend((self.output / "inputs", self.output / ".platform_state"))
        return values, {path: path.stat().st_ino for path in paths}

    def test_open_platform_log_and_exit_state_keep_their_inodes_and_paths_on_success(self):
        results = self.fixture([{}])
        values, inodes = self.platform_files()
        with (self.output / "output.log").open("ab") as log:
            packager.finalize_outputs(self.output, structures=results)
            log.write(b"after packaging through original open handle\n")
            log.flush()
        self.assertEqual((self.output / "output.log").read_bytes(),
                         values["output.log"] + b"after packaging through original open handle\n")
        for path, inode in inodes.items():
            self.assertEqual(path.stat().st_ino, inode, str(path))
        for name, contents in values.items():
            if name != "output.log":
                self.assertEqual((self.output / name).read_bytes(), contents)
        with zipfile.ZipFile(self.output / "ZincSight_results.zip") as archive:
            self.assertTrue(set(values).isdisjoint(archive.namelist()))
            self.assertFalse(any(name.startswith(("inputs/", ".platform_state/")) for name in archive.namelist()))

    def test_open_platform_log_and_service_roots_are_untouched_during_failed_publication(self):
        results = self.fixture([{}])
        values, inodes = self.platform_files()
        original = self.snapshot()
        replace = os.replace
        published = []
        def fail_after_one_public_file(source, destination):
            if Path(source).parent.name == "public":
                if published:
                    raise OSError("publication failure with platform files")
                published.append(Path(destination).name)
            return replace(source, destination)
        with (self.output / "output.log").open("ab") as log:
            with patch.object(packager.os, "replace", side_effect=fail_after_one_public_file):
                with self.assertRaisesRegex(OSError, "publication failure with platform files"):
                    packager.finalize_outputs(self.output, structures=results)
            self.assert_unchanged(original)
            log.write(b"after failed packaging through original open handle\n")
            log.flush()
        self.assertEqual((self.output / "output.log").read_bytes(),
                         values["output.log"] + b"after failed packaging through original open handle\n")
        for path, inode in inodes.items():
            self.assertEqual(path.stat().st_ino, inode, str(path))
        for name, contents in values.items():
            if name != "output.log":
                self.assertEqual((self.output / name).read_bytes(), contents)

    def test_invalid_or_missing_artifacts_leave_original_outputs_unchanged(self):
        cases = [lambda r: r[0]["artifacts"].pop("mmcif"),
                 lambda r: r[0]["artifacts"].update(pdb="../outside.pdb"),
                 lambda r: r[0]["artifacts"].update(pdb="/tmp/outside.pdb"),
                 lambda r: r[0]["artifacts"].update(pdb="structures/missing.pdb"),
                 lambda r: r[0]["artifacts"].update(pdb="output.log"),
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

    def test_symlink_output_or_app_owned_member_cannot_copy_external_files(self):
        results = self.fixture([{}])
        outside = self.output.parent / "private.txt"
        outside.write_text("private")
        linked = self.output / results[0]["artifacts"]["pdb"]
        linked.unlink()
        linked.symlink_to(outside)
        with self.assertRaises(ValueError):
            packager.finalize_outputs(self.output, structures=results)
        self.assertTrue(linked.is_symlink())
        alias = self.output.parent / "alias"
        alias.symlink_to(self.output, target_is_directory=True)
        with self.assertRaises(ValueError):
            packager.finalize_outputs(alias, structures=results)
        self.assertEqual(outside.read_text(), "private")

    def test_unknown_platform_symlink_is_untouched_and_not_followed(self):
        results = self.fixture([{}])
        outside = self.output.parent / "platform-data"
        outside.mkdir()
        (outside / "private.txt").write_text("private platform data")
        linked = self.output / "platform-linked"
        linked.symlink_to(outside, target_is_directory=True)
        inode = linked.lstat().st_ino
        packager.finalize_outputs(self.output, structures=results)
        self.assertTrue(linked.is_symlink())
        self.assertEqual(linked.lstat().st_ino, inode)
        with zipfile.ZipFile(self.output / "ZincSight_results.zip") as archive:
            self.assertFalse(any(name.startswith("platform-linked") for name in archive.namelist()))
        self.assertEqual((outside / "private.txt").read_text(), "private platform data")

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

    def test_directory_inode_and_permissions_survive_when_directory_rename_is_forbidden(self):
        results = self.fixture([{}])
        self.output.chmod(0o750)
        before = self.output.stat()
        replace = os.replace
        moves = []
        def forbid_output_rename(source, destination):
            source, destination = Path(source), Path(destination)
            if source == self.output or destination == self.output:
                raise PermissionError("mounted output directory cannot be renamed")
            self.assertTrue(source.is_relative_to(self.output))
            self.assertTrue(destination.is_relative_to(self.output))
            moves.append((source, destination))
            return replace(source, destination)
        with patch.object(packager.os, "replace", side_effect=forbid_output_rename):
            packager.finalize_outputs(self.output, structures=results)
        after = self.output.stat()
        self.assertEqual((before.st_ino, before.st_dev, before.st_mode, before.st_uid, before.st_gid),
                         (after.st_ino, after.st_dev, after.st_mode, after.st_uid, after.st_gid))
        self.assertTrue(moves)
        self.assertTrue((self.output / "ZincSight_results.zip").is_file())
        self.assertFalse(list(self.output.glob(".packaging-*")))

    def test_mid_file_publication_failure_restores_all_original_children(self):
        results = self.fixture([{}])
        original = self.snapshot()
        inode = self.output.stat().st_ino
        replace = os.replace
        published = []
        def fail_publication(source, destination):
            if Path(source).parent.name == "public":
                if len(published) == 2:
                    raise OSError("publication failure")
                published.append(Path(destination).name)
            return replace(source, destination)
        with patch.object(packager.os, "replace", side_effect=fail_publication):
            with self.assertRaisesRegex(OSError, "publication failure"):
                packager.finalize_outputs(self.output, structures=results)
        self.assertEqual(len(published), 2)
        self.assertEqual(self.output.stat().st_ino, inode)
        self.assert_unchanged(original)

    def test_mid_original_backup_failure_restores_previously_moved_children(self):
        results = self.fixture([{}])
        original = self.snapshot()
        replace = os.replace
        backed_up = []
        def fail_backup(source, destination):
            if Path(destination).parent.name == "original":
                if len(backed_up) == 2:
                    raise OSError("backup failure")
                backed_up.append(Path(source).name)
            return replace(source, destination)
        with patch.object(packager.os, "replace", side_effect=fail_backup):
            with self.assertRaisesRegex(OSError, "backup failure"):
                packager.finalize_outputs(self.output, structures=results)
        self.assertEqual(len(backed_up), 2)
        self.assert_unchanged(original)

    def test_failed_rollback_retains_original_files_for_recovery(self):
        results = self.fixture([{}])
        original = self.snapshot()
        replace = os.replace
        def fail_publication_and_one_restore(source, destination):
            source = Path(source)
            if source.parent.name == "public":
                raise OSError("publication failure")
            if source.parent.name == "original" and source.name == "CITATION.txt":
                raise OSError("restore failure")
            return replace(source, destination)
        with patch.object(packager.os, "replace", side_effect=fail_publication_and_one_restore):
            with self.assertRaisesRegex(RuntimeError, "rollback was incomplete"):
                packager.finalize_outputs(self.output, structures=results)
        temporary, = self.output.glob(".packaging-*")
        for name, contents in original.items():
            public = self.output / name
            saved = temporary / "original" / name
            self.assertTrue(public.is_file() or saved.is_file(), name)
            self.assertEqual((public if public.is_file() else saved).read_bytes(), contents)


if __name__ == "__main__":
    unittest.main()
