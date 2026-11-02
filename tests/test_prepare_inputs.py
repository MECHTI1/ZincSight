"""Network-free preparation, provenance, and archive-boundary regressions."""
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import stat
import tarfile
import tempfile
import unittest
from unittest.mock import Mock, patch
import urllib.error
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("prepare_inputs", ROOT / "integrations/tamarind/prepare_inputs.py")
prepare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prepare)
PDB = b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \nEND\n"
CIF = b"data_test\nloop_\n_atom_site.id\n1\n"


class Response(io.BytesIO):
    def __init__(self, data=PDB, url="https://files.rcsb.org/download/1KLS.cif", headers=None):
        super().__init__(data)
        self.url, self.headers = url, headers or {}

    def geturl(self):
        return self.url


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def write(self, name, data=PDB):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        return path

    def test_identifier_templates_preserve_explicit_fragments_versions_and_ted_domains(self):
        cases = {
            "1kls": ("1KLS", "https://files.rcsb.org/download/1KLS.cif", "cif"),
            "2a0s-assembly1": ("2A0S-assembly1", "https://files.rcsb.org/download/2A0S-assembly1.cif", "cif"),
            "p12345": ("AF-P12345-F1-model_v6", "https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif", "cif"),
            "AF-P12345-F3-model_v4": ("AF-P12345-F3-model_v4", "https://alphafold.ebi.ac.uk/files/AF-P12345-F3-model_v4.cif", "cif"),
            "MGYP002537940442": ("MGYP002537940442", "https://api.esmatlas.com/fetchPredictedStructure/MGYP002537940442.pdb", "pdb"),
            "AF-A0A287EUI1-F1-model_v4_TED01": ("AF-A0A287EUI1-F1-model_v4_TED01", "https://ted.cathdb.info/api/v1/files/AF-A0A287EUI1-F1-model_v4_TED01.pdb", "pdb"),
        }
        for identifier, expected in cases.items():
            with self.subTest(identifier=identifier):
                result = prepare.resolve_identifier(identifier)
                self.assertEqual((result["canonical_id"], result["source_url"], result["file_format"]), expected)
                prepare.validate_provider_url(result["source_url"])
        for identifier in ("https://user:secret@example.com/file", "../../private", "1KLS-assembly0", "AF-P12345-F0-model_v4", "AF-P12345-F1-model_v4_TED00"):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError) as error:
                prepare.resolve_identifier(identifier)
            self.assertNotIn("secret", str(error.exception))

    def test_colab_text_and_id_list_formats(self):
        self.assertEqual(prepare.parse_ids("1KLS, P12345\n1KLS\tMGYP002537940442 # note\n"), ["1KLS", "P12345", "MGYP002537940442"])
        self.assertEqual(prepare.read_id_list(b"structure_id,description\n1KLS,example protein\nP12345,second protein\n", ".csv"), ["1KLS", "P12345"])
        self.assertEqual(prepare.read_id_list(b"accession\tnote\nP12345\tignore this\n", ".tsv"), ["P12345"])
        self.assertEqual(prepare.read_id_list(b'"1KLS","P12345"\n', ".csv"), ["1KLS", "P12345"])

    def test_downloaded_aliases_are_deduplicated_and_manifest_hashes_cover_archived_bytes(self):
        ids_file = self.write("ids.csv", b"id,note\nP12345,one\nAF-P12345-F1-model_v6,same model\n")
        downloader = Mock(return_value=CIF)
        output = self.root / "prepared.zip"
        manifest = prepare.prepare_inputs(output, ids="1kls,1KLS", id_files=[ids_file], downloader=downloader)
        self.assertEqual(downloader.call_count, 2)
        self.assertEqual(manifest["structure_count"], 2)
        self.assertFalse(manifest["requires_network_at_runtime"])
        self.assertEqual(manifest["structures"][1]["canonical_id"], "AF-P12345-F1-model_v6")
        with zipfile.ZipFile(output) as archive:
            self.assertEqual(json.loads(archive.read("manifest.json")), manifest)
            for entry in manifest["structures"]:
                data = archive.read(entry["path"])
                self.assertEqual(entry["sha256"], hashlib.sha256(data).hexdigest())
                self.assertEqual(entry["size_bytes"], len(data))
        self.assertEqual(json.loads(output.with_suffix(".manifest.json").read_text()), manifest)
        self.assertNotIn(str(self.root), json.dumps(manifest))

    def test_folder_gzip_and_archives_preserve_duplicate_basenames_without_overwriting(self):
        first = self.write("folder/a/model.pdb")
        second = self.write("folder/b/model.pdb", PDB.replace(b"0.000", b"1.000"))
        compressed = self.write("single.ent.gz", gzip.compress(PDB))
        archive = self.root / "structures.zip"
        with zipfile.ZipFile(archive, "w") as stream:
            stream.writestr("nested/model.mmcif", CIF)
            stream.writestr("README.md", b"Ignored documentation")
        tar_path = self.root / "more.tgz"
        with tarfile.open(tar_path, "w:gz") as stream:
            info = tarfile.TarInfo("other/model.pdb.gz")
            data = gzip.compress(PDB)
            info.size = len(data)
            stream.addfile(info, io.BytesIO(data))
        manifest = prepare.prepare_inputs(self.root / "batch.zip", sources=[first.parent.parent, compressed, archive, tar_path])
        self.assertEqual(manifest["structure_count"], 5)
        paths = [entry["path"] for entry in manifest["structures"]]
        self.assertEqual(len(set(paths)), 5)
        self.assertEqual([entry["original_name"] for entry in manifest["structures"]].count("model.pdb"), 2)
        self.assertEqual([entry["file_format"] for entry in manifest["structures"]], ["pdb", "pdb", "pdb", "cif", "pdb"])
        self.assertEqual(manifest["structures"][2]["sha256"], hashlib.sha256(PDB).hexdigest())
        self.assertEqual(manifest["structures"][3]["archive_member"], "nested/model.mmcif")

    def test_zip_paths_links_and_portable_path_collisions_are_rejected(self):
        for index, bad_name in enumerate(("../escape.pdb", "/absolute.pdb", "C:/file.pdb", "dir\\escape.pdb", "a/./file.pdb")):
            with self.subTest(name=bad_name):
                path = self.root / f"bad-{index}.zip"
                with zipfile.ZipFile(path, "w") as archive:
                    archive.writestr(bad_name, PDB)
                output = self.root / f"output-{index}.zip"
                with self.assertRaises(ValueError):
                    prepare.prepare_inputs(output, sources=[path])
                self.assertFalse(output.exists())
        link_archive = self.root / "link.zip"
        with zipfile.ZipFile(link_archive, "w") as archive:
            info = zipfile.ZipInfo("link.pdb")
            info.external_attr = (stat.S_IFLNK | 0o777) << 16
            archive.writestr(info, "elsewhere.pdb")
        with self.assertRaises(ValueError):
            prepare.prepare_inputs(self.root / "link-output.zip", sources=[link_archive])
        duplicate_archive = self.root / "duplicate.zip"
        with zipfile.ZipFile(duplicate_archive, "w") as archive:
            archive.writestr("A/model.pdb", PDB)
            archive.writestr("a/MODEL.pdb", PDB)
        with self.assertRaises(ValueError):
            prepare.prepare_inputs(self.root / "duplicate-output.zip", sources=[duplicate_archive])

    def test_tar_hardlinks_and_traversal_are_rejected_without_extraction(self):
        for index, name in enumerate(("../escape.pdb", "link.pdb")):
            archive = self.root / f"bad-{index}.tar"
            with tarfile.open(archive, "w") as stream:
                info = tarfile.TarInfo(name)
                if index:
                    info.type, info.linkname = tarfile.LNKTYPE, "target.pdb"
                    stream.addfile(info)
                else:
                    info.size = len(PDB)
                    stream.addfile(info, io.BytesIO(PDB))
            with self.assertRaises(ValueError):
                prepare.prepare_inputs(self.root / f"tar-output-{index}.zip", sources=[archive])
        self.assertFalse((self.root.parent / "escape.pdb").exists())

    def test_standard_tar_dot_prefix_is_normalized_but_traversal_is_not(self):
        archive = self.root / "standard.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            directory = tarfile.TarInfo(".")
            directory.type = tarfile.DIRTYPE
            stream.addfile(directory)
            info = tarfile.TarInfo("./nested/model.pdb")
            info.size = len(PDB)
            stream.addfile(info, io.BytesIO(PDB))
        manifest = prepare.prepare_inputs(self.root / "standard.zip", sources=[archive])
        self.assertEqual(manifest["structures"][0]["archive_member"], "nested/model.pdb")
        with self.assertRaises(ValueError):
            prepare.safe_archive_name("./../escape.pdb")

    def test_item_file_total_and_gzip_expansion_limits_fail_before_output(self):
        first, second = self.write("one.pdb"), self.write("two.pdb")
        cases = [
            ({"sources": [first, second], "max_structures": 1}),
            ({"sources": [first], "max_file_bytes": len(PDB) - 1}),
            ({"sources": [first, second], "max_total_bytes": len(PDB) * 2 - 1}),
            ({"sources": [self.write("large.pdb.gz", gzip.compress(PDB * 5))], "max_file_bytes": len(PDB) * 2}),
        ]
        for index, kwargs in enumerate(cases):
            output = self.root / f"bounded-{index}.zip"
            with self.subTest(index=index), self.assertRaises(ValueError):
                prepare.prepare_inputs(output, **kwargs)
            self.assertFalse(output.exists())
        downloader = Mock(return_value=CIF)
        with self.assertRaises(ValueError):
            prepare.prepare_inputs(self.root / "too-many-ids.zip", ids="1KLS,P12345", max_structures=1, downloader=downloader)
        downloader.assert_not_called()

    def test_archive_and_folder_entry_limits_are_enforced_even_for_unrelated_files(self):
        archive = self.root / "ignored.zip"
        with zipfile.ZipFile(archive, "w") as stream:
            stream.writestr("large.bin", b"x" * 100)
        with self.assertRaises(ValueError):
            prepare.prepare_inputs(self.root / "archive-limit.zip", sources=[archive], max_total_bytes=50)
        folder = self.root / "many"
        for index in range(11):
            self.write(f"many/{index}.md", b"ignore")
        with patch.object(prepare, "MAX_ARCHIVE_ENTRIES", 1), self.assertRaises(ValueError):
            prepare.prepare_inputs(self.root / "folder-limit.zip", sources=[folder])

    def test_failure_preserves_existing_outputs_and_successful_archives_are_deterministic(self):
        valid = self.write("input.pdb")
        first, second = self.root / "first.zip", self.root / "second.zip"
        prepare.prepare_inputs(first, sources=[valid])
        prepare.prepare_inputs(second, sources=[valid])
        self.assertEqual(first.read_bytes(), second.read_bytes())
        previous_zip = first.read_bytes()
        previous_manifest = first.with_suffix(".manifest.json").read_bytes()
        bad = self.write("invalid.pdb", b"not a structure")
        with self.assertRaises(ValueError):
            prepare.prepare_inputs(first, sources=[bad], overwrite=True)
        self.assertEqual(first.read_bytes(), previous_zip)
        self.assertEqual(first.with_suffix(".manifest.json").read_bytes(), previous_manifest)
        with self.assertRaises(FileExistsError):
            prepare.prepare_inputs(first, sources=[valid])
        self.assertEqual(list(self.root.glob(".zincsight-*")), [])

    def test_cli_accepts_local_files_without_network_or_credentials(self):
        path = self.write("sample.pdb")
        output = self.root / "cli.zip"
        with contextlib.redirect_stdout(io.StringIO()), patch.dict("os.environ", {"TAMARIND_API_KEY": "must-not-be-read"}):
            self.assertEqual(prepare.main([str(path), "--output", str(output)]), 0)
        self.assertNotIn(b"must-not-be-read", output.with_suffix(".manifest.json").read_bytes())


class DownloadTests(unittest.TestCase):
    URL = "https://files.rcsb.org/download/1KLS.cif"

    def test_provider_allowlist_applies_before_requests_and_to_redirects(self):
        opener = Mock()
        for url in ("http://files.rcsb.org/download/1KLS.cif", "https://files.rcsb.org.evil.example/download/1KLS.cif",
                    "https://user:password@files.rcsb.org/download/1KLS.cif", self.URL + "?token=secret"):
            with self.subTest(url=url), self.assertRaises(ValueError):
                prepare.download_structure(url, opener=opener)
        opener.open.assert_not_called()
        with self.assertRaises(ValueError):
            prepare.ProviderRedirectHandler().redirect_request(urllib.request.Request(self.URL), None, 302, "Found", {}, "https://example.com/file")
        opener.open.return_value = Response(CIF, url="https://example.com/redirect")
        with self.assertRaises(ValueError):
            prepare.download_structure(self.URL, opener=opener)

    def test_transient_errors_have_bounded_retries_and_hard_http_errors_do_not_retry(self):
        opener, sleep = Mock(), Mock()
        opener.open.side_effect = [urllib.error.URLError("temporary"), Response(CIF, self.URL)]
        self.assertEqual(prepare.download_structure(self.URL, opener=opener, sleep=sleep), CIF)
        self.assertEqual(opener.open.call_count, 2)
        sleep.assert_called_once_with(0.25)
        opener.open.side_effect = urllib.error.HTTPError(self.URL, 404, "Not found", {}, None)
        opener.open.reset_mock()
        with self.assertRaisesRegex(ValueError, "404"):
            prepare.download_structure(self.URL, opener=opener, sleep=sleep)
        self.assertEqual(opener.open.call_count, 1)
        opener.open.side_effect = urllib.error.URLError("temporary")
        opener.open.reset_mock()
        with self.assertRaisesRegex(ValueError, "bounded retries"):
            prepare.download_structure(self.URL, opener=opener, sleep=Mock())
        self.assertEqual(opener.open.call_count, 3)

    def test_download_size_and_html_guards(self):
        for response in (Response(CIF, self.URL, {"Content-Length": "5000001"}),
                         Response(b"x" * 11, self.URL), Response(b"<html>error</html>", self.URL, {"Content-Type": "text/html"})):
            opener = Mock()
            opener.open.return_value = response
            with self.subTest(headers=response.headers), self.assertRaises(ValueError):
                prepare.download_structure(self.URL, max_bytes=10, opener=opener)


if __name__ == "__main__":
    unittest.main()
