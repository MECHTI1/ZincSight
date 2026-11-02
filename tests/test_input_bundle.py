"""Bounded upload/archive reader checks without scientific dependencies."""
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
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('tamarind_input_bundle', ROOT / 'integrations/tamarind/input_bundle.py')
bundle_reader = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bundle_reader)
PDB = b'ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n'
CIF = b'data_input\n_atom_site.id 1\n'


class InputBundleTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.app = Path(temporary.name)
        self.inputs = self.app / 'inputs'
        self.inputs.mkdir()

    def write(self, name, raw=PDB):
        path = self.inputs / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
        return path

    def zip(self, name, entries):
        path = self.inputs / name
        with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as archive:
            for member, raw in entries:
                archive.writestr(member, raw)
        return path

    def collect(self, path=None, **limits):
        environment = {} if path is None else {'batch_archive': str(path)}
        return bundle_reader.collect_structures(self.app, environment, **limits)

    def test_pdb_ent_cif_mmcif_and_uppercase_extensions_are_collected(self):
        for name, raw in [('a.PDB', PDB), ('b.ent', PDB), ('c.CIF', CIF), ('d.mmcif', CIF)]:
            self.write(name, raw)
        structures = self.collect()
        self.assertEqual([item['file_format'] for item in structures], ['pdb', 'pdb', 'cif', 'cif'])
        self.assertEqual([item['source'] for item in structures], ['a.PDB', 'b.ent', 'c.CIF', 'd.mmcif'])
        self.assertEqual(structures[0]['input_sha256'], hashlib.sha256(PDB).hexdigest())

    def test_explicit_multiple_paths_and_legacy_input_are_resolved_and_deduplicated(self):
        first = self.write('first.pdb')
        second = self.write('nested/second.cif', CIF)
        environment = {'structures': json.dumps(['first.pdb', str(first), 'inputs/nested/second.cif'])}
        structures = bundle_reader.collect_structures(self.app, environment)
        self.assertEqual(len(structures), 2)
        self.assertEqual([item['source'] for item in structures], ['first.pdb', 'second.cif'])
        self.assertEqual(bundle_reader.staged_paths(self.app, {'structure': 'inputs/first.pdb'}), [first])
        self.assertEqual(bundle_reader.staged_paths(self.app, {'structures': [str(second)]}), [second])

    def test_inputs_json_is_supported_without_environment_path_values(self):
        self.write('selected.pdb')
        self.write('ignored.cif', CIF)
        (self.app / 'inputs.json').write_text(json.dumps({'structures': ['selected.pdb']}))
        self.assertEqual([item['source'] for item in self.collect()], ['selected.pdb'])

    def test_gzip_preserves_decompressed_hash_and_compression_provenance(self):
        compressed = gzip.compress(PDB)
        path = self.write('sample.pdb.gz', compressed)
        structure = self.collect(path)[0]
        self.assertEqual(structure['structure_text'], PDB.decode())
        self.assertEqual(structure['source'], 'sample.pdb')
        self.assertEqual(structure['input_bytes'], len(PDB))
        self.assertEqual(structure['input_sha256'], hashlib.sha256(PDB).hexdigest())
        self.assertEqual(structure['provenance']['compressed_sha256'], hashlib.sha256(compressed).hexdigest())

    def test_zip_handles_nested_files_compressed_members_and_duplicate_basenames(self):
        path = self.zip('batch.zip', [('a/sample.pdb', PDB), ('b/sample.cif', CIF),
                                      ('c/third.ent.gz', gzip.compress(PDB)), ('README.txt', b'ignored')])
        structures = self.collect(path)
        self.assertEqual(len(structures), 3)
        self.assertEqual(len({item['key'] for item in structures}), 3)
        self.assertEqual([item['source'] for item in structures], ['a/sample.pdb', 'b/sample.cif', 'c/third.ent'])
        self.assertEqual(structures[2]['provenance']['archive'], 'batch.zip')
        self.assertFalse((self.app / 'a').exists())

    def test_tar_tar_gz_and_tgz_are_read_without_extracting_members(self):
        for suffix, mode in [('tar', 'w'), ('tar.gz', 'w:gz'), ('tgz', 'w:gz')]:
            path = self.inputs / ('batch.' + suffix)
            with tarfile.open(path, mode) as archive:
                directory = tarfile.TarInfo('folder')
                directory.type = tarfile.DIRTYPE
                archive.addfile(directory)
                info = tarfile.TarInfo('folder/protein.pdb')
                info.size = len(PDB)
                archive.addfile(info, io.BytesIO(PDB))
            with self.subTest(suffix=suffix):
                result = self.collect(path)
                self.assertEqual(result[0]['source'], 'folder/protein.pdb')
                self.assertEqual(result[0]['structure_text'], PDB.decode())
        self.assertFalse((self.app / 'folder').exists())

    def test_prepared_manifest_preserves_identifier_and_checks_exact_bytes(self):
        digest = hashlib.sha256(PDB).hexdigest()
        manifest = {'structures': [{'path': 'structures/one.pdb', 'sha256': digest,
                                    'identifier': 'AF-P05067-F1-model_v6', 'original_name': 'P05067'}]}
        path = self.zip('prepared.zip', [('manifest.json', json.dumps(manifest)), ('structures/one.pdb', PDB)])
        item = self.collect(path)[0]
        self.assertEqual(item['source'], 'P05067')
        self.assertEqual(item['provenance']['identifier'], 'AF-P05067-F1-model_v6')
        self.assertEqual(item['input_sha256'], digest)
        manifest['structures'][0]['sha256'] = '0' * 64
        path = self.zip('corrupt.zip', [('manifest.json', json.dumps(manifest)), ('structures/one.pdb', PDB)])
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            self.collect(path)

    def test_manifest_and_members_use_the_same_normalized_archive_names(self):
        manifest = {'structures': [{'path': 'folder/one.pdb'}]}
        path = self.zip('normalized.zip', [('./manifest.json', json.dumps(manifest)), ('folder\\one.pdb', PDB)])
        self.assertEqual(self.collect(path)[0]['source'], 'folder/one.pdb')

    def test_invalid_manifests_missing_members_and_unsupported_entries_are_rejected(self):
        cases = [([], []), ({'structures': 'bad'}, []), ({'structures': [None]}, []),
                 ({'structures': [{'path': 'missing.pdb'}]}, []),
                 ({'structures': [{'path': 'notes.txt'}]}, [('notes.txt', b'notes')]),
                 ({'structures': [{'path': 'one.pdb'}, {'path': 'one.pdb'}]}, [('one.pdb', PDB)])]
        for number, (manifest, extra) in enumerate(cases):
            path = self.zip(f'invalid-{number}.zip', [('manifest.json', json.dumps(manifest))] + extra)
            with self.subTest(manifest=manifest), self.assertRaises(ValueError):
                self.collect(path)

    def test_zip_traversal_duplicate_normalized_names_and_symlinks_are_rejected(self):
        for name in ('../escape.pdb', '/absolute.pdb', 'C:\\escape.pdb'):
            path = self.zip('unsafe.zip', [(name, PDB)])
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'Archive paths'):
                self.collect(path)
        path = self.zip('duplicates.zip', [('one.pdb', PDB), ('./one.pdb', PDB)])
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.collect(path)
        path = self.inputs / 'symlink.zip'
        info = zipfile.ZipInfo('one.pdb')
        info.create_system = 3
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        with zipfile.ZipFile(path, 'w') as archive:
            archive.writestr(info, '../outside.pdb')
        with self.assertRaisesRegex(ValueError, 'Symbolic links'):
            self.collect(path)

    def test_tar_traversal_and_links_are_rejected(self):
        for kind in ('traversal', 'symlink', 'hardlink'):
            path = self.inputs / (kind + '.tar')
            with tarfile.open(path, 'w') as archive:
                member = tarfile.TarInfo('../escape.pdb' if kind == 'traversal' else 'one.pdb')
                if kind == 'traversal':
                    member.size = len(PDB)
                    archive.addfile(member, io.BytesIO(PDB))
                else:
                    member.type = tarfile.SYMTYPE if kind == 'symlink' else tarfile.LNKTYPE
                    member.linkname = '../outside.pdb'
                    archive.addfile(member)
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                self.collect(path)

    def test_staged_files_cannot_escape_inputs_even_through_symlinks(self):
        outside = self.app / 'outside.pdb'
        outside.write_bytes(PDB)
        link = self.inputs / 'linked.pdb'
        link.symlink_to(outside)
        for value in (str(outside), '../outside.pdb', 'linked.pdb', 'missing.pdb'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'staged files'):
                bundle_reader.collect_structures(self.app, {'structure': value})
        self.assertEqual(outside.read_bytes(), PDB)

    def test_structure_count_per_file_total_size_and_gzip_expansion_are_bounded(self):
        path = self.zip('two.zip', [('one.pdb', b'123456'), ('two.pdb', b'123456')])
        for limits in ({'max_structures': 1}, {'max_file_bytes': 5},
                       {'max_file_bytes': 10, 'max_total_bytes': 11}):
            with self.subTest(limits=limits), self.assertRaises(ValueError):
                self.collect(path, **limits)
        path = self.write('expanded.pdb.gz', gzip.compress(b'x' * 1000))
        with self.assertRaisesRegex(ValueError, 'uncompressed size'):
            self.collect(path, max_file_bytes=50)
        path = self.zip('too-many.zip', [(f'entry-{i}.txt', b'x') for i in range(25)])
        with self.assertRaisesRegex(ValueError, 'too many entries'):
            self.collect(path, max_structures=1)

    def test_empty_binary_nul_and_unsupported_inputs_are_rejected(self):
        for raw in (b'', b'\xff', b'ATOM\x00'):
            path = self.write('invalid.pdb', raw)
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                self.collect(path)
        path = self.write('notes.txt', b'not a structure')
        with self.assertRaisesRegex(ValueError, 'No supported'):
            self.collect(path)

    def test_invalid_upload_lists_and_inputs_json_fail_clearly(self):
        self.write('valid.pdb')
        for value in ('[bad JSON', '{"not": "a filename"}', 7, [None]):
            with self.subTest(value=value), self.assertRaises(ValueError):
                bundle_reader.collect_structures(self.app, {'structures': value})
        (self.app / 'inputs.json').write_text('[]')
        with self.assertRaisesRegex(ValueError, 'named input values'):
            self.collect()


if __name__ == '__main__':
    unittest.main()
