"""Stdlib-only adapter and real result-packaging contracts; structures are checked in-image."""
import copy
import builtins
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('tamarind_adapter', ROOT / 'integrations/tamarind/predict.py')
adapter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(adapter)
SITE = {'probability': 0.91, 'coordinates_angstrom': [1.1, 2.2, 3.3],
        'coordinating_residues': ['A_3', 'A_5', 'A_9'],
        'residue_combination': 'EHH', 'template_id': 'example', 'template_rmsd_angstrom': 0.2}
PDB = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n'


class AdapterTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.app = Path(temporary.name)
        (self.app / 'inputs').mkdir()
        self.structure = self.app / 'inputs/sample.pdb'
        self.structure.write_text(PDB)
        self.calls = []
        self.inspected_inputs = []
        self.allocation = {
            'requested_cpu_workers': 8, 'cpu_workers': 2, 'declared_cpu_budget': 8,
            'declared_memory_gib': 24, 'detected_cpu_limit': 2, 'effective_memory_gib': 6,
            'memory_budget_gib_per_worker': 3, 'worker_limit_applied': True,
        }
        patcher = patch.object(adapter, 'worker_allocation', return_value=self.allocation)
        self.addCleanup(patcher.stop)
        patcher.start()

    def predictor(self, structures, **settings):
        self.calls.append((copy.deepcopy(structures), settings))
        self.inspected_inputs.append([Path(source['structure_path']).read_text()
                                      for source in structures])
        return [{'key': source['key'], 'sites': [copy.deepcopy(SITE)], 'site_count': 1,
                 'truncated': False, 'protein_residue_count': 3,
                 'histidine_rotamers': settings['histidine_rotamers'],
                 'minimum_probability': settings['minimum_probability']} for source in structures]

    def exporter(self, text, file_format, result, directory, stem):
        directory.mkdir(parents=True, exist_ok=True)
        paths = {}
        for key, suffix in [('pdb', '.pdb'), ('mmcif', '.cif'),
                            ('pymol_session', '.pse'), ('site_annotations', '_sites.json')]:
            path = directory / (stem + suffix)
            path.write_bytes(f'{stem}:{key}:export test artifact\n'.encode())
            paths[key] = path
        return {'artifact_paths': paths, 'warnings': [], 'predicted_marker_chain': 'Z'}

    def run_adapter(self, environment=None, predictor=None, exporter=None):
        return adapter.run(self.app, environment or {},
                           self.predictor if predictor is None else predictor,
                           self.exporter if exporter is None else exporter)

    def rows(self, filename):
        with (self.app / 'out' / filename).open(newline='') as stream:
            return list(csv.DictReader(stream))

    def bundled_bytes(self, filename):
        bundles = list((self.app / 'out').glob('*.zip'))
        self.assertEqual(len(bundles), 1)
        with zipfile.ZipFile(bundles[0]) as archive:
            return archive.read(filename)

    def bundled_rows(self, filename):
        return list(csv.DictReader(self.bundled_bytes(filename).decode().splitlines()))

    def assert_clean_public_output(self, result):
        self.assertEqual(result['output_layout'], 'summary_and_bundle')
        self.assertTrue((self.app / 'out' / result['result_bundle']).is_file())
        expected = {'00_Summary.csv', 'sites.csv', 'CITATION.txt',
                    'prediction.json', result['result_bundle']}
        for item in result['structures']:
            self.assertNotIn('prediction_json', item)
            self.assertIn('prediction_json_in_bundle', item)
            self.assertIn('bundle_members', item)
            artifacts = item['artifacts']
            self.assertEqual(artifacts['details_bundle'], result['result_bundle'])
            self.assertEqual(artifacts['structure'], artifacts['pdb'] if 'pdb' in artifacts else artifacts['mmcif'])
            for name in ('structure', 'pymol_session'):
                path = artifacts[name]
                self.assertEqual(Path(path).name, path)
                expected.add(path)
                self.assertTrue((self.app / 'out' / path).is_file())
            self.assertEqual((self.app / 'out' / artifacts['structure']).read_bytes(),
                             self.bundled_bytes(item['bundle_members']['pdb' if 'pdb' in artifacts else 'mmcif']))
            self.assertEqual((self.app / 'out' / artifacts['pymol_session']).read_bytes(),
                             self.bundled_bytes(item['bundle_members']['pymol_session']))
        self.assertEqual({path.name for path in (self.app / 'out').iterdir()}, expected)

    def test_single_input_settings_resources_and_download_contract(self):
        result = self.run_adapter({'structure': str(self.structure), 'histidine_rotamers': 'false',
                                   'minimum_probability': '0.4', 'cpu_workers': '8'})
        inputs, settings = self.calls[0]
        self.assertEqual(self.inspected_inputs[0][0], PDB)
        self.assertNotIn('structure_text', inputs[0])
        self.assertTrue(Path(inputs[0]['structure_path']).is_relative_to(Path('/tmp')))
        self.assertFalse(Path(inputs[0]['structure_path']).is_relative_to(self.app))
        self.assertEqual(settings, {'histidine_rotamers': False, 'minimum_probability': 0.4,
                                   'max_sites': None, 'cpu_workers': 2,
                                   'engine_root': adapter.ENGINE_ROOT, 'timeout': 7140})
        adapter.worker_allocation.assert_called_once()
        self.assertEqual(adapter.worker_allocation.call_args.args[0], 8)
        self.assertEqual(result['structure_count'], 1)
        self.assertEqual(result['sites'], [SITE])
        self.assertEqual(result['source'], 'sample.pdb')
        self.assertEqual(result['input_sha256'], hashlib.sha256(PDB.encode()).hexdigest())
        self.assertEqual(result['software_revision'], adapter.ENGINE_REVISION)
        self.assertEqual(result['integration_version'], adapter.INTEGRATION_VERSION)
        self.assertEqual(result['resources']['cpu_workers'], 2)
        self.assertEqual(result['database_transport'], 'unix_socket')
        self.assertNotIn('max_sites', adapter.read_parameters({}))
        self.assertNotIn('max_sites', result)
        self.assertNotIn('max_sites_per_structure', result)
        self.assertTrue(all('max_sites' not in item for item in result['structures']))
        self.assertEqual(json.loads((self.app / 'out/prediction.json').read_text()), result)
        self.assert_clean_public_output(result)
        self.assertNotIn('prediction_json', result)
        summary = self.rows('00_Summary.csv')[0]
        self.assertEqual(summary['structure'], result['display_name'])
        self.assertEqual(summary['structure_file'], result['artifacts']['structure'])
        self.assertEqual(summary['sites'], '1')
        self.assertEqual(float(summary['top_probability']), SITE['probability'])
        self.assertEqual(summary['status'], 'Complete')
        row = self.bundled_rows('sites.csv')[0]
        self.assertEqual(json.loads(row['coordinating_residues']), SITE['coordinating_residues'])
        self.assertEqual(float(row['z_angstrom']), 3.3)
        site = self.rows('sites.csv')[0]
        self.assertEqual(tuple(site), tuple(adapter.SITE_FIELDS))
        self.assertEqual(site['structure_id'], result['key'])
        self.assertEqual(site['source'], result['source'])
        self.assertEqual((site['site_rank'], float(site['probability'])), ('1', SITE['probability']))
        self.assertEqual([float(site[axis + '_angstrom']) for axis in 'xyz'], SITE['coordinates_angstrom'])
        self.assertEqual((self.app / 'out/sites.csv').read_bytes(), self.bundled_bytes('sites.csv'))
        for path in result['artifacts'].values():
            self.assertTrue((self.app / 'out' / path).is_file())
        per_input = json.loads(self.bundled_bytes(result['prediction_json_in_bundle']))
        self.assertEqual(per_input['sites'], [SITE])
        self.assertEqual(per_input['artifacts'],
                         {key: result['bundle_members'][key] for key in per_input['artifacts']})
        self.assertEqual(result['bundle_members']['prediction_json'], result['prediction_json_in_bundle'])
        with zipfile.ZipFile(self.app / 'out' / result['result_bundle']) as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn(result['bundle_members']['pymol_session'], archive.namelist())
            self.assertIn(result['bundle_members']['pdb'], archive.namelist())
            self.assertIn('CITATION.bib', archive.namelist())
            self.assertNotIn(result['result_bundle'], archive.namelist())
            original = json.loads(archive.read('prediction.json'))
            self.assertEqual(original['structures'][0]['artifacts'], per_input['artifacts'])

    def test_citation_is_zincsight_paper_in_all_machine_and_download_formats(self):
        result = self.run_adapter()
        self.assertEqual(result['reference']['doi'], '10.1002/pro.70350')
        self.assertIn('Mechtinger G', result['citation'])
        self.assertIn('Protein Science. 2025;34(11):e70350', result['citation'])
        self.assertIn('10.1002/pro.70350', (self.app / 'out/CITATION.txt').read_text())
        self.assertIn('doi = {10.1002/pro.70350}', self.bundled_bytes('CITATION.bib').decode())
        self.assertEqual(json.loads(self.bundled_bytes('citation.json')), result['reference'])
        self.assertFalse((self.app / 'out/CITATION.bib').exists())
        self.assertFalse((self.app / 'out/citation.json').exists())

    def test_mixed_upload_keeps_duplicate_basenames_distinct_and_zero_site_inputs_visible(self):
        archive = self.app / 'inputs/nested.zip'
        second_text = PDB.replace('0.000', '1.000')
        with zipfile.ZipFile(archive, 'w') as bundle:
            bundle.writestr('nested/sample.pdb', second_text)
        def mixed(inputs, **settings):
            results = self.predictor(inputs, **settings)
            results[1].update(sites=[], site_count=0)
            return results
        result = self.run_adapter({'structures': json.dumps(['sample.pdb', 'nested.zip'])}, mixed)
        self.assertEqual(result['structure_count'], 2)
        self.assertEqual(self.inspected_inputs, [[PDB, second_text]])
        self.assertNotIn('sites', result)
        self.assertEqual(result['returned_site_count'], 1)
        self.assertEqual([item['key'] for item in result['structures']], ['s0001_sample', 's0002_sample'])
        self.assert_clean_public_output(result)
        summaries = self.rows('00_Summary.csv')
        self.assertEqual(len(summaries), 2)
        self.assertEqual(summaries[1]['sites'], '0')
        self.assertEqual(summaries[1]['status'], 'No sites predicted')
        self.assertEqual(summaries[1]['top_probability'], '')
        self.assertNotEqual(summaries[0]['structure_file'], summaries[1]['structure_file'])
        self.assertTrue((self.app / 'out' / result['structures'][1]['artifacts']['pymol_session']).is_file())
        self.assertEqual(len(self.rows('sites.csv')), 1)
        self.assertEqual(len(self.bundled_rows('sites.csv')), 1)

    def test_zero_sites_still_exports_structure_and_session_with_header_only_sites(self):
        def empty(inputs, **settings):
            results = self.predictor(inputs, **settings)
            results[0].update(sites=[], site_count=0)
            return results
        result = self.run_adapter(predictor=empty)
        self.assertEqual(result['returned_site_count'], 0)
        self.assertFalse(result['truncated'])
        self.assert_clean_public_output(result)
        self.assertEqual(self.rows('sites.csv'), [])
        self.assertEqual(self.bundled_rows('sites.csv'), [])
        self.assertEqual(len(self.rows('00_Summary.csv')), 1)
        self.assertEqual(self.rows('00_Summary.csv')[0]['status'], 'No sites predicted')
        self.assertTrue((self.app / 'out' / result['artifacts']['pymol_session']).is_file())

    def test_all_sites_above_threshold_are_returned_without_a_hidden_per_structure_limit(self):
        (self.app / 'inputs/second.cif').write_text('data_x\n_atom_site.id 1\n')
        counts = (137, 211)
        def all_sites(inputs, **settings):
            self.assertIsNone(settings['max_sites'])
            results = self.predictor(inputs, **settings)
            for result, total in zip(results, counts):
                sites = [{**copy.deepcopy(SITE), 'probability': 0.5 + index / (2 * (total + 1)),
                          'coordinates_angstrom': [index + 0.123456789, index / 3, -index / 7],
                          'template_id': f'template-{index + 1}'} for index in range(total)]
                result.update(site_count=total, sites=sites, truncated=False)
            return results
        result = self.run_adapter({'minimum_probability': '0.4'}, all_sites)
        self.assertEqual((result['site_count'], result['returned_site_count']), (sum(counts), sum(counts)))
        self.assertFalse(result['truncated'])
        self.assertNotIn('max_sites_per_structure', result)
        self.assertTrue(all('max_sites' not in item for item in result['structures']))
        self.assertEqual([int(row['sites']) for row in self.rows('00_Summary.csv')], list(counts))
        self.assertEqual([row['status'] for row in self.rows('00_Summary.csv')], ['Complete', 'Complete'])
        site_rows = self.rows('sites.csv')
        self.assertEqual(len(site_rows), sum(counts))
        self.assertEqual(tuple(site_rows[0]), tuple(adapter.SITE_FIELDS))
        for item in result['structures']:
            rows = [row for row in site_rows if row['structure_id'] == item['key']]
            for rank, (row, site) in enumerate(zip(rows, item['sites']), 1):
                self.assertEqual(int(row['site_rank']), rank)
                self.assertEqual(float(row['probability']), site['probability'])
                self.assertEqual([float(row[axis + '_angstrom']) for axis in 'xyz'], site['coordinates_angstrom'])
                self.assertEqual(json.loads(row['coordinating_residues']), site['coordinating_residues'])
        for item, count in zip(result['structures'], counts):
            self.assertEqual(len(item['sites']), count)
            archived = json.loads(self.bundled_bytes(item['prediction_json_in_bundle']))
            self.assertEqual(archived['sites'], item['sites'])
            self.assertNotIn('max_sites', archived)
        self.assertEqual(len(self.bundled_rows('sites.csv')), sum(counts))

    def test_spreadsheet_labels_are_escaped_while_json_keeps_original_source(self):
        raw = PDB.encode()
        archive = self.app / 'inputs/bundle.zip'
        source_name = '=HYPERLINK(unsafe)'
        manifest = {'structures': [{'path': 'input.pdb', 'original_name': source_name,
                                      'sha256': hashlib.sha256(raw).hexdigest(), 'identifier': '1KLS'}]}
        with zipfile.ZipFile(archive, 'w') as bundle:
            bundle.writestr('manifest.json', json.dumps(manifest))
            bundle.writestr('input.pdb', raw)
        result = self.run_adapter({'batch_archive': 'bundle.zip'})
        self.assertEqual(result['source'], source_name)
        self.assertEqual(result['input_provenance']['identifier'], '1KLS')
        self.assertEqual(self.bundled_rows('sites.csv')[0]['source'], "'" + source_name)
        self.assertEqual(self.bundled_rows('structures.csv')[0]['source'], "'" + source_name)
        self.assertFalse(self.rows('00_Summary.csv')[0]['structure'].startswith(('=', '+', '-', '@')))
        self.assertEqual(self.rows('sites.csv')[0]['source'], "'" + source_name)

    def test_engine_and_export_failures_never_write_final_success_bundle(self):
        with self.assertRaisesRegex(RuntimeError, 'engine failed'):
            self.run_adapter(predictor=Mock(side_effect=RuntimeError('engine failed')))
        self.assertFalse((self.app / 'out/prediction.json').exists())
        with self.assertRaisesRegex(RuntimeError, 'export failed'):
            self.run_adapter(exporter=Mock(side_effect=RuntimeError('export failed')))
        self.assertFalse((self.app / 'out/prediction.json').exists())
        self.assertEqual(list((self.app / 'out').glob('*.zip')), [])

    def test_missing_or_reordered_engine_identities_are_rejected(self):
        for output in ([], [{'key': 'wrong'}], [{'key': 's0001_sample'}, {'key': 'extra'}]):
            with self.subTest(output=output), self.assertRaisesRegex(RuntimeError, 'identities'):
                self.run_adapter(predictor=Mock(return_value=output))
        self.assertFalse((self.app / 'out/prediction.json').exists())

    def test_missing_required_exports_and_paths_outside_output_are_rejected(self):
        for paths in ({'pdb': self.app / 'out/only.pdb'},
                      {'mmcif': self.app / 'outside.cif', 'pymol_session': self.app / 'outside.pse'}):
            with self.subTest(paths=paths), self.assertRaises((RuntimeError, ValueError)):
                self.run_adapter(exporter=Mock(return_value={'artifact_paths': paths}))
        self.assertFalse((self.app / 'out/prediction.json').exists())

    def test_invalid_scalar_settings_are_rejected_before_engine_execution(self):
        invalid = {'histidine_rotamers': ['maybe', ''],
                   'minimum_probability': ['nan', 'inf', '-0.1', '1.1'],
                   'cpu_workers': ['nan', 'inf', '0', '9', '2.5']}
        predictor = Mock()
        for key, values in invalid.items():
            for value in values:
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    self.run_adapter({key: value}, predictor)
        predictor.assert_not_called()

    def test_missing_or_mismatched_engine_pin_is_rejected(self):
        engine = self.app / 'engine'
        engine.mkdir()
        with patch.object(adapter, 'ENGINE_ROOT', engine):
            with self.assertRaisesRegex(RuntimeError, 'missing or mismatched'):
                adapter.load_predictor()
            (engine / 'SOURCE_REVISION').write_text('wrong revision\n')
            with self.assertRaisesRegex(RuntimeError, 'missing or mismatched'):
                adapter.load_predictor()

    def platform_outputs(self):
        output = self.app / 'out'
        output.mkdir()
        service_files = {'output.log': b'platform log before prediction\n', '.exit_code': b'running',
                         'settings.json': b'{"platform":"settings"}',
                         '.platform_state/state.json': b'{"platform":"service"}'}
        for name, contents in service_files.items():
            path = output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(contents)
        inodes = {output / name: (output / name).stat().st_ino for name in service_files}
        inodes[output / '.platform_state'] = (output / '.platform_state').stat().st_ino
        inodes[output] = output.stat().st_ino
        return output, service_files, inodes

    def test_details_stay_private_and_public_files_use_create_writes_without_platform_mutations(self):
        output, service_files, inodes = self.platform_outputs()
        export_directories, created = [], set()
        def recording_exporter(text, file_format, result, directory, stem):
            directory = Path(directory).resolve()
            self.assertFalse(directory.is_relative_to(output))
            self.assertTrue(directory.is_relative_to(Path('/tmp')))
            export_directories.append(directory)
            return self.exporter(text, file_format, result, directory, stem)

        def public_path(value, directory_fd=None):
            if isinstance(value, int):
                return None
            path = Path(os.fsdecode(value))
            if not path.is_absolute() and directory_fd is not None:
                path = Path(os.readlink(f'/proc/self/fd/{directory_fd}')) / path
            path = path.resolve()
            return path if path.is_relative_to(output) else None

        replace, rename, unlink, rmdir = os.replace, os.rename, os.unlink, os.rmdir
        def forbid_public_move(operation):
            def guarded(source, destination, *args, **kwargs):
                self.assertIsNone(public_path(source, kwargs.get('src_dir_fd')), 'Public output was renamed')
                self.assertIsNone(public_path(destination, kwargs.get('dst_dir_fd')), 'Public output was moved into')
                return operation(source, destination, *args, **kwargs)
            return guarded
        def forbid_public_delete(operation):
            def guarded(path, *args, **kwargs):
                self.assertIsNone(public_path(path, kwargs.get('dir_fd')), 'Public output was deleted')
                return operation(path, *args, **kwargs)
            return guarded
        def record_create(operation):
            def guarded(file, mode='r', *args, **kwargs):
                path = public_path(file)
                if path is not None and any(flag in mode for flag in ('w', 'a', 'x')):
                    self.assertEqual(path.parent, output, 'Private detail tree leaked into public output')
                    self.assertNotIn(path.name, service_files, 'Platform file was opened for writing')
                    self.assertFalse(path.exists(), 'Public publication overwrote an existing file')
                    created.add(path.name)
                return operation(file, mode, *args, **kwargs)
            return guarded
        with (output / 'output.log').open('ab') as log:
            with patch.object(os, 'replace', side_effect=forbid_public_move(replace)), \
                    patch.object(os, 'rename', side_effect=forbid_public_move(rename)), \
                    patch.object(os, 'unlink', side_effect=forbid_public_delete(unlink)), \
                    patch.object(os, 'rmdir', side_effect=forbid_public_delete(rmdir)), \
                    patch.object(builtins, 'open', side_effect=record_create(builtins.open)), \
                    patch.object(io, 'open', side_effect=record_create(io.open)):
                result = self.run_adapter(exporter=recording_exporter)
            log.write(b'platform log after prediction through original handle\n')
            log.flush()
        self.assertTrue(export_directories)
        self.assertTrue(all(not path.exists() for path in export_directories))
        expected = {'00_Summary.csv', 'sites.csv', 'CITATION.txt', 'prediction.json', result['result_bundle']}
        for item in result['structures']:
            expected.update((item['artifacts']['structure'], item['artifacts']['pymol_session']))
        self.assertEqual(created, expected)
        self.assertEqual({path.name for path in output.iterdir()},
                         expected | {Path(name).parts[0] for name in service_files})
        for path, inode in inodes.items():
            self.assertEqual(path.stat().st_ino, inode, str(path))
        for name, contents in service_files.items():
            if name != 'output.log':
                self.assertEqual((output / name).read_bytes(), contents)
        self.assertEqual((output / 'output.log').read_bytes(), service_files['output.log'] +
                         b'platform log after prediction through original handle\n')
        self.assertEqual(json.loads((output / 'prediction.json').read_text()), result)
        with zipfile.ZipFile(output / result['result_bundle']) as bundle:
            self.assertTrue(set(service_files).isdisjoint(bundle.namelist()))
            self.assertFalse(any(name.startswith('.platform_state/') for name in bundle.namelist()))

    def test_failed_private_export_leaves_platform_output_untouched_and_removes_private_work(self):
        output, service_files, inodes = self.platform_outputs()
        export_directories = []
        def failed_exporter(text, file_format, result, directory, stem):
            directory = Path(directory).resolve()
            export_directories.append(directory)
            self.assertFalse(directory.is_relative_to(output))
            directory.mkdir(parents=True, exist_ok=True)
            (directory / 'partial.pdb').write_text('incomplete private export')
            raise RuntimeError('private exporter failed')
        with (output / 'output.log').open('ab') as log:
            with self.assertRaisesRegex(RuntimeError, 'private exporter failed'):
                self.run_adapter(exporter=failed_exporter)
            log.write(b'platform log after failed private export\n')
            log.flush()
        self.assertTrue(export_directories)
        self.assertTrue(all(not directory.exists() for directory in export_directories))
        self.assertEqual({path.name for path in output.iterdir()}, {Path(name).parts[0] for name in service_files})
        for path, inode in inodes.items():
            self.assertEqual(path.stat().st_ino, inode, str(path))
        for name, contents in service_files.items():
            if name != 'output.log':
                self.assertEqual((output / name).read_bytes(), contents)
        self.assertEqual((output / 'output.log').read_bytes(),
                         service_files['output.log'] + b'platform log after failed private export\n')


if __name__ == '__main__':
    unittest.main()
