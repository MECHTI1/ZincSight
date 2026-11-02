"""Stdlib-only batch adapter contract tests; structures are checked in-image."""
import copy
import csv
import hashlib
import importlib.util
import json
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
        self.allocation = {
            'requested_cpu_workers': 4, 'cpu_workers': 2, 'declared_cpu_budget': 4,
            'declared_memory_gib': 12, 'detected_cpu_limit': 2, 'effective_memory_gib': 6,
            'memory_budget_gib_per_worker': 3, 'worker_limit_applied': True,
        }
        patcher = patch.object(adapter, 'worker_allocation', return_value=self.allocation)
        self.addCleanup(patcher.stop)
        patcher.start()

    def predictor(self, structures, **settings):
        self.calls.append((copy.deepcopy(structures), settings))
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
            path.write_text('export test artifact\n')
            paths[key] = path
        return {'artifact_paths': paths, 'warnings': [], 'predicted_marker_chain': 'Z'}

    def run_adapter(self, environment=None, predictor=None, exporter=None):
        return adapter.run(self.app, environment or {},
                           self.predictor if predictor is None else predictor,
                           self.exporter if exporter is None else exporter)

    def rows(self, filename):
        with (self.app / 'out' / filename).open(newline='') as stream:
            return list(csv.DictReader(stream))

    def test_single_input_settings_resources_and_download_contract(self):
        result = self.run_adapter({'structure': str(self.structure), 'histidine_rotamers': 'false',
                                   'minimum_probability': '0.4', 'max_sites': '12.0', 'cpu_workers': '8'})
        inputs, settings = self.calls[0]
        self.assertEqual(inputs[0]['structure_text'], PDB)
        self.assertEqual(settings, {'histidine_rotamers': False, 'minimum_probability': 0.4,
                                   'max_sites': 12, 'cpu_workers': 2,
                                   'engine_root': adapter.ENGINE_ROOT, 'timeout': 1200})
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
        self.assertEqual(json.loads((self.app / 'out/prediction.json').read_text()), result)
        summary = self.rows('structures.csv')[0]
        self.assertEqual(summary['structure_file'], 'structures/s0001_sample.pdb')
        self.assertEqual(summary['pymol_session'], 'structures/s0001_sample.pse')
        self.assertEqual(summary['returned_site_count'], '1')
        row = self.rows('sites.csv')[0]
        self.assertEqual(json.loads(row['coordinating_residues']), SITE['coordinating_residues'])
        self.assertEqual(float(row['z_angstrom']), 3.3)
        for path in result['artifacts'].values():
            self.assertTrue((self.app / 'out' / path).is_file())
        per_input = json.loads((self.app / 'out' / result['prediction_json']).read_text())
        self.assertEqual(per_input['sites'], [SITE])
        with zipfile.ZipFile(self.app / 'out/zincsight_results.zip') as archive:
            self.assertIsNone(archive.testzip())
            self.assertIn(summary['pymol_session'], archive.namelist())
            self.assertIn(summary['structure_file'], archive.namelist())
            self.assertIn('CITATION.bib', archive.namelist())
            self.assertNotIn('zincsight_results.zip', archive.namelist())

    def test_citation_is_zincsight_paper_in_all_machine_and_download_formats(self):
        result = self.run_adapter()
        self.assertEqual(result['reference']['doi'], '10.1002/pro.70350')
        self.assertIn('Mechtinger G', result['citation'])
        self.assertIn('Protein Science. 2025;34(11):e70350', result['citation'])
        self.assertIn('10.1002/pro.70350', (self.app / 'out/CITATION.txt').read_text())
        self.assertIn('doi = {10.1002/pro.70350}', (self.app / 'out/CITATION.bib').read_text())
        self.assertEqual(json.loads((self.app / 'out/citation.json').read_text()), result['reference'])

    def test_batch_keeps_duplicate_basenames_distinct_and_zero_site_inputs_visible(self):
        second = self.app / 'inputs/nested/sample.pdb'
        second.parent.mkdir()
        second.write_text(PDB.replace('0.000', '1.000'))
        def mixed(inputs, **settings):
            results = self.predictor(inputs, **settings)
            results[1].update(sites=[], site_count=0)
            return results
        result = self.run_adapter({'structures': json.dumps(['sample.pdb', 'nested/sample.pdb'])}, mixed)
        self.assertEqual(result['structure_count'], 2)
        self.assertNotIn('sites', result)
        self.assertEqual(result['returned_site_count'], 1)
        self.assertEqual([item['key'] for item in result['structures']], ['s0001_sample', 's0002_sample'])
        summaries = self.rows('structures.csv')
        self.assertEqual(len(summaries), 2)
        self.assertEqual(summaries[1]['returned_site_count'], '0')
        self.assertNotEqual(summaries[0]['structure_file'], summaries[1]['structure_file'])
        self.assertTrue((self.app / 'out' / summaries[1]['pymol_session']).is_file())
        self.assertEqual(len(self.rows('sites.csv')), 1)

    def test_zero_sites_still_exports_structure_and_session_with_header_only_sites(self):
        def empty(inputs, **settings):
            results = self.predictor(inputs, **settings)
            results[0].update(sites=[], site_count=0)
            return results
        result = self.run_adapter(predictor=empty)
        self.assertEqual(result['returned_site_count'], 0)
        self.assertFalse(result['truncated'])
        self.assertEqual(self.rows('sites.csv'), [])
        self.assertEqual(len(self.rows('structures.csv')), 1)
        self.assertTrue((self.app / 'out' / result['artifacts']['pymol_session']).is_file())

    def test_site_limits_and_truncation_are_recorded_per_structure(self):
        (self.app / 'inputs/second.cif').write_text('data_x\n_atom_site.id 1\n')
        def truncated(inputs, **settings):
            results = self.predictor(inputs, **settings)
            for result, total in zip(results, (5, 8)):
                result.update(site_count=total, truncated=True)
            return results
        result = self.run_adapter({'max_sites': '1'}, truncated)
        self.assertEqual((result['site_count'], result['returned_site_count']), (13, 2))
        self.assertTrue(result['truncated'])
        self.assertEqual(result['max_sites_per_structure'], 1)
        self.assertEqual([item['max_sites'] for item in result['structures']], [1, 1])

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
        self.assertEqual(self.rows('sites.csv')[0]['source'], "'" + source_name)
        self.assertEqual(self.rows('structures.csv')[0]['source'], "'" + source_name)

    def test_engine_and_export_failures_never_write_final_success_bundle(self):
        with self.assertRaisesRegex(RuntimeError, 'engine failed'):
            self.run_adapter(predictor=Mock(side_effect=RuntimeError('engine failed')))
        self.assertFalse((self.app / 'out/prediction.json').exists())
        with self.assertRaisesRegex(RuntimeError, 'export failed'):
            self.run_adapter(exporter=Mock(side_effect=RuntimeError('export failed')))
        self.assertFalse((self.app / 'out/prediction.json').exists())
        self.assertFalse((self.app / 'out/zincsight_results.zip').exists())

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
                   'max_sites': ['nan', 'inf', '0', '101', '2.5'],
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


if __name__ == '__main__':
    unittest.main()
