"""Tamarind batch adapter with native structure previews and PyMOL downloads."""
from __future__ import annotations
import csv
import json
import math
import os
from pathlib import Path
import sys
import time
import zipfile

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from input_bundle import collect_structures
from resources import worker_allocation

ENGINE_ROOT = Path('/opt/zincsight')
ENGINE_REVISION = '7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03'
INTEGRATION_VERSION = '0.2.0'
CITATION = (
    'Mechtinger G, Axel G, Kolodny R, Ben-Tal N. Interpretable prediction of zinc ion '
    'location in proteins with ZincSight. Protein Science. 2025;34(11):e70350. '
    'https://doi.org/10.1002/pro.70350'
)
BIBTEX = '''@article{Mechtinger2025ZincSight,
  author = {Mechtinger, G. and Axel, G. and Kolodny, R. and Ben-Tal, N.},
  title = {Interpretable prediction of zinc ion location in proteins with ZincSight},
  journal = {Protein Science},
  year = {2025},
  volume = {34},
  number = {11},
  pages = {e70350},
  doi = {10.1002/pro.70350},
  url = {https://doi.org/10.1002/pro.70350}
}
'''
VALIDATION_STATUS = (
    'Development preview with reproducibility fixes; historical outputs can differ. '
    'Published accuracy benchmarks have not been re-established for this revision.'
)
SITE_FIELDS = ['structure_id', 'source', 'site_rank', 'probability', 'x_angstrom',
               'y_angstrom', 'z_angstrom', 'coordinating_residues', 'residue_combination',
               'template_id', 'template_rmsd_angstrom']
STRUCTURE_FIELDS = ['structure_id', 'source', 'site_count', 'returned_site_count',
                    'highest_probability', 'truncated', 'protein_residue_count',
                    'structure_file', 'mmcif_file', 'pymol_session', 'prediction_json',
                    'site_annotations', 'warnings']


def read_parameters(environment):
    flag = str(environment.get('histidine_rotamers', 'true')).strip().lower()
    if flag not in ('true', 'false', '1', '0'):
        raise ValueError('histidine_rotamers must be true or false.')
    try:
        probability = float(environment.get('minimum_probability', '0.0001'))
        count = float(environment.get('max_sites', '100'))
        workers = float(environment.get('cpu_workers', '4'))
    except (ValueError, TypeError) as exc:
        raise ValueError('Probability, site count and CPU workers must be numbers.') from exc
    if not math.isfinite(probability) or not 0 <= probability <= 1:
        raise ValueError('minimum_probability must be between 0 and 1.')
    if not math.isfinite(count) or not count.is_integer() or not 1 <= count <= 100:
        raise ValueError('max_sites must be a whole number between 1 and 100 per structure.')
    if not math.isfinite(workers) or not workers.is_integer() or not 1 <= workers <= 8:
        raise ValueError('cpu_workers must be a whole number between 1 and 8.')
    return {'histidine_rotamers': flag in ('true', '1'), 'minimum_probability': probability,
            'max_sites': int(count), 'cpu_workers': int(workers)}


def load_predictor():
    revision = ENGINE_ROOT / 'SOURCE_REVISION'
    if not revision.is_file() or revision.read_text().strip() != ENGINE_REVISION:
        raise RuntimeError('The pinned ZincSight engine is missing or mismatched; rebuild the image.')
    from batch_engine import run_predictions
    return run_predictions


def _write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def _write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _csv_label(value):
    return "'" + value if value.startswith(('=', '+', '-', '@')) else value


def run(app: Path, environment, predictor=None, exporter=None) -> dict:
    started = time.perf_counter()
    app = Path(app)
    parameters = read_parameters(environment)
    resources = worker_allocation(parameters.pop('cpu_workers'), environment)
    structures = collect_structures(app, environment)
    predictor = predictor or load_predictor()
    if exporter is None:
        from structure_outputs import export_structure_outputs
        exporter = export_structure_outputs
    prediction_started = time.perf_counter()
    results = predictor(structures, **parameters, cpu_workers=resources['cpu_workers'],
                        engine_root=ENGINE_ROOT, timeout=1200)
    prediction_seconds = time.perf_counter() - prediction_started
    if (not isinstance(results, list) or len(results) != len(structures)
            or not all(isinstance(result, dict) for result in results)
            or [r.get('key') for r in results] != [s['key'] for s in structures]):
        raise RuntimeError('The engine returned inconsistent structure identities.')
    output = app / 'out'
    output.mkdir(parents=True, exist_ok=True)
    reference = {'citation': CITATION, 'doi': '10.1002/pro.70350',
                 'url': 'https://doi.org/10.1002/pro.70350', 'bibtex': BIBTEX}
    (output / 'CITATION.txt').write_text('Please cite the ZincSight method:\n\n' + CITATION + '\n', encoding='utf-8')
    (output / 'CITATION.bib').write_text(BIBTEX, encoding='utf-8')
    _write_json(output / 'citation.json', reference)
    all_site_rows, summary_rows = [], []
    for source, result in zip(structures, results):
        result.update(
            source=source['source'], input_sha256=source['input_sha256'],
            input_provenance=source['provenance'], input_bytes=source['input_bytes'],
            software_revision=ENGINE_REVISION, integration='tamarind',
            integration_version=INTEGRATION_VERSION, database_transport='unix_socket',
            validation_status=VALIDATION_STATUS, citation=CITATION, reference=reference,
            returned_site_count=len(result['sites']), max_sites=parameters['max_sites'],
            resources=resources,
        )
        artifacts = exporter(source['structure_text'], source['file_format'], result,
                             output / 'structures', stem=source['key'])
        paths = {name: Path(path).resolve().relative_to(output.resolve()).as_posix()
                 for name, path in artifacts['artifact_paths'].items()}
        if 'mmcif' not in paths or 'pymol_session' not in paths:
            raise RuntimeError('The exporter did not produce the required structure and PyMOL session.')
        if any(not (output / path).is_file() for path in paths.values()):
            raise RuntimeError('The exporter returned a missing structure artifact.')
        result['artifacts'] = paths
        result['visualization'] = {key: value for key, value in artifacts.items() if key != 'artifact_paths'}
        prediction_path = f"predictions/{source['key']}.json"
        result['prediction_json'] = prediction_path
        _write_json(output / prediction_path, result)
        rows = []
        for rank, site in enumerate(result['sites'], 1):
            x, y, z = site['coordinates_angstrom']
            rows.append({'structure_id': source['key'], 'source': _csv_label(source['source']),
                         'site_rank': rank, 'probability': site['probability'],
                         'x_angstrom': x, 'y_angstrom': y, 'z_angstrom': z,
                         'coordinating_residues': json.dumps(site['coordinating_residues']),
                         'residue_combination': site['residue_combination'],
                         'template_id': site['template_id'],
                         'template_rmsd_angstrom': site['template_rmsd_angstrom']})
        all_site_rows.extend(rows)
        _write_csv(output / 'sites' / f"{source['key']}.csv", SITE_FIELDS, rows)
        summary_rows.append({
            'structure_id': source['key'], 'source': _csv_label(source['source']),
            'site_count': result['site_count'], 'returned_site_count': len(result['sites']),
            'highest_probability': result['sites'][0]['probability'] if result['sites'] else '',
            'truncated': result['truncated'], 'protein_residue_count': result['protein_residue_count'],
            'structure_file': paths.get('pdb') or paths['mmcif'], 'mmcif_file': paths['mmcif'],
            'pymol_session': paths['pymol_session'], 'prediction_json': prediction_path,
            'site_annotations': paths['site_annotations'],
            'warnings': '; '.join(artifacts.get('warnings', [])),
        })
    _write_csv(output / 'structures.csv', STRUCTURE_FIELDS, summary_rows)
    _write_csv(output / 'sites.csv', SITE_FIELDS, all_site_rows)
    result = {
        'method': 'ZincSight', 'structure_count': len(results), 'structures': results,
        'site_count': sum(item['site_count'] for item in results),
        'returned_site_count': sum(item['returned_site_count'] for item in results),
        'truncated': any(item['truncated'] for item in results), 'max_sites_per_structure': parameters['max_sites'],
        'histidine_rotamers': parameters['histidine_rotamers'],
        'minimum_probability': parameters['minimum_probability'],
        'citation': CITATION, 'reference': reference, 'software_revision': ENGINE_REVISION,
        'integration': 'tamarind', 'integration_version': INTEGRATION_VERSION,
        'database_transport': 'unix_socket', 'validation_status': VALIDATION_STATUS,
        'resources': resources,
        'performance': {'prediction_seconds': prediction_seconds,
                        'total_seconds': time.perf_counter() - started,
                        'database_instances_per_batch': 1},
    }
    if len(results) == 1:
        result.update({key: value for key, value in results[0].items()
                       if key not in ('resources', 'reference')})
    _write_json(output / 'prediction.json', result)
    (output / 'RESULTS_README.txt').write_text(
        'ZincSight results\n\n'
        'structures.csv: one row and viewable structure per input, including zero-site results.\n'
        'sites.csv: sites with residues, coordinates, probabilities and templates.\n'
        'structures/: combined first-model PDB/mmCIF, PyMOL .pse, and marker annotations.\n'
        'predictions/: per-structure prediction JSON.\n'
        'CITATION.txt / CITATION.bib: the ZincSight paper citation.\n\n'
        'Predicted Zn markers locate potential binding sites; Zn labels do not assign bound metal identity.\n'
        'Original atoms and experimental ions are retained. Marker annotations distinguish predictions.\n'
        'PDB coordinates have 0.001 Angstrom precision; mmCIF and JSON retain higher precision.\n'
        'The first model is analyzed and exported; later models are not analyzed.\n'
        'The site limit applies independently to each input structure.\n\n' + CITATION + '\n', encoding='utf-8')
    with zipfile.ZipFile(output / 'zincsight_results.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(output.rglob('*')):
            if path.is_file() and path.name != 'zincsight_results.zip':
                archive.write(path, path.relative_to(output).as_posix())
    print(f"Processed {len(results)} structures with {resources['cpu_workers']} CPU workers; "
          f"returned {result['returned_site_count']} sites. PDB/mmCIF and PyMOL sessions are ready.")
    return result


if __name__ == '__main__':
    try:
        run(Path('/app'), os.environ)
    except Exception as exc:
        print(f'ZincSight failed: {exc}', file=sys.stderr)
        sys.exit(1)
