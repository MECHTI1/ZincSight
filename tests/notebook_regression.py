"""Reproduce notebook examples on original/native and current/embedded code.

Historical differences are reported, not silently accepted as new goldens.
Current native/embedded and single/batch differences fail this check.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_MINIMUM_PROBABILITY = 0.0001
SOURCES = {
    '1KLS': 'https://files.rcsb.org/download/1KLS.cif',
    '2A0S-assembly1': 'https://files.rcsb.org/download/2A0S-assembly1.cif',
    '8QEP': 'https://files.rcsb.org/download/8QEP.cif',
    'AF-P0A6G5-F1-model_v6': 'https://alphafold.ebi.ac.uk/files/AF-P0A6G5-F1-model_v6.cif',
    'AF-A0A068N621-F1-model_v6': 'https://alphafold.ebi.ac.uk/files/AF-A0A068N621-F1-model_v6.cif',
    'MGYP002718891411': 'https://api.esmatlas.com/fetchPredictedStructure/MGYP002718891411.pdb',
}


def download_inputs(directory):
    from Bio.PDB import PDBParser, MMCIFIO
    directory.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, url in SOURCES.items():
        extension = url.rsplit('.', 1)[-1]
        path = directory / f'{name}.{extension}'
        if not path.exists():
            with urllib.request.urlopen(url, timeout=40) as response:
                content = response.read(20_000_001)
            if len(content) > 20_000_000:
                raise ValueError(f'Unexpectedly large input: {name}')
            path.write_bytes(content)
        manifest.append({'id': name, 'url': url, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        if extension == 'pdb':
            structure = PDBParser(QUIET=True).get_structure(name, str(path))
            writer = MMCIFIO()
            writer.set_structure(structure)
            writer.save(str(path.with_suffix('.cif')))
    (directory / 'sources.json').write_text(json.dumps(manifest, indent=2))


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def notebook_rows(data):
    return [row for row in data if float(row['prob']) >= NOTEBOOK_MINIMUM_PROBABILITY]


def site_key(row):
    return (row['structure_id'], tuple(sorted(row[f'chain_resi_{i}'] for i in range(1, 5)
                                           if row[f'chain_resi_{i}'])))


def compare(expected, actual):
    e, a = {site_key(r): r for r in expected}, {site_key(r): r for r in actual}
    if len(e) != len(expected) or len(a) != len(actual):
        raise AssertionError('Duplicate coordinating-residue sets need explicit matching.')
    common = e.keys() & a.keys()
    return {
        'expected': len(e), 'actual': len(a), 'shared_residue_sets': len(common),
        'missing': sorted(e.keys() - a.keys()), 'extra': sorted(a.keys() - e.keys()),
        'max_probability_delta': max((abs(float(e[k]['prob']) - float(a[k]['prob'])) for k in common), default=0),
        'max_position_delta': max((math.dist(json.loads(e[k]['predicted_ion_pos']),
                                            json.loads(a[k]['predicted_ion_pos'])) for k in common), default=0),
        'changed_templates': sum(e[k]['template_id'] != a[k]['template_id'] for k in common),
    }


def predict(source, inputs, names, backend, output, label, cores=2):
    """Run the full engine in an isolated process; retain its CSV and log."""
    with tempfile.TemporaryDirectory(prefix='zincsight-regression-') as directory:
        work = Path(directory)
        shutil.copytree(source / 'src', work / 'src', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(source / 'main_execute.py', work)
        shutil.copy2(ROOT / 'integrations/mcp/embedded_database.py', work)
        (work / 'query').mkdir()
        for name in names:
            shutil.copy2(inputs / f'{name}.cif', work / 'query')
        if backend == 'pglite':
            code = (
                'from pathlib import Path\nfrom embedded_database import embedded_database\n'
                f'with embedded_database(Path.cwd(), Path({str(ROOT / "integrations/mcp/pglite/server.mjs")!r})):\n'
                '    from src.primary_script import main\n'
                f'    main(sorted(str(p) for p in Path("query").glob("*.cif")), True, "results", {cores}, False)\n'
            )
        else:
            code = ('from main_execute import execute_zincsight\n'
                    f'execute_zincsight(True, "", "query", "results", {cores}, False)\n')
        (work / 'run.py').write_text(code)
        env = {**os.environ, 'PYTHON_DOTENV_DISABLED': '1', 'DEBUGGING_TEMP': 'false',
               'KEEP_TEMP_TABLES_TEMP': 'false', 'FILTER_PROB_AFTER_COMPRESSION': 'false',
               'EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD': '0.0', 'OMP_NUM_THREADS': '1',
               'OPENBLAS_NUM_THREADS': '1', 'MPLBACKEND': 'Agg'}
        print(f'RUN {label} ({backend}, {len(names)} structures, {cores} cores)', flush=True)
        with (output / f'{label}.log').open('w') as log:
            process = subprocess.Popen([sys.executable, 'run.py'], cwd=work, env=env,
                                       stdout=log, stderr=log, start_new_session=True)
            try:
                code = process.wait(timeout=420)
            finally:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
        if code:
            raise RuntimeError(f'{label} failed; see its saved log')
        result = output / f'{label}.csv'
        shutil.copy2(work / 'results/table/table_predicted_zn_binding_sites.csv', result)
        return rows(result)


def predict_through_wrapper(name, output):
    sys.path.insert(0, str(ROOT))
    from integrations.mcp.runner import run_prediction
    # Exercise the UniProt shorthand as well as explicit AlphaFold model IDs.
    query = 'P0A6G5' if name == 'AF-P0A6G5-F1-model_v6' else name
    result = run_prediction(query_id=query, histidine_rotamers=True,
                            max_sites=100, database_backend='pglite')
    assert result['minimum_probability'] == NOTEBOOK_MINIMUM_PROBABILITY
    assert not result['truncated']
    (output / f'wrapper-{name}.json').write_text(json.dumps(result, indent=2))
    actual = []
    for site in result['sites']:
        row = {'structure_id': name, 'prob': str(site['probability']),
               'predicted_ion_pos': json.dumps(site['coordinates_angstrom']),
               'template_id': site['template_id']}
        for i in range(1, 5):
            residues = site['coordinating_residues']
            row[f'chain_resi_{i}'] = residues[i - 1] if i <= len(residues) else ''
        actual.append(row)
    return actual


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-root', type=Path)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--embedded-only', action='store_true')
    parser.add_argument('--expected', type=Path, help='Optional local reference CSV; never uploaded by this script.')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = args.input_dir.resolve()
    download_inputs(inputs)
    results = {}
    if args.original_root:
        results['original-native'] = predict(args.original_root.resolve(), inputs, SOURCES, 'native', output, 'original-native')
    if not args.embedded_only:
        results['current-native'] = predict(ROOT, inputs, SOURCES, 'native', output, 'current-native')
    results['current-embedded'] = predict(ROOT, inputs, SOURCES, 'pglite', output, 'current-embedded')
    expected = rows(args.expected) if args.expected else results.get('original-native')
    report = {'reference': 'local_csv' if args.expected else 'original_native_fresh_run',
              'reference_minimum_probability': NOTEBOOK_MINIMUM_PROBABILITY,
              'raw_engine_counts': {label: len(actual) for label, actual in results.items()},
              'notebook_counts': {label: len(notebook_rows(actual)) for label, actual in results.items()},
              'reference_comparisons': {}, 'consistency': {}}
    if expected is not None:
        for label, actual in results.items():
            report['reference_comparisons'][label] = {
                name: compare([r for r in notebook_rows(expected) if r['structure_id'] == name],
                              [r for r in notebook_rows(actual) if r['structure_id'] == name]) for name in SOURCES
            }
            print(label, json.dumps(report['reference_comparisons'][label]), flush=True)
    if 'current-native' in results:
        report['consistency']['native_vs_embedded'] = compare(results['current-native'], results['current-embedded'])
    for name in SOURCES:
        print(f'RUN wrapper-{name}', flush=True)
        single = predict_through_wrapper(name, output)
        batch = [r for r in notebook_rows(results['current-embedded']) if r['structure_id'] == name]
        report['consistency'][f'batch_vs_single-{name}'] = compare(batch, single)
    (output / 'report.json').write_text(json.dumps(report, indent=2))
    print('CONSISTENCY', json.dumps(report['consistency']), flush=True)
    for label, result in report['consistency'].items():
        assert not result['missing'] and not result['extra'], (label, result)
        assert result['max_probability_delta'] <= 1e-7, (label, result)
        assert result['max_position_delta'] <= 1e-6, (label, result)
        assert result['changed_templates'] == 0, (label, result)


if __name__ == '__main__':
    main()
