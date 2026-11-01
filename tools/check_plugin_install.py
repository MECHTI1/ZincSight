"""Extract a built plugin and run its actual launcher, with a clean cache by default."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--cache-dir', type=Path)
    parser.add_argument('--structure', type=Path, help='Optional public 1KLS mmCIF for offline prediction testing')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='zincsight-plugin-check-') as temporary:
        root = Path(temporary)
        with zipfile.ZipFile(args.archive) as archive:
            for entry in archive.infolist():
                assert (root / entry.filename).resolve().is_relative_to(root.resolve())
            archive.extractall(root)
        script = root / 'skills/zincsight-preview/scripts/predict.py'
        environment = os.environ.copy()
        environment.update(ZINCSIGHT_CACHE_DIR=str(args.cache_dir.resolve() if args.cache_dir else root / 'cache'),
                           OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MPLBACKEND='Agg')
        command = [sys.executable, str(script)]
        if args.structure:
            command += ['--structure', str(args.structure.resolve()), '--format', 'cif']
        else:
            command += ['--query-id', '1KLS']
        completed = subprocess.run(command, env=environment, text=True, capture_output=True, timeout=1000)
        if completed.returncode:
            print(completed.stderr[-6000:], file=sys.stderr)
            raise RuntimeError('Extracted plugin launcher failed.')
        result = json.loads(completed.stdout)
        assert result['database_backend'] == 'pglite'
        assert result['histidine_rotamers'] is True
        assert result['minimum_probability'] == 0.0001
        assert result['site_count'] == len(result['sites']) == 1, '1KLS smoke control must contain one site.'
        assert result['software_revision'] == 'ac1fe1e96e8b2ffe4609ab33f9beff579c76ca07'
        assert 'Experimental' in result['validation_status']
        site = result['sites'][0]
        assert 0 <= site['probability'] <= 1
        assert len(site['coordinates_angstrom']) == 3
        assert all(math.isfinite(value) for value in site['coordinates_angstrom'])
        assert set(site['coordinating_residues']) == {'A_5', 'A_8', 'A_21', 'A_26'}
        bad = subprocess.run([sys.executable, str(script), '--structure', str(root / 'missing.cif')],
                             env=environment, text=True, capture_output=True, timeout=10)
        assert bad.returncode != 0 and not bad.stdout.strip(), 'An invalid input must not produce predictions.'
        report = {'package': args.archive.name, 'clean_cache': args.cache_dir is None,
                  'positive_prediction': result, 'invalid_input_rejected': True}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps({'passed': True, 'clean_cache': args.cache_dir is None,
                          'site_count': result['site_count'], 'invalid_input_rejected': True}))


if __name__ == '__main__':
    main()
