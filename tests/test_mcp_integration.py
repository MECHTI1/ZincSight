"""Opt-in real predictions. Requires a non-root Linux host and RCSB access."""
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from integrations.mcp.runner import ROOT, get_structure, parse_sites, run_prediction

pytestmark = pytest.mark.skipif(
    os.environ.get('ZINCSIGHT_RUN_INTEGRATION') != '1',
    reason='Set ZINCSIGHT_RUN_INTEGRATION=1 on a non-root Linux host.',
)


@pytest.fixture(scope='module')
def structure_text():
    return get_structure('1A1F', None, 'pdb')


def cif_text(pdb_text, scale=1):
    from Bio.PDB import PDBParser, MMCIFIO
    model = next(PDBParser(QUIET=True).get_structure('query', io.StringIO(pdb_text)).get_models())
    for atom in model.get_atoms():
        atom.set_coord(atom.coord * scale)
    writer = MMCIFIO()
    writer.set_structure(model)
    stream = io.StringIO()
    writer.save(stream)
    return stream.getvalue()


def direct_prediction(tmp_path, text, rotamers):
    """Call the original entry point directly, without the MCP runner/worker."""
    shutil.copytree(ROOT / 'src', tmp_path / 'src', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(ROOT / 'main_execute.py', tmp_path / 'main_execute.py')
    (tmp_path / 'query').mkdir()
    (tmp_path / 'query/query.cif').write_text(cif_text(text))
    (tmp_path / 'baseline.py').write_text(
        'from main_execute import execute_zincsight\n'
        'from pathlib import Path\n'
        'import json\n'
        'if __name__ == "__main__":\n'
        f'    result = execute_zincsight({rotamers!r}, "", "query", "results", 2, False)\n'
        '    Path("status.json").write_text(json.dumps({"has_sites": result is not False}))\n'
    )
    env = {**os.environ, 'PYTHON_DOTENV_DISABLED': '1', 'FILTER_PROB_AFTER_COMPRESSION': 'false',
           'DEBUGGING_TEMP': 'false', 'KEEP_TEMP_TABLES_TEMP': 'false',
           'EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD': '0.0'}
    result = subprocess.run([sys.executable, 'baseline.py'], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout[-1500:] + result.stderr[-4500:]
    status = json.loads((tmp_path / 'status.json').read_text())
    return parse_sites(tmp_path / 'results/table/table_predicted_zn_binding_sites.csv', 0) if status['has_sites'] else []


def comparable(sites):
    return sorted(sites, key=lambda s: tuple(sorted(s['coordinating_residues'])))


@pytest.mark.parametrize('rotamers', [False, True])
def test_matches_direct_predictor(tmp_path, structure_text, rotamers):
    expected = direct_prediction(tmp_path, structure_text, rotamers)
    assert expected, 'Positive control 1A1F should have predicted zinc sites.'
    result = run_prediction(structure_text=structure_text, histidine_rotamers=rotamers, max_sites=100)
    assert not result['truncated']
    assert comparable(result['sites']) == comparable(expected)


def test_mmcif_matches_pdb(structure_text):
    pdb = run_prediction(structure_text=structure_text, max_sites=100)
    cif = run_prediction(structure_text=cif_text(structure_text), file_format='cif', max_sites=100)
    assert comparable(pdb['sites']) == comparable(cif['sites'])


def test_no_sites_control(structure_text):
    # Artificial geometry control, not a biological nonbinder: scaling all atom
    # coordinates separates potential ligands beyond normal coordination distance.
    result = run_prediction(structure_text=cif_text(structure_text, scale=10), file_format='cif')
    assert result['sites'] == []
    assert result['site_count'] == 0
