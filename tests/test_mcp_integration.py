"""Opt-in real predictions; native comparisons need a non-root Linux host."""
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from integrations.mcp.runner import ROOT, get_structure, parse_sites, run_prediction

requires_native = pytest.mark.skipif(
    os.environ.get('ZINCSIGHT_RUN_INTEGRATION') != '1',
    reason='Set ZINCSIGHT_RUN_INTEGRATION=1 on a non-root Linux host.',
)
requires_pglite = pytest.mark.skipif(
    os.environ.get('ZINCSIGHT_RUN_PGLITE') != '1',
    reason='Set ZINCSIGHT_RUN_PGLITE=1 after installing the optional npm dependencies.',
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
    """Call this checkout's standard entry point without the MCP runner/worker."""
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


@requires_native
@pytest.mark.parametrize('backend', ['native', 'pglite'])
@pytest.mark.parametrize('rotamers', [False, True])
def test_matches_direct_predictor(tmp_path, structure_text, rotamers, backend):
    expected = direct_prediction(tmp_path, structure_text, rotamers)
    assert expected, 'Positive control 1A1F should have predicted zinc sites.'
    result = run_prediction(structure_text=structure_text, histidine_rotamers=rotamers,
                            max_sites=100, database_backend=backend)
    assert not result['truncated']
    assert comparable(result['sites']) == comparable(expected)


@requires_native
def test_mmcif_matches_pdb(structure_text):
    pdb = run_prediction(structure_text=structure_text, max_sites=100, database_backend='native')
    cif = run_prediction(structure_text=cif_text(structure_text), file_format='cif', max_sites=100, database_backend='native')
    assert comparable(pdb['sites']) == comparable(cif['sites'])


@requires_native
def test_no_sites_control(structure_text):
    # Artificial geometry control, not a biological nonbinder: scaling all atom
    # coordinates separates potential ligands beyond normal coordination distance.
    result = run_prediction(structure_text=cif_text(structure_text, scale=10), file_format='cif', database_backend='native')
    assert result['sites'] == []
    assert result['site_count'] == 0


@requires_pglite
@pytest.mark.parametrize('rotamers', [False, True])
def test_embedded_prediction_and_formats(structure_text, rotamers):
    pdb = run_prediction(structure_text=structure_text, histidine_rotamers=rotamers, database_backend='pglite')
    cif = run_prediction(structure_text=cif_text(structure_text), file_format='cif',
                         histidine_rotamers=rotamers, database_backend='pglite')
    assert pdb['site_count'] > 0
    assert sum(site['probability'] > .9 for site in pdb['sites']) == 3
    assert comparable(pdb['sites']) == comparable(cif['sites'])


@requires_pglite
def test_embedded_jobs_are_isolated(structure_text):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as pool:
        positive = pool.submit(run_prediction, structure_text=structure_text, database_backend='pglite')
        negative = pool.submit(run_prediction, structure_text=cif_text(structure_text, scale=10),
                               file_format='cif', database_backend='pglite')
        assert positive.result()['site_count'] > 0
        assert negative.result()['sites'] == []


@requires_pglite
def test_embedded_prediction_through_mcp(structure_text):
    import asyncio
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    async def exercise():
        parameters = StdioServerParameters(command=sys.executable,
            args=['-m', 'integrations.mcp.server'], cwd=str(ROOT),
            env={**os.environ, 'ZINCSIGHT_DATABASE_BACKEND': 'pglite'})
        async with stdio_client(parameters) as (reader, writer):
            async with ClientSession(reader, writer) as session:
                await session.initialize()
                result = await session.call_tool('predict_zinc_binding_sites', {'structure_text': structure_text})
                assert not result.isError
                payload = result.structuredContent or json.loads(result.content[0].text)
                assert payload['site_count'] > 0
                assert payload['database_backend'] == 'pglite'
    asyncio.run(exercise())


@requires_native
@pytest.mark.parametrize('pdb_id', ['2CBA', '1UBQ'])
@pytest.mark.parametrize('rotamers', [False, True])
def test_embedded_matches_native_additional_structures(pdb_id, rotamers):
    text = get_structure(pdb_id, None, 'pdb')
    native = run_prediction(structure_text=text, histidine_rotamers=rotamers, database_backend='native')
    embedded = run_prediction(structure_text=text, histidine_rotamers=rotamers, database_backend='pglite')
    assert not native['truncated'] and not embedded['truncated']
    assert comparable(native['sites']) == comparable(embedded['sites'])
