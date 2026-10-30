import csv
import json
from pathlib import Path
import subprocess
from unittest.mock import patch

import pytest
from integrations.mcp import runner


@pytest.mark.parametrize('kwargs', [
    {}, {'pdb_id': '../etc/passwd'}, {'pdb_id': '1A1F;echo'},
    {'pdb_id': '1A1F', 'structure_text': 'ATOM  '},
    {'structure_text': 'not a structure'},
    {'structure_text': 'ATOM  ' + 'x' * runner.MAX_BYTES},
    {'pdb_id': '1A1F', 'minimum_probability': float('nan')},
    {'pdb_id': '1A1F', 'max_sites': 0},
    {'query_id': 'https://example.com/structure.cif'},
    {'query_id': '../1KLS'}, {'query_id': '1KLS?token=abc'},
    {'query_id': '1KLS', 'structure_text': 'ATOM  '},
])
def test_invalid_input(kwargs):
    with pytest.raises(ValueError):
        runner.run_prediction(**kwargs)


def test_parse_real_repository_results():
    sites = runner.parse_sites(runner.ROOT / 'Data/RMSD_plDDT_dataset/table_predicted_zn_binding_sites.csv', .99)
    assert sites
    assert all(.99 <= s['probability'] <= 1 for s in sites)
    assert sites[0]['coordinating_residues'] == ['A_412', 'A_453', 'A_415', 'A_450']
    assert len(sites[0]['coordinates_angstrom']) == 3


def test_cif_first_model_preserves_coordinate_precision(tmp_path):
    from Bio.PDB import MMCIFParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from integrations.mcp.worker import write_first_model
    original = tmp_path / 'original.cif'
    original.write_text('''data_precision
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 12.123456 2.234567 3.345678 1.0 90.0 1 A 1
ATOM 2 C CA . ALA A 1 1 ? 99.987654 8.876543 7.765432 1.0 90.0 1 A 2
#
''')
    model = next(MMCIFParser(QUIET=True).get_structure('query', original).get_models())
    target = tmp_path / 'first.cif'
    write_first_model(model, original, 'cif', target)
    data = MMCIF2Dict(str(target))
    assert data['_atom_site.Cartn_x'] == ['12.123456']
    assert data['_atom_site.Cartn_y'] == ['2.234567']
    assert data['_atom_site.Cartn_z'] == ['3.345678']
    assert data['_atom_site.pdbx_PDB_model_num'] == ['1']


@pytest.mark.parametrize('accession,host,path,file_format', [
    ('1KLS', 'files.rcsb.org', '/download/1KLS.cif', 'cif'),
    ('2A0S-assembly1', 'files.rcsb.org', '/download/2A0S-assembly1.cif', 'cif'),
    ('8QEP', 'files.rcsb.org', '/download/8QEP.cif', 'cif'),
    ('P0A6G5', 'alphafold.ebi.ac.uk', '/files/AF-P0A6G5-F1-model_v6.cif', 'cif'),
    ('AF-A0A068N621-F1-model_v6', 'alphafold.ebi.ac.uk', '/files/AF-A0A068N621-F1-model_v6.cif', 'cif'),
    ('MGYP002718891411', 'api.esmatlas.com', '/fetchPredictedStructure/MGYP002718891411.pdb', 'pdb'),
])
def test_notebook_identifiers_use_expected_fixed_providers(accession, host, path, file_format):
    from urllib.parse import urlparse
    with patch.object(runner, 'download_structure', return_value='structure') as download:
        text, actual_format, canonical = runner.resolve_query_id(accession)
    url = urlparse(download.call_args.args[0])
    assert (url.scheme, url.netloc, url.path) == ('https', host, path)
    assert text == 'structure' and actual_format == file_format
    assert canonical in path


def test_failure_is_not_reported_as_no_sites():
    with patch.object(runner.subprocess, 'Popen') as popen, patch.object(runner.os, 'killpg'):
        popen.return_value.wait.return_value = 1
        with pytest.raises(RuntimeError, match='failed'):
            runner.run_prediction(structure_text='ATOM  placeholder')


def test_successful_empty_result_and_cleanup():
    directories = []
    def start(*args, **kwargs):
        directory = Path(kwargs['cwd'])
        directories.append(directory)
        assert (directory / 'src/settings.py').exists()
        assert not (directory / '.env').exists()
        (directory / 'status.json').write_text(json.dumps({'has_sites': False}))
        from unittest.mock import Mock
        return Mock(wait=lambda **kw: 0, pid=123456)
    with patch.object(runner.subprocess, 'Popen', side_effect=start), patch.object(runner.os, 'killpg'):
        result = runner.run_prediction(structure_text='ATOM  placeholder')
    assert result['sites'] == [] and result['site_count'] == 0
    assert not directories[0].exists()


def test_timeout_kills_worker_group():
    with patch.object(runner.subprocess, 'Popen') as popen, patch.object(runner.os, 'killpg') as kill:
        popen.return_value.wait.side_effect = [subprocess.TimeoutExpired('worker', 180), 0]
        with pytest.raises(RuntimeError, match='180-second'):
            runner.run_prediction(structure_text='ATOM  placeholder')
        kill.assert_called_once()


def test_mcp_tool_metadata():
    import asyncio
    from integrations.mcp.server import mcp
    tools = asyncio.run(mcp.list_tools())
    assert len(tools) == 1
    tool = tools[0]
    assert tool.name == 'predict_zinc_binding_sites'
    assert tool.annotations.readOnlyHint is True
    assert 'pdb_id' in tool.inputSchema['properties']
    assert tool.inputSchema['properties']['minimum_probability']['default'] == 0.0001


@pytest.mark.parametrize('cutoff,expected_count', [(None, 2), (0, 3), (0.9, 1)])
def test_notebook_cutoff_and_explicit_override(cutoff, expected_count):
    def start(*args, **kwargs):
        directory = Path(kwargs['cwd'])
        table = directory / 'results/table/table_predicted_zn_binding_sites.csv'
        table.parent.mkdir(parents=True)
        with table.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=['prob', 'predicted_ion_pos',
                'chain_resi_1', 'chain_resi_2', 'chain_resi_3', 'chain_resi_4',
                'resi_comb_ordered', 'template_id', 'template_rmsd'])
            writer.writeheader()
            for probability in [0.000099, 0.0001, 0.99]:
                writer.writerow(dict(prob=probability, predicted_ion_pos='[1,2,3]',
                    chain_resi_1='A_1', chain_resi_2='A_2', chain_resi_3='A_3',
                    chain_resi_4='', resi_comb_ordered='{C,C,C}', template_id=14,
                    template_rmsd=0.1))
        (directory / 'status.json').write_text(json.dumps({'has_sites': True}))
        from unittest.mock import Mock
        return Mock(wait=lambda **kw: 0, pid=123456)
    options = {} if cutoff is None else {'minimum_probability': cutoff}
    with patch.object(runner.subprocess, 'Popen', side_effect=start), patch.object(runner.os, 'killpg'):
        result = runner.run_prediction(structure_text='ATOM  placeholder', max_sites=1, **options)
    assert result['minimum_probability'] == (0.0001 if cutoff is None else cutoff)
    assert result['site_count'] == expected_count
    assert result['truncated'] == (expected_count > 1)
    assert [site['probability'] for site in result['sites']] == [0.99]


def test_stdio_protocol_and_invalid_tool_call():
    import asyncio
    import sys
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    async def exercise():
        parameters = StdioServerParameters(command=sys.executable,
            args=['-m', 'integrations.mcp.server'], cwd=str(runner.ROOT))
        async with stdio_client(parameters) as (reader, writer):
            async with ClientSession(reader, writer) as session:
                await session.initialize()
                listing = await session.list_tools()
                assert listing.tools[0].name == 'predict_zinc_binding_sites'
                result = await session.call_tool('predict_zinc_binding_sites', {'pdb_id': '../bad'})
                assert result.isError
    asyncio.run(exercise())
