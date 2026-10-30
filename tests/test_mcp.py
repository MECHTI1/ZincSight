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
