#!/usr/bin/env python3
"""Pinned, isolated setup and execution for the experimental ZincSight skill."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

REVISION = 'ac1fe1e96e8b2ffe4609ab33f9beff579c76ca07'
SOURCE_SHA256 = 'df3ca7a4dd0fdd66a2efbde49d906c1aa09e9c636a90ed7ad20db16f15bc32dd'
SOURCE_URL = f'https://codeload.github.com/MECHTI1/ZincSight/tar.gz/{REVISION}'
LOCK = Path(__file__).with_name('requirements-lock.txt')


def run(command, **kwargs):
    subprocess.run(command, check=True, stdout=sys.stderr, stderr=sys.stderr, **kwargs)


def setup():
    python = sys.executable if sys.version_info[:2] == (3, 12) else shutil.which('python3.12')
    if not python:
        raise RuntimeError('Python 3.12 is required for this tested dependency set.')
    node, npm = shutil.which('node'), shutil.which('npm')
    if not node or not npm:
        raise RuntimeError('Node.js 22 or newer and npm are required.')
    version = subprocess.check_output([node, '--version'], text=True).strip()
    if int(version.removeprefix('v').split('.')[0]) < 22:
        raise RuntimeError(f'Node.js 22 or newer is required; found {version}.')
    base = Path(os.environ.get('ZINCSIGHT_CACHE_DIR', str(Path.home() / '.cache/zincsight-preview')))
    base.mkdir(parents=True, exist_ok=True)
    root = base / REVISION
    root.mkdir(exist_ok=True)
    source, venv = root / 'source', root / 'venv'
    ready = root / 'ready.json'
    wanted = {'revision': REVISION, 'requirements_sha256': hashlib.sha256(LOCK.read_bytes()).hexdigest()}
    with (root / 'setup.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if ready.exists() and json.loads(ready.read_text()) == wanted and (venv / 'bin/python').exists():
            return source, venv / 'bin/python'
        print('Preparing the isolated ZincSight Preview environment...', file=sys.stderr, flush=True)
        if not source.exists():
            with tempfile.TemporaryDirectory(prefix='source-', dir=root) as temporary:
                stage = Path(temporary)
                with urllib.request.urlopen(SOURCE_URL, timeout=60) as response:
                    data = response.read(30_000_001)
                if len(data) > 30_000_000 or hashlib.sha256(data).hexdigest() != SOURCE_SHA256:
                    raise RuntimeError('Pinned source archive failed size/hash verification.')
                archive = stage / 'source.tar.gz'
                archive.write_bytes(data)
                with tarfile.open(archive) as package:
                    package.extractall(stage, filter='data')
                shutil.move(str(stage / f'ZincSight-{REVISION}'), source)
        if not (venv / 'bin/python').exists():
            # A restored cloud cache may contain dangling or cyclic interpreter
            # links. Recreate only this generated environment before installing.
            run([python, '-m', 'venv', '--clear', str(venv)])
        environment = os.environ.copy()
        if shutil.which('gcc'):
            environment.setdefault('CC', 'gcc')
        run([str(venv / 'bin/python'), '-m', 'pip', 'install', '-r', str(LOCK)], env=environment, timeout=600)
        run([npm, 'ci', '--prefix', str(source / 'integrations/mcp/pglite'),
             '--ignore-scripts', '--no-audit', '--no-fund'], timeout=180)
        run([str(venv / 'bin/python'), '-c', 'import pymol, psycopg2, Bio, scipy, sklearn'], timeout=30)
        ready.write_text(json.dumps(wanted))
    return source, venv / 'bin/python'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--query-id')
    group.add_argument('--structure', type=Path)
    parser.add_argument('--format', choices=['pdb', 'cif'])
    parser.add_argument('--histidine-rotamers', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--minimum-probability', type=float, default=0.0001,
                        help='Notebook cutoff; set 0 to retain all computed sites')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--setup-only', action='store_true')
    args = parser.parse_args()
    if not args.setup_only and not (args.query_id or args.structure):
        parser.error('Supply --query-id or --structure.')
    if args.structure and (not args.structure.is_file() or args.structure.stat().st_size > 5_000_000):
        parser.error('The input must be an existing structure file of at most 5 MB.')
    source, python = setup()
    if args.setup_only:
        result = {'ready': True, 'software_revision': REVISION}
    else:
        command = [str(python), '-m', 'integrations.mcp.predict', '--database', 'pglite', '--max-sites', '100',
                   '--minimum-probability', str(args.minimum_probability)]
        if args.query_id:
            command += ['--query-id', args.query_id]
        else:
            file_format = args.format or ('cif' if args.structure.suffix.lower() in ('.cif', '.mmcif') else 'pdb')
            command += ['--structure', str(args.structure.resolve()), '--format', file_format]
        if args.histidine_rotamers:
            command += ['--histidine-rotamers']
        completed = subprocess.run(command, cwd=source, text=True, capture_output=True, timeout=240)
        if completed.returncode:
            print(completed.stderr[-5000:], file=sys.stderr)
            raise RuntimeError('ZincSight prediction failed; no scientific result was produced.')
        result = json.loads(completed.stdout)
        result['software_revision'] = REVISION
        result['validation_status'] = 'Experimental integration; historical notebook-output differences remain unresolved.'
    content = json.dumps(result, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content + '\n')
    print(content)


if __name__ == '__main__':
    main()
