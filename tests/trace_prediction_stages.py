"""Diagnostic ablation of state and ordering fixes on public notebook inputs.

This writes scientific intermediate tables, not a biological accuracy benchmark.
All source modifications are confined to disposable copies. No private reference
is read or uploaded. Run native PostgreSQL on a supported non-root host.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
STATE_FILE = 'src/dis_angles_calculations/with_rotations/return_nonredundant_dicts_his_comb.py'
ORDER_FILES = [
    'src/ii_search.py', 'src/motif_search_primary.py', 'src/primary_script.py',
    'src/scoring_and_compression/compress_table_by_proximity.py',
    'src/dis_angles_calculations/with_rotations/dis_angles_calc_primary_generator_with_his_rotation_180deg.py',
    'src/dis_angles_calculations/without_rotations/dis_angles_calc_primary_generator_without_his_rotation.py',
]

CHILD = '''
import contextlib, json, os
from pathlib import Path

config = json.loads(Path('config.json').read_text())
output = Path(config['output'])
if config['backend'] == 'pglite':
    from embedded_database import embedded_database
    database = embedded_database(Path.cwd(), Path(config['launcher']))
else:
    database = contextlib.nullcontext()

with database:
    from src import primary_script as pipeline
    from src.settings import get_db_connection
    template = [None]
    from src.dis_angles_calculations.with_rotations import (
        dis_angles_calc_primary_generator_with_his_rotation_180deg as rotamer_pipeline,
        return_nonredundant_dicts_his_comb as exclusions,
    )
    worker_inputs = output / 'worker-inputs'
    worker_inputs.mkdir()
    original_worker = rotamer_pipeline.worker
    def traced_worker(item):
        match_id, residues = item
        record = {'template_id': template[0], 'match_id': match_id,
                  'residue_order': list(residues), 'residues': residues,
                  'excluded_indices_before': sorted(getattr(
                      exclusions, 'processed_indices_to_move_out_dict_1', set()))}
        with (worker_inputs / (str(template[0]) + '-' + str(os.getpid()) + '.jsonl')).open('a') as stream:
            stream.write(json.dumps(record) + '\\n')
        return original_worker(item)
    rotamer_pipeline.worker = traced_worker

    def snapshot(stage, query):
        conn = get_db_connection()
        cur = conn.cursor()
        cur.execute(query, (tuple(config['names']),))
        names = [column[0] for column in cur.description]
        rows = [dict(zip(names, row)) for row in cur.fetchall()]
        with (output / 'snapshot-counts.jsonl').open('a') as stream:
            stream.write(json.dumps({'stage': stage, 'template_id': template[0], 'rows': len(rows)}) + '\\n')
        with (output / (stage + '.jsonl')).open('a') as stream:
            for row in rows:
                row['trace_template_id'] = template[0]
                stream.write(json.dumps(row) + '\\n')
        cur.close()
        conn.close()

    first = pipeline.main_first_step
    def search(site_id):
        template[0] = site_id
        first(site_id)
        snapshot('01-candidates', """SELECT f.*, r.rmsd_overall, r.rmsd_close_atom
            FROM af_dataset_final_motif_search_table_v2 f
            LEFT JOIN af_dataset_rmsd_of_matches_table_v2 r USING (match_id)
            WHERE f.pdbid_alphafoldmodel IN %s""")
    pipeline.main_first_step = search

    second = pipeline.main_second_step
    def reduce():
        second()
        snapshot('02-deduplicated', """SELECT * FROM af_dataset_final_motif_search_table_v2
            WHERE pdbid_alphafoldmodel IN %s""")
    pipeline.main_second_step = reduce

    third = pipeline.main_third_step_his_rotation_180deg
    def geometry(cores):
        third(cores)
        snapshot('03-geometry', """SELECT * FROM af_dataset_final_motif_search_table_v2
            WHERE pdbid_alphafoldmodel IN %s""")
    pipeline.main_third_step_his_rotation_180deg = geometry

    third_without_rotamers = pipeline.main_third_step_without_his_rotation
    def geometry_without_rotamers():
        third_without_rotamers()
        snapshot('03-geometry', """SELECT * FROM af_dataset_final_motif_search_table_v2
            WHERE pdbid_alphafoldmodel IN %s""")
    pipeline.main_third_step_without_his_rotation = geometry_without_rotamers

    score = pipeline.final_scoring_and_insertion_to_table
    def scoring(paths):
        score(paths)
        template[0] = None
        snapshot('04-scored', """SELECT * FROM scored_af_dataset_with_aggregated_final_tables
            WHERE pdbid_alphafoldmodel IN %s""")
    pipeline.final_scoring_and_insertion_to_table = scoring

    compress = pipeline.table_compression
    def compression(conn):
        compress(conn)
        conn.commit()  # Make the completed stage visible to the read-only trace.
        snapshot('05-compressed', """SELECT * FROM final_compressed_table_with_scored_binding_sites
            WHERE pdbid_alphafoldmodel IN %s""")
    pipeline.table_compression = compression

    paths = [str(Path('query') / (name + '.cif')) for name in config['names']]
    if config['backend'] == 'pglite':
        pipeline.main(paths, config['histidine_rotamers'], str(output), 2, False)
    else:
        from main_execute import execute_zincsight
        execute_zincsight(config['histidine_rotamers'], '', 'query', str(output), 2, False)
'''


def trace(original, current, inputs, output, backend, variant, histidine_rotamers=True):
    from notebook_regression import SOURCES
    destination = output / variant
    destination.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='zincsight-stage-trace-') as directory:
        work = Path(directory)
        shutil.copytree(original / 'src', work / 'src', ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(original / 'main_execute.py', work)
        if variant in ('state-fix', 'current'):
            shutil.copy2(current / STATE_FILE, work / STATE_FILE)
        if variant in ('ordering', 'current'):
            for name in ORDER_FILES:
                shutil.copy2(current / name, work / name)
        if backend == 'pglite':
            # Required lifecycle adaptation only; no ordering/state change here.
            for name in ORDER_FILES[-2:]:
                path = work / name
                content = path.read_text()
                if 'DROP TABLE temp_final_motif' not in content:
                    content = content.replace('    # Commit changes\n',
                        '    cur.execute("DROP TABLE temp_final_motif")\n    # Commit changes\n', 1)
                    assert 'DROP TABLE temp_final_motif' in content
                    path.write_text(content)
        shutil.copy2(current / 'integrations/mcp/embedded_database.py', work)
        (work / 'query').mkdir()
        manifest = []
        for name in SOURCES:
            source = inputs / (name + '.cif')
            shutil.copy2(source, work / 'query')
            manifest.append({'id': name, 'sha256': hashlib.sha256(source.read_bytes()).hexdigest()})
        config = {'output': str(destination), 'backend': backend, 'names': list(SOURCES),
                  'histidine_rotamers': histidine_rotamers,
                  'launcher': str(current / 'integrations/mcp/pglite/server.mjs')}
        (work / 'config.json').write_text(json.dumps(config))
        (work / 'run.py').write_text(CHILD)
        (destination / 'inputs.json').write_text(json.dumps(manifest, indent=2))
        (destination / 'settings.json').write_text(json.dumps({
            'variant': variant, 'backend': backend, 'histidine_rotamers': histidine_rotamers,
            'minimum_probability': 0.0001, 'cores': 2,
            'source_sha256': {str(p.relative_to(work)): hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (work / 'src').rglob('*.py')},
        }, indent=2))
        env = {**os.environ, 'PYTHON_DOTENV_DISABLED': '1', 'DEBUGGING_TEMP': 'false',
               'KEEP_TEMP_TABLES_TEMP': 'false', 'FILTER_PROB_AFTER_COMPRESSION': 'true',
               'EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD': '0.0001', 'OMP_NUM_THREADS': '1',
               'OPENBLAS_NUM_THREADS': '1', 'MPLBACKEND': 'Agg'}
        print('TRACE', backend, variant, flush=True)
        with (destination / 'run.log').open('w') as log:
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
            raise RuntimeError(f'Trace failed: {destination / "run.log"}')
        from collections import Counter
        with (destination / 'table/table_predicted_zn_binding_sites.csv').open() as stream:
            counts = dict(Counter(row['structure_id'] for row in csv.DictReader(stream)))
        (destination / 'counts.json').write_text(json.dumps(counts, indent=2))
        print(variant, counts, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-root', type=Path, required=True)
    parser.add_argument('--current-root', type=Path, default=ROOT)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--backend', choices=['native', 'pglite'], default='native')
    parser.add_argument('--histidine-rotamers', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--variants', nargs='+', choices=['original', 'state-fix', 'ordering', 'current'],
                        default=['original', 'state-fix', 'ordering', 'current'])
    args = parser.parse_args()
    for variant in args.variants:
        trace(args.original_root.resolve(), args.current_root.resolve(), args.input_dir.resolve(),
              args.output_dir.resolve(), args.backend, variant, args.histidine_rotamers)


if __name__ == '__main__':
    main()
