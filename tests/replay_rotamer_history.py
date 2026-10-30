"""Replay public trace candidates to check dependence on earlier calculations.

The geometry code must be identical between source roots. Only the original
versus corrected exclusion helper is exchanged. This uses no archived result
or expected scientific score and requires no database after trace collection.
"""
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_FILE = 'src/dis_angles_calculations/with_rotations/create_dict_resi_dis_hisangles_zncoord_rot180deg.py'
EXCLUSIONS_FILE = 'src/dis_angles_calculations/with_rotations/return_nonredundant_dicts_his_comb.py'


def difference(first, second):
    if not isinstance(first, dict) or not isinstance(second, dict):
        return {'error': True, 'changed': first != second}
    position = math.dist(first['metalcoord'], second['metalcoord'])
    first_angles = first['HIS_anlges_stats']['candidate_point_angles_stast']
    second_angles = second['HIS_anlges_stats']['candidate_point_angles_stast']
    angle = max(abs(float(first_angles[k]) - float(second_angles[k])) for k in first_angles)
    return {'error': False, 'changed': position > 1e-8 or angle > 1e-8,
            'position_delta_angstrom': position, 'max_angle_delta': angle}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trace-dir', type=Path, required=True)
    parser.add_argument('--original-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--template-id', type=int, default=15)
    parser.add_argument('--max-warmups', type=int, default=20)
    args = parser.parse_args()
    assert (ROOT / GEOMETRY_FILE).read_bytes() == (args.original_root / GEOMETRY_FILE).read_bytes(), \
        'This experiment isolates the exclusion helper; geometry implementations must match.'
    sys.path.insert(0, str(ROOT))
    from src.dis_angles_calculations.with_rotations import create_dict_resi_dis_hisangles_zncoord_rot180deg as geometry
    from src.dis_angles_calculations.with_rotations import return_nonredundant_dicts_his_comb as corrected
    spec = importlib.util.spec_from_file_location('original_exclusions', args.original_root / EXCLUSIONS_FILE)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    records = {}
    for path in (args.trace_dir / 'worker-inputs').glob(f'{args.template_id}-*.jsonl'):
        for line in path.read_text().splitlines():
            record = json.loads(line)
            records[record['match_id']] = record
    assert records, 'No worker inputs were captured for the selected template.'
    candidates = {}
    for line in (args.trace_dir / '02-deduplicated.jsonl').read_text().splitlines():
        record = json.loads(line)
        if record['trace_template_id'] == args.template_id:
            candidates[record['match_id']] = record
    warmups = [records[key] for key in sorted(records)[:args.max_warmups]]
    report = {'template_id': args.template_id, 'warmup_count': len(warmups), 'candidates': []}
    original_callback = geometry.create_list_only_new_dicts_his_rot
    try:
        for match_id, candidate in sorted(candidates.items()):
            target = records[match_id]['residues']
            result = {'structure': candidate['pdbid_alphafoldmodel'], 'match_id': match_id,
                      'residues': sorted(target), 'modes': {}}
            for label, callback in [('original', original.main), ('corrected', corrected.main)]:
                geometry.create_list_only_new_dicts_his_rot = callback
                original.processed_indices_to_move_out_dict_1.clear()
                baseline = geometry.caclulate_dis_CoordinationAngles_HISangles_stats(copy.deepcopy(target))
                trials = []
                for warmup in warmups:
                    original.processed_indices_to_move_out_dict_1.clear()
                    geometry.caclulate_dis_CoordinationAngles_HISangles_stats(copy.deepcopy(warmup['residues']))
                    before = sorted(original.processed_indices_to_move_out_dict_1)
                    actual = geometry.caclulate_dis_CoordinationAngles_HISangles_stats(copy.deepcopy(target))
                    trials.append({'warmup_match_id': warmup['match_id'],
                                   'excluded_indices_before': before, **difference(baseline, actual)})
                result['modes'][label] = trials
            report['candidates'].append(result)
    finally:
        geometry.create_list_only_new_dicts_his_rot = original_callback
    summary = {}
    for mode in ('original', 'corrected'):
        trials = [trial for candidate in report['candidates'] for trial in candidate['modes'][mode]]
        summary[mode] = {'candidate_count': len(report['candidates']), 'trial_count': len(trials),
                         'changed_trials': sum(t['changed'] for t in trials),
                         'error_trials': sum(t['error'] for t in trials),
                         'max_position_delta_angstrom': max(t.get('position_delta_angstrom', 0) for t in trials)}
    report['summary'] = summary
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))
    print(json.dumps(summary, indent=2))
    assert summary['original']['changed_trials'] > 0, 'Original state dependence was not reproduced.'
    assert summary['corrected']['changed_trials'] == 0 and summary['corrected']['error_trials'] == 0, \
        'Corrected geometry must not depend on the preceding candidate.'


if __name__ == '__main__':
    main()
