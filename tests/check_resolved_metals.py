"""A small coordinate sanity check, not an accuracy or calibration benchmark.

Match predictions one-to-one to the resolved transition-metal ions in three
public notebook examples. Unmatched predictions are not classified as false
positives: absence of a resolved ion does not establish absence of binding.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def check(input_dir, table, minimum_probability=0.0001, maximum_distance=1.0):
    with table.open() as stream:
        predictions = list(csv.DictReader(stream))
    report = {'minimum_probability': minimum_probability,
              'maximum_distance_angstrom': maximum_distance, 'structures': []}
    for name in ('1KLS', '2A0S-assembly1', '8QEP'):
        model = MMCIFParser(QUIET=True).get_structure(name, input_dir / f'{name}.cif')[0]
        ions = [atom for atom in model.get_atoms()
                if atom.element in {'ZN', 'FE', 'MN', 'CO', 'CU', 'NI'}]
        sites = [row for row in predictions if row['structure_id'] == name
                 and float(row['prob']) >= minimum_probability]
        if not ions or len(sites) < len(ions):
            raise AssertionError(f'{name}: missing resolved ions or too few predictions')
        distances = cdist(np.asarray([ion.coord for ion in ions]),
                          np.asarray([json.loads(site['predicted_ion_pos']) for site in sites]))
        ion_indices, site_indices = linear_sum_assignment(distances)
        matches = []
        for i, j in zip(ion_indices, site_indices):
            atom, site = ions[i], sites[j]
            matches.append({'chain': atom.get_parent().get_parent().id,
                            'residue': atom.get_parent().id[1], 'element': atom.element,
                            'distance_angstrom': float(distances[i, j]),
                            'probability': float(site['prob'])})
        report['structures'].append({'structure': name, 'predictions': len(sites),
                                     'resolved_ions': len(ions), 'matches': matches})
    report['passed'] = all(match['distance_angstrom'] <= maximum_distance
                           for item in report['structures'] for match in item['matches'])
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = check(args.input_dir, args.table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report))
    if not report['passed']:
        raise SystemExit('A resolved-metal coordinate check exceeded 1 angstrom')
