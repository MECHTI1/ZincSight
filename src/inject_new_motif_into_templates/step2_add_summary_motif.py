#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
auto_add_motif.py

Append a new motif across four representative CSV tables by parsing an mmCIF
using specific atom-selection rules. Pure Python only—no database.

This script runs with test inputs:
  CSV directory: editable_template_files_sql
  mmCIF file path: ../mmcif_with_ZN/5ROB.cif
  PDBID: 5ROB
  Binding residues: B_8,B_9,B_19

Outputs CSVs in the same format as the existing tables:
  - minimized_training_cluster_information.csv
  - training_representative_metal_sites_kruskal_v2.csv
  - motif_representative_coordinates_table_v2.csv
  - motif_representative_detailed_coordinates_table_v2.csv
"""
import os
import sys
import math
import json
from pathlib import Path

import pandas as pd
from functools import lru_cache
from Bio.PDB import MMCIFParser
from src.inject_new_motif_into_templates.Kruskal_Algorithm_V2 import main as kruskal_main

# Helper to parse Postgres-style array literals "{x,y,z}" to Python list of floats
def parse_array(s: str) -> list[float]:
    return [float(v) for v in s.strip('{}').split(',') if v]

# One-letter conversion
AA3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Atoms to include for detailed coordinates
ATOM_DICT = {
    'HIS': ['CG', 'CD2', 'ND1', 'CE1', 'NE2', 'CB'],
    'GLU': ['OE1', 'OE2', 'CD', 'CG'],
    'CYS': ['SG', 'CB'],
    'ASP': ['OD1', 'OD2', 'CG', 'CB'],
    'MET': ['SD', 'CG'],
    'THR': ['OG1', 'CG2', 'CA'],
    'TYR': ['OH', 'CZ'],
    'SER': ['OG', 'CB'],
    'ASN': ['ND2', 'OD1', 'CB']
}

# Close-atom rules
CLOSE_RULES = {
    'HIS': ('NE2', 'ND1'), 'GLU': ('OE1', 'OE2'), 'CYS': 'SG',
    'THR': 'OG1', 'ASP': ('OD1', 'OD2'), 'MET': 'SD',
    'TYR': 'OH', 'SER': 'OG', 'ASN': 'OD1'
}

# Far-atom rules
FAR_RULES = {
    'HIS': 'CG', 'GLU': 'CG', 'CYS': 'CB', 'THR': 'CA',
    'ASP': 'CB', 'MET': 'CG', 'TYR': 'CZ', 'SER': 'CB', 'ASN': 'CB'
}

@lru_cache(maxsize=None)
def get_atom_close_coordinates(residue):
    choice = CLOSE_RULES.get(residue.get_resname())
    if not choice:
        raise ValueError(f"No close-atom rule for {residue.get_resname()}")
    if isinstance(choice, tuple):
        points = [residue[c].get_coord() for c in choice if residue.has_id(c)]
        avg = sum(points) / len(points)
        bf = residue[choice[0]].get_bfactor()
        return avg.tolist(), bf
    atom = residue[choice]
    return atom.get_coord().tolist(), atom.get_bfactor()

@lru_cache(maxsize=None)
def get_far_atom_coordinates(residue):
    choice = FAR_RULES.get(residue.get_resname())
    if not choice:
        raise ValueError(f"No far-atom rule for {residue.get_resname()}")
    atom = residue[choice]
    return atom.get_coord().tolist()


def create_all_relevant_atoms_from_residues_list(structure):
    sel = {}
    for model in structure:
        for chain in model:
            for res in chain:
                key = f"{chain.id}_{res.get_id()[1]}"
                name = res.get_resname()
                if name in ATOM_DICT:
                    coords = {}
                    for atm in ATOM_DICT[name]:
                        if res.has_id(atm):
                            coords[atm] = res[atm].get_coord().tolist()
                    sel[key] = coords
    return sel


def main(pdb_path,pdbid,residues):
    residues = residues.split(',')
    csv_dir='editable_template_files_sql'
    # PathsP
    sum_fp   =Path(__file__).parent /"editable_template_files_sql"/ 'minimized_training_cluster_information.csv'
    ii_fp    = Path(__file__).parent /"editable_template_files_sql"/'training_representative_metal_sites_kruskal_v2.csv'
    coord_fp =Path(__file__).parent /"editable_template_files_sql"/'motif_representative_coordinates_table_v2.csv'
    det_fp   = Path(__file__).parent /"editable_template_files_sql"/ 'motif_representative_detailed_coordinates_table_v2.csv'

    # Load tables
    df_sum   = pd.read_csv(sum_fp)
    df_ii    = pd.read_csv(ii_fp)
    df_coord = pd.read_csv(coord_fp)
    df_det   = pd.read_csv(det_fp)

    # Parse structure
    structure = MMCIFParser(QUIET=True).get_structure(pdbid, pdb_path)
    detailed_atoms = create_all_relevant_atoms_from_residues_list(structure)
    zn_atoms = [a for a in structure.get_atoms() if a.element == 'ZN']
    if not zn_atoms:
        print("No Zn atoms found", file=sys.stderr)
        sys.exit(1)

    # New motif site_id
    new_site = df_sum['site_id'].max() + 1
    new_site_training_representative =df_ii['site_id'].max() + 1
    # Build coordinate rows
    codes, coord_rows, det_rows = [], [], []
    for cr in residues:
        chain, resi = cr.split('_')
        resi = int(resi)
        res_obj = None
        for model in structure:
            if chain in model:
                for r in model[chain]:
                    if r.get_id()[1] == resi:
                        res_obj = r
                        break
            if res_obj:
                break
        if not res_obj:
            print(f"Residue {cr} not found", file=sys.stderr)
            sys.exit(1)
        aa1 = AA3TO1[res_obj.get_resname()]
        codes.append(aa1)
        close_c, bf = get_atom_close_coordinates(res_obj)
        far_c = get_far_atom_coordinates(res_obj)
        coord_rows.append({
            'pdbid': pdbid.lower(),
            'chain_resi': cr,
            'close_atom_coord': '{' + ','.join(f"{x:.3f}" for x in close_c) + '}',
            'far_atom_coord': '{' + ','.join(f"{x:.3f}" for x in far_c) + '}',
            'b_factor': int(bf)
        })
        det_rows.append({
            'pdbid': pdbid.lower(),
            'chain_resi': cr,
            'resi_type': res_obj.get_resname(),
            'dict_atom_coord': json.dumps(detailed_atoms.get(cr, {}))
        })

    # Append to summary
    resi_comb = '{' + ','.join(codes) + '}'
    df_sum = pd.concat([df_sum, pd.DataFrame([{
        'site_id': new_site,
        'pdbid': pdbid.lower(),
        'resi_comb': resi_comb,
        'num_of_residues': len(codes),
        'cluster_size': 1
    }])], ignore_index=True)

    # Append coords
    df_coord = pd.concat([df_coord, pd.DataFrame(coord_rows)], ignore_index=True)
    df_det   = pd.concat([df_det,   pd.DataFrame(det_rows)], ignore_index=True)

    # Step 1: Build all pairwise edges
    all_edges = []
    for i in range(len(residues)):
        for j in range(i+1, len(residues)):
            r1, r2 = residues[i], residues[j]
            c1 = parse_array(df_coord[df_coord['chain_resi']==r1].iloc[-1]['close_atom_coord'])
            c2 = parse_array(df_coord[df_coord['chain_resi']==r2].iloc[-1]['close_atom_coord'])
            f1 = parse_array(df_coord[df_coord['chain_resi']==r1].iloc[-1]['far_atom_coord'])
            f2 = parse_array(df_coord[df_coord['chain_resi']==r2].iloc[-1]['far_atom_coord'])
            cd = round(math.dist(c1, c2), 3)
            fd = round(math.dist(f1, f2), 3)
            aa1 = AA3TO1[df_det[df_det['chain_resi']==r1].iloc[-1]['resi_type']]
            aa2 = AA3TO1[df_det[df_det['chain_resi']==r2].iloc[-1]['resi_type']]
            aa_pair = ''.join(sorted([aa1, aa2]))
            all_edges.append([r1, r2, cd, fd, aa_pair])

    # Step 2: Apply Kruskal algorithm
    kruskal_input = [[e[0], e[1], e[2]] for e in all_edges]  # use close_atom_dis as weight
    kruskal_result = kruskal_main(residues, kruskal_input)

    # Step 3: Filter only selected edges from full list
    kruskal_edges = []
    for r1, r2, cd in kruskal_result:
        for e in all_edges:
            if {e[0], e[1]} == {r1, r2}:
                kruskal_edges.append((
                    new_site_training_representative,
                    None,
                    pdbid.lower(),
                    e[4],  # aa_pair
                    e[2],  # close_atom_dis
                    e[3],  # far_atom_dis
                    e[0],
                    e[1]
                ))
                break

    # Step 4: Reorder with overlapping path

    # Step 5: Add ii_order_id
    kruskal_edges = [
        dict(zip(['site_id', 'ii_order_id', 'pdbid', 'aa_pair', 'close_atom_dis', 'far_atom_dis', 'chain_resi_1', 'chain_resi_2'], [*e[:1], i+1, *e[2:]]))
        for i, e in enumerate(kruskal_edges)
    ]

    # Append training rows
    df_ii = pd.concat([df_ii, pd.DataFrame(kruskal_edges)], ignore_index=True)


    # Append training rows

    # Save
    df_sum.to_csv(sum_fp, index=False)
    df_ii.to_csv(ii_fp, index=False)
    df_coord.to_csv(coord_fp, index=False)
    df_det.to_csv(det_fp, index=False)

    print(f"Motif site {new_site} added to {csv_dir}")

if __name__ == '__main__':
    csv_dir: str = 'editable_template_files_sql'
    cif_dir: str = '../mmcif_with_ZN/5ROB.cif'
    pdbid: str = '5ROB'
    residues: str = 'B_8,B_9,B_19'

    main(csv_dir,cif_dir, pdbid,residues)
