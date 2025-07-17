import argparse
from src.inject_new_motif_into_templates.main_add_templates import main as add_template

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run full motif injection pipeline'
    )
    parser.add_argument('pdb_id', help='PDB ID with a new motif')
    parser.add_argument('residues', help='Comma-separated CHAIN_RESIDUE pairs, e.g. B_8,B_9,B_19')
    args = parser.parse_args()
    add_template(
        pdbid=args.pdb_id,
        residues=args.residues,
    )
