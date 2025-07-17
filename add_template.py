import argparse
from src.inject_new_motif_into_templates.main_add_templates import main as add_template

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run full motif injection pipeline'
    )
    parser.add_argument('--pdbid', default='4f3w', help='PDB ID for new motif')
    parser.add_argument('--residues', default='B_89,B_56,B_92',
                        help='Comma-separated residues, e.g. B_8,B_9,B_19')
    parser.add_argument('--output', default='dump.sql', help='Final SQL dump path')
    args = parser.parse_args()
    add_template(
        pdbid=args.pdbid,
        residues=args.residues,
    )
