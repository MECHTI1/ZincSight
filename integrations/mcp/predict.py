"""Run one prediction directly, without an MCP client or a hosted service."""
import argparse
import json
from pathlib import Path
from .runner import DEFAULT_MINIMUM_PROBABILITY, MAX_BYTES, run_prediction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--pdb-id')
    source.add_argument('--query-id', help='PDB/assembly, UniProt, full AlphaFold model or MGYP identifier')
    source.add_argument('--structure', type=Path)
    parser.add_argument('--format', choices=['pdb', 'cif'], default='pdb')
    parser.add_argument('--database', choices=['native', 'pglite'], default='native')
    parser.add_argument('--histidine-rotamers', action='store_true')
    parser.add_argument('--max-sites', type=int, default=100)
    parser.add_argument('--minimum-probability', type=float, default=DEFAULT_MINIMUM_PROBABILITY,
                        help='Post-prediction cutoff (default: 0.0001, as in the public notebook; 0 keeps all sites)')
    args = parser.parse_args()
    text = None
    if args.structure:
        if args.structure.stat().st_size > MAX_BYTES:
            parser.error('Structure exceeds the 5 MB limit.')
        text = args.structure.read_text()
    result = run_prediction(pdb_id=args.pdb_id, structure_text=text, file_format=args.format,
                            histidine_rotamers=args.histidine_rotamers, database_backend=args.database,
                            query_id=args.query_id, max_sites=args.max_sites,
                            minimum_probability=args.minimum_probability)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
