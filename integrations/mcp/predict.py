"""Run one prediction directly, without an MCP client or a hosted service."""
import argparse
import json
from pathlib import Path
from .runner import MAX_BYTES, run_prediction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--pdb-id')
    source.add_argument('--structure', type=Path)
    parser.add_argument('--format', choices=['pdb', 'cif'], default='pdb')
    parser.add_argument('--database', choices=['native', 'pglite'], default='native')
    parser.add_argument('--histidine-rotamers', action='store_true')
    args = parser.parse_args()
    text = None
    if args.structure:
        if args.structure.stat().st_size > MAX_BYTES:
            parser.error('Structure exceeds the 2 MB limit.')
        text = args.structure.read_text()
    result = run_prediction(pdb_id=args.pdb_id, structure_text=text, file_format=args.format,
                            histidine_rotamers=args.histidine_rotamers, database_backend=args.database)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
