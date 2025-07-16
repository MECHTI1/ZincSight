import sys
from argparse import ArgumentParser
from pathlib import Path

from pydash import py_, is_empty

from main_execute import execute_zincsight
from src.settings import QUERY_STRUCTURES_DIR

if __name__ == '__main__':
    parser = ArgumentParser(epilog='At least one of IDENTIFIERS or FILE is required')
    parser.add_argument('-r', '--rotamers', action='store_true', help='Include histidine rotamers')
    parser.add_argument('-i', '--identifiers', help='Comma-delimited list of identifiers')
    parser.add_argument('-f', '--file', help='Path to text file containing identifier per line')

    args = parser.parse_args()

    if not args.identifiers and not args.file:
        parser.print_help()
        sys.exit(1)

    ids = [] if not args.identifiers else args.identifiers.split(',')

    if args.file:
        path = Path(args.file)
        if not path.exists():
            print(f'Input file "{args.file}" not found')
            sys.exit(1)
        ids.extend(py_(path.read_text('utf8').split('\n')).map(lambda x: x.strip()).reject(is_empty).value())

    include_histidine_rotamers = args.rotamers
    project_dir = Path(__file__).parent
    query_structures_path = QUERY_STRUCTURES_DIR
    output_path = str((project_dir / 'results').absolute())

    structure_ids_for_download = ','.join(ids)

    compressed_results_path = execute_zincsight(
        include_histidine_rotamers,
        structure_ids_for_download,
        query_structures_path,
        output_path
    )

    print(f" Compressed predictions path: {compressed_results_path}")