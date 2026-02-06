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
    parser.add_argument('-c','--cores', type=int, help='Number of cores', default=2)
    # NEW: add boolean flag
    parser.add_argument('-p','--pse_output', action='store_true',help='If set, generate PyMOL .pse session output')



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
    create_pse_output=args.pse_output
    project_dir = Path(__file__).parent
    query_structures_path = QUERY_STRUCTURES_DIR
    output_path = str((project_dir / 'results').absolute())

    structure_ids_for_download = ','.join(ids)

    num_cores= args.cores

    compressed_results_path = execute_zincsight(
        include_histidine_rotamers,
        structure_ids_for_download,
        query_structures_path,
        output_path,
        num_cores,
        create_pse_output
    )

    print(f" Compressed predictions path: {compressed_results_path}")