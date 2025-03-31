import sys
sys.path.append('/content/ZincSight')
from main_execute import execute_zincsight


structure_ids_for_download = ' ' # as 8QEP, P0A6G5, AF-A0A068N621-F1-v4, MGYP002718891411'
query_structures_path = '' # Path to input structures
include_histidine_rotamers = True #True is default
output_path = ''  # Path to save output


compressed_results_path = execute_zincsight(
    include_histidine_rotamers,       # Boolean flag
    structure_ids_for_download,       # String of structure IDs
    query_structures_path,  # Input directory path
    output_path             # Output directory path
)

print (f" Compressed predictions path: {compressed_results_path}")