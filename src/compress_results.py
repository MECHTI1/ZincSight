import os
import subprocess
import tarfile
from datetime import datetime
from pathlib import Path


def create_results_tarfile(output_file_path, dir_paths, reference_csv_path):
    """
    Create a .tar.gz file containing only the specified folders within path_output.

    Parameters:
        output_file_path (str): Path for the output .tar.gz file.
        folder_names (list of str): List of folder names within path_output to include in the archive.
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with tarfile.open(output_file_path, 'w:gz') as tf:
        for dir_path in dir_paths:
            tf.add(dir_path, dir_path.name)

        tf.add(reference_csv_path, reference_csv_path.name)

    print(f'Tar file created: {output_file_path}')


def compress_unified_results(sample_id, his_rot_sampling, path_output):
    """
    Compress specified input directories within path_output and save to the output directory.

    Parameters:
        sample_id (str): Identifier for the sample.
        his_rot_sampling (bool): Whether rotation sampling is enabled.
    """
    # Define input directories to include within path_output
    dirs_to_compress = [Path(path_output) / 'structures_with_predicted_zn', Path(path_output) / 'table']
    dirs_to_compress = [p for p in dirs_to_compress if p.exists()]
    # Define output directory
    output_dir = os.path.join(path_output, "compressed_results")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Format date safely for filenames
    date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    str_his_rot_sampling = "r-sampled" if his_rot_sampling else "r-unsampled"

    if sample_id:
        sample_id += "_"

    # Define path for the compressed resulted file , such that the file will be within the output directory
    output_file_path = os.path.join(
        output_dir,
        f"ZincSight_{sample_id}{str_his_rot_sampling}_{date}.tar.gz"
    )

    # **Path to the reference CSV file (assumed to be in the same directory as the script)**
    reference_csv_path = Path(__file__).parent / 'reference_motif_templates.csv'

    # Call the function to create the tar file with specified folder names (not full paths)
    create_results_tarfile(output_file_path, dirs_to_compress, reference_csv_path)

    # Return the full path to the created tar.gz file
    return output_file_path

if __name__ == "__main__":
    from src.settings import RESULTS_DIR
    path_output=RESULTS_DIR
    sample_id = "input_id"
    his_rot_sampling_boolean = True
    compress_unified_results(sample_id, his_rot_sampling_boolean, path_output)