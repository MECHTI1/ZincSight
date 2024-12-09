import os
import subprocess
from datetime import datetime

def create_results_tarfile(output_file_path, folder_names, path_output):
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

    # Create the tar.gz file using tar and pigz for the specified folders within path_output
    with open(output_file_path, 'wb') as f_out:
        # Switch to path_output and archive each specified folder by its name, avoiding full paths
        tar_command = ['tar', 'cf', '-', '-C', path_output] + folder_names
        p1 = subprocess.Popen(tar_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['pigz'], stdin=p1.stdout, stdout=f_out)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
    print(f'Tar file created: {output_file_path}')


def compress_unified_results(sample_id, his_rot_sampling, path_output):
    """
    Compress specified input directories within path_output and save to the output directory.

    Parameters:
        sample_id (str): Identifier for the sample.
        his_rot_sampling (bool): Whether rotation sampling is enabled.
    """
    # Define input directories to include within path_output
    dirs_to_compress = ['structures_with_predicted_zn', 'table']
    folder_suffix = "_".join(dirs_to_compress)  # Join folder names for filename suffix

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
        f"ZincSight_{sample_id}{str_his_rot_sampling}_{folder_suffix}_{date}.tar.gz"
    )

    # Call the function to create the tar file with specified folder names (not full paths)
    create_results_tarfile(output_file_path, dirs_to_compress, path_output)

    # Return the full path to the created tar.gz file
    return output_file_path

if __name__ == "__main__":
    from src.settings import RESULTS_DIR
    path_output=RESULTS_DIR
    sample_id = "input_id"
    his_rot_sampling_boolean = True
    compress_unified_results(sample_id, his_rot_sampling_boolean, path_output)