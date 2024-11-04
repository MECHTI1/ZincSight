import os
import subprocess
from src.settings import RESULTS_DIR
from datetime import datetime

def create_results_tarfile(output_filename, folder_names):
    """
    Create a .tar.gz file containing only the specified folders within RESULTS_DIR.

    Parameters:
        output_filename (str): Path for the output .tar.gz file.
        folder_names (list of str): List of folder names within RESULTS_DIR to include in the archive.
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_filename)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create the tar.gz file using tar and pigz for the specified folders within RESULTS_DIR
    with open(output_filename, 'wb') as f_out:
        # Switch to RESULTS_DIR and archive each specified folder by its name, avoiding full paths
        tar_command = ['tar', 'cf', '-', '-C', RESULTS_DIR] + folder_names
        p1 = subprocess.Popen(tar_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['pigz'], stdin=p1.stdout, stdout=f_out)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
    print(f'Tar file created: {output_filename}')


def compress_unified_results(sample_id, his_rot_sampling):
    """
    Compress specified input directories within RESULTS_DIR and save to the output directory.

    Parameters:
        sample_id (str): Identifier for the sample.
        his_rot_sampling (bool): Whether rotation sampling is enabled.
    """
    # Define input directories to include within RESULTS_DIR
    input_dirs = ['structures_with_predicted_zn', 'table']
    folder_suffix = "_".join(input_dirs)  # Join folder names for filename suffix

    # Define output directory
    output_dir = os.path.join(RESULTS_DIR, "compressed_results")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Format date safely for filenames
    date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    str_his_rot_sampling = "r-sampled" if his_rot_sampling else "r-unsampled"

    if sample_id:
        sample_id += "_"

    # Define output filename within the output directory
    output_filename = os.path.join(
        output_dir,
        f"ZincSight_{sample_id}{str_his_rot_sampling}_{folder_suffix}_{date}.tar.gz"
    )

    # Call the function to create the tar file with specified folder names (not full paths)
    create_results_tarfile(output_filename, input_dirs)


if __name__ == "__main__":
    sample_id = "input_id"
    his_rot_sampling_boolean = True
    compress_unified_results(sample_id, his_rot_sampling_boolean)