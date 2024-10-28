import os
import subprocess
from src.settings import STRUCTURES_WITH_PREDICTED_ZN, RESULTS_DIR

def create_tarfile(output_filename, source_dir):
    # Ensure output file is not inside the source_dir to avoid recursion issues
    output_dir = os.path.dirname(output_filename)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create the tar.gz file using tar and pigz
    with open(output_filename, 'wb') as f_out:
        p1 = subprocess.Popen(['tar', 'cf', '-', '-C', source_dir, '.'], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['pigz'], stdin=p1.stdout, stdout=f_out)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
    print(f'Tar file created: {output_filename}')



def main():# Create a new directory within RESULTS_DIR
    compressed_dir = os.path.join(RESULTS_DIR, "compressed_structures_with_predicted_zn")
    os.makedirs(compressed_dir, exist_ok=True)
    
    # Define the output file path within the new directory
    output_filename = os.path.join(compressed_dir, "structures_with_predicted_zn.tar.gz")
    source_dir = STRUCTURES_WITH_PREDICTED_ZN
    
    create_tarfile(output_filename, source_dir)




if __name__ == "__main__":
    main()




