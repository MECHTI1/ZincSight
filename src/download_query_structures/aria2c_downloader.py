import subprocess
import time
import aria2p
from dotenv import load_dotenv
from pathlib import Path
from src.settings import QUERY_STRUCTURES_DIR
import os
import shutil
import gzip
import tarfile
import json
import concurrent.futures


# Mount Google Drive or Google Cloud Storage
def mount_storage(service):
    # Define the common mount point for GCS
    MOUNT_POINT = '/content/mounted_storage'  # Clearer and more descriptive

    if service.lower() == 'drive':
        from google.colab import drive
        MOUNT_POINT = '/content/drive'  # Google Drive uses the default /content/drive in Colab
        print(f"Mounting Google Drive at {MOUNT_POINT}...")
        drive.mount(MOUNT_POINT)
        print(f"Google Drive mounted successfully at {MOUNT_POINT}.")
        return MOUNT_POINT  # Return the mount point for Google Drive

    elif service.lower() == 'gcs':
        # Create the directory for mounting GCS
        if not os.path.exists(MOUNT_POINT):
            os.makedirs(MOUNT_POINT)

        # Install gcsfuse for mounting GCS
        print("Installing gcsfuse for Google Cloud Storage...")
        os.system('apt-get install -y gcsfuse')

        from google.colab import auth
        print("Authenticating Google Cloud...")
        auth.authenticate_user()

        # Set the project ID
        project_id = input("Enter your Google Cloud project ID: ")
        os.system(f'gcloud config set project {project_id}')

        # Ask for the bucket name and mount it
        bucket_name = input("Enter your Google Cloud Storage bucket name: ")
        print(f"Mounting Google Cloud Storage bucket '{bucket_name}' at {MOUNT_POINT}...")
        os.system(f'gcsfuse {bucket_name} {MOUNT_POINT}')
        print(f"Google Cloud Storage bucket '{bucket_name}' mounted successfully at {MOUNT_POINT}.")
        return MOUNT_POINT  # Return the mount point for GCS

    else:
        print("Unknown service. Please choose 'drive' or 'gcs'.")
        return None  # In case of unknown service, return None


def tar_or_tar_gz_uncompress_untar(complete_tar_gz_path, completed_downloads_dir, boolean_delete_after_extract=False):
    """
    Extract files from a .tar.gz archive and log the extracted files in a JSON file to avoid re-extraction.

    Parameters:
    - tar_gz_path (str): The path to the. tar / .tar.gz file.
    - target_dir (str): The directory to extract the files to.

    Returns:
    - A list of file names that were extracted in this session.
    """
    # Get the parent directory of the tar file
    # parent_directory = os.path.dirname(os.path.dirname(complete_tar_gz_path))

    # Check the extension and open with the correct mode
    if complete_tar_gz_path.endswith(".tar.gz") or complete_tar_gz_path.endswith(".tgz"):
        mode = "r:gz"
    elif complete_tar_gz_path.endswith(".tar"):
        mode = "r"
    else:
        raise ValueError("Unsupported file format")

    try:
        with tarfile.open(complete_tar_gz_path, mode) as tar:
            # Extract all contents
            tar.extractall(path=completed_downloads_dir)
            print(f"The following tar ot tar.gz file have been untarred or untarred and uncompressed: "
                  f"{complete_tar_gz_path}")

            # If extraction was successful, remove the original file
            if boolean_delete_after_extract:
                os.remove(complete_tar_gz_path)
                print(f"Original archive {complete_tar_gz_path} has been deleted.")


    except Exception as e:
        print(f"An error occurred: {e}. The original tar/ tar.gz file will not be uncompress / deleted.")


def move_completed_downloads(temp_dir, target_dir, max_files=1000, boolean_delete_after_extract=False):
    """
       Moves completed downloads from the temp_dir to the target_dir.
       If a completed download is a .tar.gz file, it lists and extracts its contents to the target_dir.

       Parameters:
       - temp_dir (str): The temporary directory where downloads are stored.
       - target_dir (str): The directory to move completed files to.
       - max_files (int): The maximum number of files to list/extract from .tar.gz archives.
       """

    # create a folder within temp_dir of completed downloads
    completed_downloads_temp_dir = os.path.join(temp_dir,
                                                'completed_downloads')  # Define the path for the completed downloads folder
    os.makedirs(completed_downloads_temp_dir, exist_ok=True)  # Create the directory if it doesn't exist

    all_files = os.listdir(temp_dir)
    aria2_files = {f[:-6] for f in all_files if f.endswith('.aria2')}
    completed_files = [f for f in all_files if f not in aria2_files and not f.endswith(
        '.aria2')]  # A completed download does not have an aria2 file.
    completed_structure_files = [f for f in completed_files if f.endswith(('.pdb', '.mmcif'))]
    print (completed_structure_files)
    if completed_structure_files != []:
        def move_file(completed_structure_file):
            src_file = os.path.join(temp_dir, file)
            dest_file = os.path.join(completed_downloads_temp_dir, completed_structure_file)
            shutil.move(src_file, dest_file)
            print(f"Moved {file} to {completed_downloads_temp_dir}")

        # Use ThreadPoolExecutor to move files in parallel
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(move_file, completed_structure_files)

    def process_tar_gz_file(tar_gz_file_name):
        print(f"Processing archive: {tar_gz_file_name}")
        tar_gz_file_path = os.path.join(temp_dir, tar_gz_file_name)
        tar_or_tar_gz_uncompress_untar(tar_gz_file_path, completed_downloads_temp_dir, boolean_delete_after_extract)

    tar_gz_files_to_process = filter(lambda f: f.endswith(('.tar.gz', '.tar')),
                                     completed_files)  # Filter only the .tar and .tar.gz files
    print (tar_gz_files_to_process)
    # Use ThreadPoolExecutor to process files in parallel without a loop
    with concurrent.futures.ThreadPoolExecutor() as executor:
        executor.map(process_tar_gz_file, tar_gz_files_to_process)

    # Move all completed files to the completed_downloads_temp_dir

    list_overall_files_moved_to_target_directory = []
    for file in completed_downloads_temp_dir:
        source_file_path = os.path.join(completed_downloads_temp_dir, file)
        target_file_path = os.path.join(target_dir, file)
        print ("source_file_path:", source_file_path,"target_file_path:",target_file_path)
        shutil.move(source_file_path, target_file_path)
        print(f"Moved file: {file} to {target_file_path}")

        # Delete the original non .tar.gz file
        os.remove(completed_downloads_temp_dir)
        print(f"Deleted original file: {completed_downloads_temp_dir}")

        list_overall_files_moved_to_target_directory.extend(file)

        if len(list_overall_files_moved_to_target_directory) >= max_files:
            return list_overall_files_moved_to_target_directory


def get_aria2():
    return aria2p.API(aria2p.Client(host="http://localhost", port=6800))  # No secret needed in Colab


# Function to read URLs from a file
def read_urls(file_path):
    try:
        with open(file_path, 'r') as file:
            return [line.strip() for line in file if line.strip()]
    except (FileNotFoundError, Exception) as e:
        print(f"Error reading file {file_path}: {e}")
        exit(1)


# Function to download files, skipping already completed ones
def download_files(urls, aria2,
                   temp_download_directory):  # first time think? or retrn to downloads. not start what completed already.
    gids = []
    for url in urls:
        file_name = url.split('/')[-1]  # Extract file name from URL
        file_path = os.path.join(temp_download_directory, file_name)

        if os.path.exists(file_path + '.aria2'):  #Add only if aria2 file exists- (Means that it started downloading...)
            try:
                download = aria2.add_uris([url])
                gids.append(download.gid)
                print(f"Continue download: {url} (GID: {download.gid})")
            except Exception as e:
                print(f'Error Continue downloading {url}: {e}')
        else:
            try:
                download = aria2.add_uris([url])
                gids.append(download.gid)  # Store GID to track download status later
                print(f"Started download: {url} (GID: {download.gid})")
            except Exception as e:
                print(f'Error downloading {url}: {e}')
    return gids


# Function to check the status of downloads with timeout and error handling
def monitor_downloads(aria2, gids, min_size_gb, timeout=2 * 60 * 60):
    min_size_bytes = min_size_gb * 1024 * 1024 * 1024  # Convert 10GB to bytes
    start_time = time.time()
    total_completed_downloads_size = 0
    while True:
        print("total_completed_downloads_size in bytes: ", total_completed_downloads_size)
        boolean_all_downloads_completed = True
        boolean_error_occurred = False
        for gid in gids:
            download = aria2.get_download(gid)
            print(f"File: {download.name} | Progress: {download.progress:.2f}% | "
                  f"Downloaded: {download.completed_length / (1024 * 1024):.2f} MB | "
                  f"Total size: {download.total_length / (1024 * 1024):.2f} MB | "
                  f"Speed: {download.download_speed / (1024 * 1024):.2f} MB/s | "
                  f"Status: {download.status} | Download {gid} status")
            if download.status != 'complete':  # is ['active', 'waiting', 'paused']:
                boolean_all_downloads_completed = False
            elif download.status == 'complete':
                total_completed_downloads_size += download.total_length
            elif download.status == 'error':
                boolean_error_occurred = True

        if boolean_all_downloads_completed:
            print("All downloads are complete.")
            break
        elif boolean_error_occurred:
            print(" an error downloading occurred.")
            break
        elif time.time() - start_time > timeout:
            print("Timeout reached. Exiting download monitor.")
            break
        # Check if the conditions are met (10 GB total downloaded and one file completed)
        elif total_completed_downloads_size >= min_size_bytes:
            print("10 GB downloaded and at least one file has completed downloading.")
            aria2.pause_all()  # Pause all downloads
            break
        time.sleep(10)

    return boolean_all_downloads_completed


def main(directory_structures_urls_to_download, temp_download_directory, target_directory=QUERY_STRUCTURES_DIR,
         min_size_gb=10, delete_after_extract=False):
    try:
        process = subprocess.Popen(['aria2c', '--enable-rpc', '--rpc-listen-all=false', '--continue=true', '--split=4',
                                    '--max-concurrent-downloads=3',
                                    '--max-connection-per-server=2'])
        print(f"Aria2c process started with PID: {process.pid}")
        # Wait for the process to complete
        time.sleep(5)  # Give it some time to initialize the RPC interface

        # time.sleep(10)  # Give it a moment to start up
    except Exception as e:
        print(f"Error starting aria2c: {e}")
        exit(1)

    # Check if temp download directory is empty
    if not os.path.exists(temp_download_directory):
        print("Temp download directory does not exists. Starting fresh downloads.")
    else:
        print("Temp download directory is exists. resume downloads.")

    urls = read_urls(directory_structures_urls_to_download)  # Read URLs
    aria2 = get_aria2()  # Get aria
    aria2.set_global_options({"dir": temp_download_directory})
    gids = download_files(urls, aria2, temp_download_directory)  # Start downloads and get GIDs
    returned_boolean_all_downloads_completed = monitor_downloads(aria2, gids,
                                                                 min_size_gb)  # Monitor downloads until completion
    move_completed_downloads(temp_download_directory, target_directory)  # Move completed downloads to target directory

    # Shutdown aria2c gracefully after downloads complete
    if process:
        try:
            process.terminate()
            process.wait()
            print("aria2c has been shut down.")
        except Exception as e:
            print(f"Error shutting down aria2c: {e}")

    return returned_boolean_all_downloads_completed


if __name__ == "__main__":
    target_directory = QUERY_STRUCTURES_DIR
    directory_structures_urls_to_download = "structures_urls_to_download.txt"

    load_dotenv()  # Load environment variables from the .env file

    if os.getenv('USER') == "mechti":
        temp_download_directory = str(os.getenv('LOCAL_STORAGE_DIRECTORY'))

    else:
        service = input("Enter 'drive' for Google Drive or 'gcs' for Google Cloud Storage: ")
        MOUNT_POINT = mount_storage(service)  # Get the mounted directory from the mount_storage function
        if MOUNT_POINT:
            temp_download_directory = os.path.join(MOUNT_POINT,
                                                   'downloads_directory')  # Set the download directory based on the mounted directory
            Path(temp_download_directory).mkdir(parents=True,
                                                exist_ok=True)  # Create the download directory if it doesn't exist
        else:
            print("No storage was mounted. Exiting.")
            exit(1)

    main(directory_structures_urls_to_download, temp_download_directory)
