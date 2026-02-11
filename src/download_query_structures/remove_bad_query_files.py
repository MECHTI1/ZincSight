import os
import subprocess

def remove_corrupted_cifs_grep(directory, required_tag="_atom_site.id", bad_list_name="bad_cif_files.txt"):
    """
    Fast bash-style check:
      find *.cif | xargs grep -L required_tag  -> list bad files
      then rm those files
    """
    directory = os.path.abspath(directory)  # string
    parent_dir = os.path.dirname(directory)  # string
    bad_list_path = os.path.join(parent_dir, bad_list_name)

    # Create list of bad CIFs (missing required_tag)
    cmd_find_grep = (
        f'find "{directory}" -name "*.cif" -print0 | '
        f'xargs -0 grep -L "{required_tag}" > "{bad_list_path}"'
    )
    subprocess.run(cmd_find_grep, shell=True, check=False)

    # If list non-empty -> remove
    if os.path.exists(bad_list_path) and os.path.getsize(bad_list_path) > 0:
        count = int(subprocess.check_output(f'wc -l < "{bad_list_path}"', shell=True).decode().strip())
        subprocess.run(f'xargs -a "{bad_list_path}" rm -f', shell=True, check=False)
        print(f"Removed {count} corrupted CIF files. List saved to: {bad_list_path}")
    else:
        print(f"No CIF files missing '{required_tag}'.")
