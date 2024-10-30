#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 23:29:54 2024

@author: mechti
"""

"""
"main.py" mediates between user input (jupyter notebook) to the src/primary_scirpt.py
This script actions:

 (1)   
 * If the user upload structures files:
     ** '.pdb' files converted to '.mmcif formats'.
 * If the user write/ upload the structures ID list :
    ** Process the user input into a format readable by the 'src/download_query_structures' package.
    ** Download structures from web.
    ** '.pdb' files converted to '.mmcif formats'

 (2) 
    Execute the ZincSight prediction process.
    This execution run on the query structures (Which located within the 'Query AlphaFold structures' dir
"""

import os
import sys
import re
from src.download_query_structures import download_query_structures
from src.download_query_structures.aria2c_downloader import mount_storage
from src.download_query_structures import aria2c_downloader
from Bio.PDB import PDBParser, MMCIFIO
from src.settings import QUERY_STRUCTURES_DIR
from src.primary_script import main as execute
import requests
from dotenv import load_dotenv


def str_clean_parse_tolist(input_string):
    cleaned_string = re.sub(r"[^\w,-]", "",input_string)  # Remove everything that isn't a letter, number, hyphen, or comma
    clean_items_list = [item.strip() for item in cleaned_string.split(",") if item.strip()]  # Split by commas, strip leading/trailing spaces, and remove any empty strings
    return clean_items_list

def split_struct_db_sources(list_non_processed_input_names):
    af_list, pdb_list, esm_list = [], [], []
    for id_name in list_non_processed_input_names:
        if 'AF-' in id_name: # Insert AF uniprot ID into af_list
            af_list.append(id_name.split('-')[1])  # Taking only the uniprot ID if bring full name
        elif id_name.startswith("MGY"): # Insert to esm_list
            esm_list.append(id_name)
            print("need to add this 'ESM' option and 'PDB'")
        elif len(id_name) == 4 and id_name[0].isdigit(): # Insert to pdb_list
            pdb_list.append(id_name)
            print("need to add this 'ESM' option and 'PDB'")
        else:  # AF uniprot ID into af_list. If an item name not satisfied the above, function consider it as UniProt model of AF model
            af_list.append(id_name)
    return af_list, pdb_list, esm_list

def primary_download_structures_list_input(string_of_ids_to_download):
    #TODO: Check why wrote it-  primary_list_paths = []
    if string_of_ids_to_download:
        clean_id_list = str_clean_parse_tolist(string_of_ids_to_download)
        af_list, pdb_list, esm_list = split_struct_db_sources(clean_id_list)
        #print(af_list, pdb_list, esm_list)
        download_query_structures.main(af_list, pdb_list, esm_list)

# Convert all PDB files in the provided directory and subdirectories to CIF format.
def convert_all_pdb_in_dir(directory):
    # Helper function to convert PDB to CIF
    def convert_pdb_to_cif(pdb_file_path):
        try:
            parser = PDBParser(QUIET=True)  # Initialize PDBParser
            structure_id = os.path.basename(pdb_file_path).replace('.pdb', '')  # Get a unique structure ID
            structure = parser.get_structure(structure_id, pdb_file_path)  # Parse the structure from the .pdb file
            io = MMCIFIO()  # Initialize MMCIFIO to save the structure in CIF format
            io.set_structure(structure)
            cif_file_path = pdb_file_path.replace('.pdb', '.cif')  # Save as .cif file (replacing .pdb with .cif)
            io.save(cif_file_path)
            os.remove(pdb_file_path)  # Delete the original .pdb file
        except Exception as e:
            print(f"Failed: {pdb_file_path}converting {pdb_file_path}: {str(e)}")
        print(f"Converted {pdb_file_path} to {cif_file_path}")

    # Walk through directory and subdirectories
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".pdb"):
                pdb_file_path = os.path.join(root, file)
                convert_pdb_to_cif(pdb_file_path)


if __name__ == "__main__":
    #Checking_algo_with_default_input = True
    try:
        #option_downloading = None

        print("For downloading structures press: 'd'. To upload structures: 'u'. For both options press '2'")
        pressed_choice = input()

        if pressed_choice == "2":
            option_downloading = True
            option_uploading = True
        elif pressed_choice == "u":
            option_downloading = False
            option_uploading = True
        elif pressed_choice == "d":
            option_downloading = True
            option_uploading = False

        if option_downloading:

            print(""" Option 1: 
                        Please write a list of the protein structure identifiers- 'PDB','AF model' or 'ESMFold models'.
                        For AF models- write a list of UniProt identifiers.
                        For ESMFold models- write a list of their MGnify identifiers                       
                      Option 2:
                        Upload a txt file list of structures download urls within local
                      Press 1 or 2 according your preference""")

            preference_list1_or_linkstxt2 = input()

            if int(preference_list1_or_linkstxt2) == 1:
                if Checking_algo_with_default_input:
                    string_of_ids_to_download = """A0A068N621, A0A0F6AZI6, A0A292DHH8, A0A2U3D0N8, A0A3F2YM30,
                                               A0A5H1ZR49, G8ZFK7, O60232, P0A6G5, P0DUH5, P37659, P38164,
                                               Q03760, Q08281, Q2K0Z2, Q2UFA9, Q5W0Q7, Q66K64, Q68EN5,
                                               Q6CXX6, Q7MVV4, Q86T03, Q8N8R7, Q8NBJ9, Q96JC1, Q9BWG6,
                                               Q9D1N4, Q9KP27, Q9M1V3, Q9NUN7, Q9NXF7"""
                else:
                    string_of_ids_to_download = input()

                primary_download_structures_list_input(string_of_ids_to_download)

            elif int(preference_list1_or_linkstxt2) == 2:

                # Add the urls txt file
                if Checking_algo_with_default_input:
                    # Testing_with_default_input
                    url = "https://raw.githubusercontent.com/facebookresearch/esm/main/scripts/atlas/v0/full/tarballs.txt"  # GitHub raw file URL
                    response = requests.get(url)
                    directory_structures_urls_to_download = "src/download_query_structures/structures_urls_to_download.txt"
                    with open(directory_structures_urls_to_download, "w") as file:
                        file.write("\n".join(response.text.splitlines()[
                                             :3]))  # Not using all links since then time check process would be time consuming

                else:
                    pass
                    # TODO: add the colab box to upload the txt file

                # Download structures to target directory and check if finish downloads with aria2c status to continue after
                load_dotenv()  # Load environment variables from the .env file
                if os.getenv('USER') == "mechti":
                    temp_download_directory = str(os.getenv('LOCAL_STORAGE_DIRECTORY'))
                    print(f"Using local storage directory: {temp_download_directory}")
                else:
                    print("Please choose a storage method for mounting: Google Drive or GCS.")
                    service = input()
                    # Get the mounted directory from the mount_storage function
                    MOUNT_POINT = mount_storage(service)

                    if MOUNT_POINT:
                        temp_download_directory = os.path.join(MOUNT_POINT,
                                                               'downloads_directory')  # Set the download directory based on the mounted directory
                        os.makedirs(temp_download_directory,
                                    exist_ok=True)  # Create the download directory if it doesn't exist
                    else:
                        print("No storage was mounted. Exiting.")
                        exit(1)
                target_directory = QUERY_STRUCTURES_DIR

                aria2c_status = aria2c_downloader.main(directory_structures_urls_to_download, temp_download_directory)

        if option_uploading:  # need to add interface uploading, when browse and select structure files
            pass
            # TODO: add the colab box to upload the structures and need to mount drive or google cloud for tracking process

        convert_all_pdb_in_dir(QUERY_STRUCTURES_DIR)

        list_query_structures_files_paths = []
        for root, dirs, files in os.walk(QUERY_STRUCTURES_DIR):
            for filename in files:
                print(os.path.join(QUERY_STRUCTURES_DIR, filename))
                list_query_structures_files_paths.append(os.path.join(QUERY_STRUCTURES_DIR, filename))

    except Exception as e:
        print(f"An exception occurred: {e}\nPlease fill all the parameters.")
        sys.exit()

    boolean_rotamer_examination = False  # initialize variable to False
    while True:
        print("Press N/Y for rotamer examination: ")
        user_input = input().upper()
        if user_input == "Y":
            boolean_rotamer_examination = True
            break
        elif user_input == "N":
            boolean_rotamer_examination = False
            break
        else:
            print("Invalid input. Please press N or Y.")

    execute(list_query_structures_files_paths, boolean_rotamer_examination)
