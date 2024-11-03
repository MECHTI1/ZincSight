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
from Bio.PDB import PDBParser, MMCIFIO
from src.settings import QUERY_STRUCTURES_DIR
from src.primary_script import main as execute



def str_clean_parse_tolist(input_string):
    cleaned_string = re.sub(r"[^\w,-]", "",input_string)  # Remove everything that isn't a letter, number, hyphen, or comma
    clean_items_list = [item.strip() for item in cleaned_string.split(",") if item.strip()]  # Split by commas, strip leading/trailing spaces, and remove any empty strings
    return clean_items_list

# TODO: need to check if I already implemented 'ESM' and 'PDB' options
def split_struct_db_sources(list_non_processed_input_names):
    af_list, pdb_list, esm_list = [], [], []
    for id_name in list_non_processed_input_names:
        if 'AF-' in id_name:
            af_list.append(id_name.split('-')[1])  # Append the af uniprot id to af_list
        elif id_name.startswith("MGY"):
            esm_list.append(id_name)  # Append full mgnify protein id to esm_list
        elif len(id_name) == 4 and id_name[0].isdigit(): #4 letter name with a number in the start considered a PDB structure
            pdb_list.append(id_name) # Insert to pdb_list
        else:  # If an input id name is not starting with 'MGY' or is a 4 letter name (with a number in the start), it will be considered as s uniprot af model
            af_list.append(id_name)
    return af_list, pdb_list, esm_list

def primary_download_structures_list_input(string_of_ids_to_download):
    #TODO: Check why wrote it-  primary_list_paths = []
    if string_of_ids_to_download:
        clean_id_list = str_clean_parse_tolist(string_of_ids_to_download)
        af_list, pdb_list, esm_list = split_struct_db_sources(clean_id_list)
        #print(af_list, pdb_list, esm_list)
        download_query_structures.main(af_list, pdb_list, esm_list)



def convert_all_pdb_to_cif_in_dir(directory):
    """Convert all PDB files in the provided directory and subdirectories to CIF format"""
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



def execute_zincsight(boolean_his_rot):
    try:
        def get_user_uploaded_ids_if_exists(parent_dir_path, file_name):
            file_path = os.path.join(parent_dir_path, file_name)
            # Check if folder exists, is not empty, and file exists with non-empty content
            if os.path.isfile(file_path):
                with open(file_path, 'r') as f:
                    content = f.read().strip()
                    return content if content else False
            return False

        str_struct_ids_for_download = get_user_uploaded_ids_if_exists('text_structures_ids_to_download', 'structure_ids_input.txt')
        if str_struct_ids_for_download:
            # Processes and downloads structure files based on input identifiers (AlphaFold/PDB/ESM formats)
            primary_download_structures_list_input(str_struct_ids_for_download)

        # convert PDB formatted query structures to mmCIF format
        convert_all_pdb_to_cif_in_dir(QUERY_STRUCTURES_DIR)

        list_query_structures_files_paths = []
        for root, dirs, files in os.walk(QUERY_STRUCTURES_DIR):
          for filename in files:
               list_query_structures_files_paths.append(os.path.join(QUERY_STRUCTURES_DIR, filename))
        execute(list_query_structures_files_paths, boolean_his_rot)

    except Exception as e:
        print(f"An exception occurred: {e}\n Please fill all the parameters.")
        sys.exit()

