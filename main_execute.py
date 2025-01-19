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
import re
from src.download_query_structures import download_query_structures
from Bio.PDB import PDBParser, MMCIFIO
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

def primary_download_structures_list_input(string_of_ids_to_download, path_query_structures):
    if string_of_ids_to_download:
        clean_id_list = str_clean_parse_tolist(string_of_ids_to_download)
        af_list, pdb_list, esm_list = split_struct_db_sources(clean_id_list)
        print (af_list, pdb_list, esm_list)
        download_query_structures.main(af_list, pdb_list, esm_list , path_query_structures)

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



def execute_zincsight(boolean_his_rot, structure_ids_for_download, path_query_structures, path_output):
    if structure_ids_for_download:
        # Processes and downloads structure files based on input identifiers (AlphaFold/PDB/ESM formats)
        primary_download_structures_list_input(structure_ids_for_download,path_query_structures)

    # convert PDB formatted query structures to mmCIF format
    convert_all_pdb_to_cif_in_dir(path_query_structures)

    # List all entries in the directory
    list_query_structures_files_paths = []
    for root, dirs, files in os.walk(path_query_structures):
        for filename in files:
            print(os.path.join(path_query_structures, filename))
            list_query_structures_files_paths.append(os.path.join(path_query_structures, filename))

    compressed_results_path = execute(list_query_structures_files_paths, boolean_his_rot, path_output)
    return compressed_results_path

if __name__=="__main__": #Behave like a test
    import time
    zincsight_start_execution = time.time()
    
    from src.settings import QUERY_STRUCTURES_DIR, RESULTS_DIR
    path_query_structures = QUERY_STRUCTURES_DIR
    path_output= RESULTS_DIR

    """
    Testing both Structure inputs: 
        1) Download structures:
            a) Manually written structure ids 
            b) Define path of txt file with structure ids 
        2) Upload structures
        """
    manually_written_structure_ids_for_download = ""
    structure_ids_from_txt_file= ""

    d_manual_ids, d_path_file_ids, u_structures =[True,False,False]

    if d_manual_ids or d_path_file_ids:
        if d_manual_ids:  # Testing - manually written structure ids
            manually_written_structure_ids_for_download = """A0A068N621, A0A0F6AZI6, A0A292DHH8, A0A2U3D0N8, A0A3F2YM30, A0A5H1ZR49, 
            G8ZFK7, P0A6G5, P38164,Q03760, Q08281, Q2K0Z2, Q2UFA9, Q5W0Q7, Q66K64, Q68EN5, Q6CXX6, Q7MVV4, 
            Q86T03, Q8N8R7, Q8NBJ9, Q9BWG6, Q9D1N4, Q9KP27, Q9M1V3, Q9NUN7, Q9NXF7"""
        if d_path_file_ids:  # Testing - defined path of txt file include structure ids
            # # Option 1 for input:
            # path_file_with_structure_ids = os.path.join("Query_structures_ids_txt_file","structures_ids_to_download.txt")

            # Option 2 for input:
            path_file_with_structure_ids = os.path.join("Query_structures_ids_txt_file", "2-D- UniprotID_of_RepClusters_nMem70Plus.csv")

            # Read the file and convert its content to a comma-separated string
            with open(path_file_with_structure_ids, "r") as file:
                lines = file.read().splitlines()  # Read all lines and remove newline characters
            structure_ids_from_txt_file = str(",".join(lines))   # Convert the list of lines to a comma-separated string


    # Combine the raw manual and file-based IDs (Create a combined string with these two strings)
    structure_ids_for_download = manually_written_structure_ids_for_download + structure_ids_from_txt_file


    boolean_his_rot = False
    compressed_results_path = execute_zincsight(boolean_his_rot,structure_ids_for_download, path_query_structures,path_output)
    print (compressed_results_path)
    
    zincsight_total_execution_time = time.time()-zincsight_start_execution
    print ("ZincSight_total_execution_time is: ", zincsight_total_execution_time, " Sec")

