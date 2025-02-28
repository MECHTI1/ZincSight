#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 19:14:54 2023

@author: mechti
"""

"""
This script iterates through the final scored table to fetch, for each structure, the predicted metal coordinates (metalcoord) and scores.
It fetches data in batches of 1000 rows. For each batch, the script calls a function with a dictionary argument to create structure models with predicted Zn (zinc) sites.
The dictionary contains only structure IDs which have all their predicted zinc sites within the current batch.
This is achieved by omitting the last structure ID in each batch iteration, as there may be additional predicted zinc sites for this structure in the subsequent batch.
The omitted structure is then included in the next batch's dictionary, from which we gather data to create the structure models.
"""

import time
from src.create_structure_models_with_predicted_zn.create_pymol_session_structure_with_predicted_zn import create_pymol_session_structure_with_predicted_zn
import os

def locate_predicted_zn_within_structures(conn,list_query_structures_files_paths, path_output):
    list_tarred_sessions_paths=[]

    start=time.time()
    
    # Establish a connection to the PostgreSQL database
    
    # Create a cursor object
    cur= conn.cursor()
    cur.execute("SELECT structure_id, score, predicted_ion_pos, prob FROM final_compressed_table_with_scored_binding_sites ORDER BY structure_id ASC")
    conn.commit()
    
    batch_size = 1000
    
    # path_to_query_structure= "/home/mechti/Documents/PhD project/Git_algo/Not_in_git/Statistics_Alphafold_models/Choosing_accessions_for_AlphaFold_testset/Downloaded AlphaFold structures which contain after 2018 unique pfam and zinc binding/"
    
    dict_Key_structure_values_score_metalcoord={}
    last_StructureID_current_batch=None
    last_StructureID_previous_batch=None
    last_StructureID_batch_scores_metalcoords=None
    row_count=0
    AF_with_predcited_ZN=set()
    
    
    # Create a dictionary to map base names to full paths
    base_name_to_full_path = {}
    for file_path in list_query_structures_files_paths:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        base_name_to_full_path[base_name] = file_path
                      
    while True:
        rows = cur.fetchmany(batch_size)
        if last_StructureID_batch_scores_metalcoords!=None:
            dict_Key_structure_values_score_metalcoord[last_StructureID_previous_batch]=last_StructureID_batch_scores_metalcoords
       
        if not rows:    # If no rows are returned, we've processed all rows. So need to create structure of previous last structure ID. 
            break
       
        for row in rows:
            row_count += 1
            StructureID= row[0]
            score=row[1]
            metalcoord=row[2]
            prob = row[3]


            if StructureID not in dict_Key_structure_values_score_metalcoord:
                dict_Key_structure_values_score_metalcoord[StructureID]=[]
            dict_Key_structure_values_score_metalcoord[StructureID].append((score,prob,metalcoord))
        
        if row_count % 1000==0: # probably did not reached the end of the table (So there will be probably another batch after, which will inherent the last structure data from the current batch)
            last_StructureID_current_batch=StructureID
            last_StructureID_batch_scores_metalcoords=dict_Key_structure_values_score_metalcoord[StructureID]
            del dict_Key_structure_values_score_metalcoord[last_StructureID_current_batch]  # delete last StructureID from dict
            
        

        # Process the dictionary of structure IDs and their scores
        for structureID, scores_probs_metalcoords in dict_Key_structure_values_score_metalcoord.items():
            if structureID in base_name_to_full_path:
                full_path_to_structure = base_name_to_full_path[structureID]
                tarred_session_path= create_pymol_session_structure_with_predicted_zn(full_path_to_structure, scores_probs_metalcoords, path_output)
                list_tarred_sessions_paths.append(tarred_session_path)
                AF_with_predcited_ZN.add(structureID)

           
        dict_Key_structure_values_score_metalcoord= {} #reset dict for the next batch iteration
        last_StructureID_previous_batch=last_StructureID_current_batch
   
    
   # create pymol file if not predcited Zn also
    list_query_structuresIDs=[]
    for structure_path in list_query_structures_files_paths:
        if os.path.isfile(structure_path):
            base_name = os.path.basename(structure_path)   # Get the base name of the file (filename with extension)
            structure_name_without_cif_extension, _ = os.path.splitext(base_name)            # Remove the file extension
            structure_ID = structure_name_without_cif_extension
            list_query_structuresIDs.append(structure_ID)
            
    # check which AF not predicted to have zinc so need to create pymol file.
    AF_structures_without_predicted_sites = [structure for structure in list_query_structuresIDs if structure not in AF_with_predcited_ZN]
    print(AF_structures_without_predicted_sites)
    
    for structureID in AF_structures_without_predicted_sites:
        scores_probs_metalcoords=None
        full_path_to_structure=base_name_to_full_path[structureID]
        tarred_session_path= create_pymol_session_structure_with_predicted_zn(full_path_to_structure,scores_probs_metalcoords, path_output)
        list_tarred_sessions_paths.append(tarred_session_path)
    end=time.time()
    print ("time: ", end-start)