#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 09:58:57 2023

@author: mechti
"""

import math
import psycopg2
from src.settings import get_db_connection
import os
     
def final_scoring_and_insertion_to_table(list_query_structures_files_paths):
    def compute_and_rmsd_from_ideal_dis(list_dis,resi_comb_list_of_AA ):
        
      if list_dis== []: 
          return None
      number_cys_resi=resi_comb_list_of_AA.count('C')
   
     
      squared_deviations_cys = [(distance -2.32)**2 for distance in  list_dis[:number_cys_resi]]
      squared_deviations_noncys= [(distance -2.15)**2 for distance in  list_dis[number_cys_resi:]]
      all_squared_deviations=squared_deviations_cys+squared_deviations_noncys

          
      # squared_deviations = [(distance -ideal_dis)**2 for distance in list_dis]
      rmsd = math.sqrt(sum(all_squared_deviations) / len(all_squared_deviations))
      return round(rmsd, 4)

    conn = get_db_connection()
    cur1 = conn.cursor()
    
    # Create a copy of the table
    cur1.execute("CREATE TABLE IF NOT EXISTS scored_af_dataset_with_aggregated_final_tables AS SELECT * FROM af_dataset_with_metalcoord_test_sites_aggregated_final_tables;")
    conn.commit()
    
    # Alter the table to add a score column if it doesn't exist
    cur1.execute("ALTER TABLE scored_af_dataset_with_aggregated_final_tables ADD COLUMN IF NOT EXISTS score DOUBLE PRECISION;")
    conn.commit()
    
    cur1.execute("ALTER TABLE scored_af_dataset_with_aggregated_final_tables ADD COLUMN id SERIAL PRIMARY KEY;")
    conn.commit()
    
    cur1.execute("SELECT id, dif_angle_base, dif_angle_plane, distances_list, resi_comb, rmsd_overall,rmsd_close_atom FROM scored_af_dataset_with_aggregated_final_tables;")
    conn.commit()
 
    
    batch_size = 1000
    
    cur2 = conn.cursor()
    
    # Origin Hyperparameters
    a_1, a_2, b_1, b_2, b_3, b_4, z_1, z_2, z_3, y_1, y_2, y_3, y_4 = 0.3, 0.2, 12.3, 13.8, 13.8, 12.9, 5.5, 8.1, 14.4, 5.0, 0.7, 14.5, 1.6

    while True:
        rows = cur1.fetchmany(batch_size)
        if not rows:    # If no rows are returned, we've processed all rows
            break
        for row in rows:
 
                    resi_comb = row[4]
                    # rmsd_overall = row[5]
                    rmsd_close_atom = row[6]
                    if rmsd_close_atom is None:
                        continue
            
                    resi_comb_list_characters = list(resi_comb)
                    resi_comb_list_of_AA = [char for char in resi_comb_list_characters if char.isalpha()]
                    
                    angle_with_dis_score = None
                    row_id = row[0] # grabbing match's id from the selected rows
                    
                    if row[1] is not None and row[2] is not None and not math.isnan(row[1]) and not math.isnan(row[2]):
                        # Check dif_angle_base, dif_angle_plane values are not NaN (Means there is histidine)
                        angle_sum = row[1] + row[2] # Compute distance score
                    else:
                        angle_sum = None
            
                    if row[3] != []: # Compute distance score
                        dis_score = compute_and_rmsd_from_ideal_dis(row[3], resi_comb_list_of_AA)
                    else:
                        dis_score = None
            
                    if (angle_sum is not None) and (dis_score is not None):
                        if len(resi_comb_list_of_AA) > 3:
                            angle_with_dis_score = a_1 * angle_sum + b_1 * dis_score + y_1 * rmsd_close_atom
                        else:
                            angle_with_dis_score = a_2 * angle_sum + b_3 * dis_score + y_2 * rmsd_close_atom + z_2
            
                    if (angle_sum is None) and (dis_score is not None):
                        if len(resi_comb_list_of_AA) > 3:
                            angle_with_dis_score = b_2 * dis_score + y_3 * rmsd_close_atom + z_1  # need add -0.2
                        else:
                            angle_with_dis_score = b_4 * dis_score + y_4 * rmsd_close_atom + z_3
            
                    # Insert/Update the score in the database for the respective row
                    if angle_with_dis_score is not None:
                        cur2.execute(f"UPDATE scored_af_dataset_with_aggregated_final_tables SET score = {angle_with_dis_score} WHERE id = {row_id};")
            
        conn.commit()
        
    #Remove rows of structures which not belong to the queried structure dataset (af test dataset), but to the minimized represetnative structures
    list_query_structuresID=[]

    # Iterate over the files in the directory
    for structure_path in list_query_structures_files_paths:
        if os.path.isfile(structure_path):
            base_name = os.path.basename(structure_path)   # Get the base name of the file (filename with extension)
            structure_name_without_cif_extension, _ = os.path.splitext(base_name)            # Remove the file extension
            structure_ID = structure_name_without_cif_extension
            list_query_structuresID.append(structure_ID)

    query = "DELETE FROM scored_af_dataset_with_aggregated_final_tables WHERE pdbid_alphafoldmodel NOT IN %s;"
    cur1.execute(query, (tuple(list_query_structuresID),))
    
    conn.commit()
    conn.close()


