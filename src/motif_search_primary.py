#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 17:15:45 2023

@author: mechti
"""

import time
import psycopg2
from src.settings import get_db_connection
import psycopg2.extras
import numpy as np
import src.ii_search
from src.qcp_superimposer import calculate_rmsd
from collections import defaultdict

def Create_Coordinates_OF_Matches_Table(conn):
    
    cur = conn.cursor()
    
    #create table with coordinates will called rmsd 
    
    # Execute the query to fetch the first row
    cur.execute("SELECT * FROM AF_DATASET_final_motif_search_table_V2 LIMIT 1;")
    
    # Fetch the result
    first_row = cur.fetchone()
    
    # Check the length of the row to get the number of columns
    num_residues_in_motif = len(first_row)-2
    print ("num_residues_in_motif", num_residues_in_motif)
    base_query = """
    SELECT 
        f.match_id, 
        {resi_order} as resi_order, 
        f.chain_resi_{resi_order} as chain_resi, 
        c.close_atom_coord, 
        c.far_atom_coord
    FROM 
        AF_DATASET_final_motif_search_table_V2 f 
    INNER JOIN 
        AF_DATASET_coordinates_table_v2 c ON f.pdbid_alphafoldmodel = c.pdbid_alphafoldmodel AND f.chain_resi_{resi_order} = c.chain_resi
    """
    
    queries = [base_query.format(resi_order=i) for i in range(1, num_residues_in_motif + 1)]
    final_query = "CREATE TABLE IF NOT EXISTS AF_DATASET_Coordinates_of_Matches_Table_V2 AS " + " UNION".join(queries) + "ORDER BY match_id, resi_order;"
    
    print(final_query)
    
    cur.execute(final_query) # Execute the query
    conn.commit()
    


def create_empty_RMSD_of_matches_table(conn):
    cur = conn.cursor()
  
    # SQL command to create Primary_metal_binding_site_table
    MOTIF_SEARCH_TABLE_TEMP_creation_command = """
    CREATE TABLE IF NOT EXISTS AF_DATASET_rmsd_of_matches_table_V2(
        match_id INTEGER REFERENCES AF_DATASET_final_motif_search_table_V2(match_id),
        rmsd_overall FLOAT,
        rmsd_close_atom FLOAT,
        calculated_weighted_RMSD FLOAT
     );
     """    
    cur.execute (MOTIF_SEARCH_TABLE_TEMP_creation_command)
    


def Create_RMSD_of_matches_table(conn):
    cur1 = conn.cursor()
    cur2= conn.cursor()
    create_empty_RMSD_of_matches_table(conn)
    
    # select coordinates of residues in the original metal binding motif that running the search with.
    cur1.execute("""
        SELECT *
        FROM AF_DATASET_Coordinates_of_Matches_Table_V2
        WHERE match_id = 1
        ORDER BY resi_order
    """)
    # Fetch all the results
    rows_original_metal_binding_site = cur1.fetchall()
    num_residues_in_metal_binding_site= len (rows_original_metal_binding_site)  # get number of residues in metal binding site, thats in order to know how many rows to select in a batch (dont want to fetch partially matching sites)
    
    print (num_residues_in_metal_binding_site)
    
    # Extract the relevant part of the data- coords
    coords_list_original_metal_binding_site = [coords for _, _, _, close_coords, far_coords in rows_original_metal_binding_site for coords in (close_coords, far_coords)]

    # Convert to a numpy array
    coords_original_metal_binding_site = np.array(coords_list_original_metal_binding_site)
    print (coords_original_metal_binding_site)
    #batch of 1000 matched_ids each time (each match have multiple residues)- means- 1000 matched sites each time.
        # execute your query
    cur1.execute("SELECT * FROM AF_DATASET_Coordinates_of_Matches_Table_V2")
    
    # fetch rows in batches
    batch_size = 1000 *  num_residues_in_metal_binding_site
    
    start_rmsd_time = time.time()  # get current time
    while True:
        thousand_sites_coordinates_without_nested_lists_of_matches = cur1.fetchmany(batch_size)
        if not thousand_sites_coordinates_without_nested_lists_of_matches:
            break
        # Group every 4 tuples into a sublist

        # Initialize a defaultdict to hold lists of rows for each match_id
        matches_dict = defaultdict(list)
        
        # Iterate through fetched rows and append each row to the list corresponding to its match_id
        for row in thousand_sites_coordinates_without_nested_lists_of_matches:
            match_id = row[0]  # Assuming the match_id is the first element of the row
            matches_dict[match_id].append(row)
        
        # Convert the dictionary values to a list to get a list of grouped matches
        thousand_sites_coordinates_with_nested_lists_of_matches = list(matches_dict.values())            
        list_tuples_to_insert_to_RMSD_of_matches_table=[]
        for matched_motif in  thousand_sites_coordinates_with_nested_lists_of_matches:
           # print(matched_motif)
            # Extract the relevant part of the data- coords
            coords_list_matched_predicted_binding_site = [coords for _, _, _, close_coords, far_coords in matched_motif for coords in (close_coords, far_coords)]

            # Convert to a numpy array
            coords_matched_predicted_metal_binding_site = np.array( coords_list_matched_predicted_binding_site)
           # print (coords_matched_predicted_metal_binding_site)
            match_id= matched_motif [0][0]
            #print ( match_id)
            #print (num_residues_in_metal_binding_site)
            
            rmsd_overall,rmsd_close_atom= calculate_rmsd (coords_original_metal_binding_site, coords_matched_predicted_metal_binding_site)
            
            if (rmsd_overall,rmsd_close_atom)!= ("error", "error"):
                single_tuple_to_insert_to_RMSD_of_matches_table= (match_id, rmsd_overall, rmsd_close_atom)
                list_tuples_to_insert_to_RMSD_of_matches_table.append (single_tuple_to_insert_to_RMSD_of_matches_table)
            else:
                print ("rmsd error")
        # SQL insert command for METAL_BINDING_SITE_II_TABLE 
        insert_command_list_tuples_to_insert_to_RMSD_of_matches_table = "INSERT INTO AF_DATASET_rmsd_of_matches_table_V2 (match_id ,rmsd_overall ,rmsd_close_atom)  VALUES %s"
      #  print ( list_tuples_to_insert_to_RMSD_of_matches_table)
        psycopg2.extras.execute_values(cur2, insert_command_list_tuples_to_insert_to_RMSD_of_matches_table, list_tuples_to_insert_to_RMSD_of_matches_table)

    end_rmsd_time=time.time()
    total_rmsd_time =end_rmsd_time-start_rmsd_time
    print(f"For RMSD The code executed in {total_rmsd_time} seconds")

def Add_resi_comb_column(conn,site_id):
    cur = conn.cursor()
    cur.execute("SELECT resi_comb FROM minimized_training_cluster_information WHERE site_id = %s;", (site_id,))
    resi_comb= cur.fetchone()[0]
    
    cur.execute("ALTER TABLE AF_DATASET_final_motif_search_table_V2 ADD COLUMN resi_comb TEXT;")
    cur.execute("UPDATE AF_DATASET_final_motif_search_table_V2 SET resi_comb = %s;", (resi_comb,))
    
    conn.commit()


def main(site_id):
           
    start_time = time.time()  # get current time
           
    Boolean_whether_delete_all_motif_search_table_temp_after_iteration= True
    conn = get_db_connection()
    src.ii_search.main(site_id, Boolean_whether_delete_all_motif_search_table_temp_after_iteration)
    Create_Coordinates_OF_Matches_Table (conn)
    Create_RMSD_of_matches_table(conn)
    Add_resi_comb_column(conn,site_id)
    cur = conn.cursor()
    
    
    # Commit the transaction
    conn.commit()
    
    # Close the connection
    cur.close()
    conn.close()
    
    end_time = time.time()  # get current time after executing your code
    execution_time = end_time - start_time  # calculate the difference in time
    print(f"The code executed in {execution_time} seconds")


if __name__ == "__main__":
    site_id= 2#slow -74, 8.  511 result of 74 
    main(site_id)






