# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os

import psycopg2
from src.settings import get_db_connection
import json
from src.dis_angles_calculations.with_rotations.create_dict_resi_dis_hisangles_zncoord_rot180deg import caclulate_dis_CoordinationAngles_HISangles_stats
import time 
import multiprocessing
import numpy as np

def create_detailed_coordinates_of_matches_table(conn):
    
    
    # Create a cursor object
    cur = conn.cursor()
    
    # Create new table "detailed_coordinates_of_matches_table" if it doesn't exist
    cur.execute("""
        CREATE TABLE IF NOT EXISTS AF_DATASET_test_detailed_coordinates_of_matches_table_V2 (
            match_id INT,
            pdbid_alphafoldmodel TEXT,
            chain_resi TEXT,
            resi_type TEXT,
            dict_atom_coord JSONB
        )
    """)
    
    # Fetch column names from the final_motif_search_table
    cur.execute("""
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'af_dataset_final_motif_search_table_v2'
    """)
    
    column_names = [row[0] for row in cur.fetchall()]
    # print ("column_names", column_names)
    # Filter out chain_resi columns
    chain_resi_columns = [name for name in column_names if name.startswith('chain_resi')]
    # Create the unnest array argument string
    unnest_arg_str = ', '.join(chain_resi_columns)
    # print (unnest_arg_str)
    # Create temporary table dynamically based on the number of chain_resi columns
    cur.execute(f"""
        CREATE TEMPORARY TABLE temp_final_motif AS
        SELECT match_id, pdbid_alphafoldmodel, UNNEST(ARRAY[{unnest_arg_str}]) AS chain_resi
        FROM AF_DATASET_final_motif_search_table_v2
    """)

    # Insert data into new table from "final_motif_search_table" and "detailed_coordinates_table_v2"
    cur.execute("""
        INSERT INTO AF_DATASET_test_detailed_coordinates_of_matches_table_V2 (match_id, pdbid_alphafoldmodel, chain_resi, resi_type, dict_atom_coord)
        SELECT 
            t.match_id, 
            t.pdbid_alphafoldmodel, 
            t.chain_resi, 
            a.resi_type,
            a.dict_atom_coord
        FROM 
            temp_final_motif t
        JOIN 
            AF_DATASET_detailed_coordinates_table_v2 a ON t.pdbid_alphafoldmodel = a.pdbid_alphafoldmodel AND t.chain_resi = a.chain_resi
    """)
    
    # Commit changes
    conn.commit()


def add_columns_of_dist_angles_stats_to_final_motif_search_table(conn):
    cur = conn.cursor()
    cur.execute("""
        ALTER TABLE AF_DATASET_final_motif_search_table_v2 
        ADD COLUMN distances_list float[],
        ADD COLUMN dif_angle_base float,
        ADD COLUMN dif_angle_plane float,
        ADD COLUMN Coordination_anlges float[],
        ADD COLUMN metalcoord float[]
    """)
    conn.commit()     # Commit changes


def add_dist_angle_stats(conn, dict_stats_batch):
    cur = conn.cursor()
    for match_id, values in dict_stats_batch.items():
        if dict_stats_batch[match_id] == "error":
            print("error")
            print(dict_stats_batch[match_id])
            print(match_id)
            continue

        # Convert NumPy arrays/lists to native Python lists of floats
        def safe_convert_array(arr):
            if isinstance(arr, (list, np.ndarray)):
                return [float(x) for x in arr]
            return None

        def safe_convert_scalar(val):
            if isinstance(val, np.generic):
                return val.item()
            return float(val) if isinstance(val, (float, int)) else None


        distances_list = safe_convert_array(values.get('distances_list', []))
        coordination_angles = safe_convert_array(values.get('Coordination_anlges', []))
        metalcoord = safe_convert_array(values.get('metalcoord', []))
        dif_angle_base = safe_convert_scalar(values['HIS_anlges_stats']['candidate_point_angles_stast'].get('dif_angle_base'))
        dif_angle_plane = safe_convert_scalar(values['HIS_anlges_stats']['candidate_point_angles_stast'].get('dif_angle_plane'))

        cur.execute("""
            UPDATE AF_DATASET_final_motif_search_table_v2 
            SET 
                distances_list = %s,
                dif_angle_base = %s,
                dif_angle_plane = %s,
                Coordination_anlges = %s,
                metalcoord= %s
            WHERE match_id = %s
        """, (distances_list, dif_angle_base, dif_angle_plane, coordination_angles, metalcoord, match_id))
    conn.commit()


def worker(item):
    match_id, dict_residues_in_candidate_motif = item
    result = caclulate_dis_CoordinationAngles_HISangles_stats(dict_residues_in_candidate_motif)
    return match_id, result

def fetch_and_process_data(conn,num_cores):
    # Create a cursor object
    cur = conn.cursor()
    
    start_time = time.time()

    # Fetch unique match_ids
    cur.execute("SELECT DISTINCT match_id FROM AF_DATASET_test_detailed_coordinates_of_matches_table_V2")
    unique_match_ids = [row[0] for row in cur.fetchall()]
    
    
    # Fetch match_ids in batches
    unique_match_ids_batch_size = 1000
    for i in range(0, len(unique_match_ids), unique_match_ids_batch_size):
        data_dict = {}

        # print (len(unique_match_ids))
        batch_match_ids = unique_match_ids[i:i + unique_match_ids_batch_size]
        
        # Convert list to string to be used in SQL command
        batch_match_ids_str = ','.join(map(str, batch_match_ids))
       # print (batch_match_ids_str)

        # Fetch rows for match_ids in current batch
        cur.execute(f"""
            SELECT match_id, chain_resi, resi_type, dict_atom_coord
            FROM AF_DATASET_test_detailed_coordinates_of_matches_table_V2
            WHERE match_id IN ({batch_match_ids_str})
        """)
        rows = cur.fetchall()
        
        # Process rows and store in dictionary
        for row in rows:
            match_id, chain_resi, resi_type, dict_atom_coord = row
            if match_id in data_dict:
                data_dict[match_id][chain_resi] = [resi_type, dict_atom_coord]
            else:
                data_dict[match_id] = {chain_resi: [resi_type, dict_atom_coord]}

      # #  print(data_dict)
      #   data_dict.clear()


        # Set the number of processes to 2 explicitly
        with multiprocessing.Pool(processes=num_cores) as pool:
            results = pool.map(worker, data_dict.items())
        dict_match_id_to_dis_CoordinationAngles_HISangles_stats_single_batch = dict(results)

        # Add the results to the database or further processing
        add_dist_angle_stats(conn, dict_match_id_to_dis_CoordinationAngles_HISangles_stats_single_batch)

    end_time = time.time()
    print ("total_time",end_time-start_time )



def main(num_cores):
    # Establish a connection to the PostgreSQL database
    conn = get_db_connection()
    create_detailed_coordinates_of_matches_table(conn)
    add_columns_of_dist_angles_stats_to_final_motif_search_table(conn)
    fetch_and_process_data(conn,num_cores)
    
    boolean_whether_delete_detailed_coordinates_of_matches_table=True
    if boolean_whether_delete_detailed_coordinates_of_matches_table:
        cur = conn.cursor()
        cur.execute("DROP TABLE  AF_DATASET_test_detailed_coordinates_of_matches_table_V2")
        conn.commit()

         
    conn.close()
 
if __name__ == "__main__":
    main()
