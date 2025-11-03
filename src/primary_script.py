#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 18:31:24 2023
@author: mechti
"""
import os
from src.settings import get_db_connection, KEEP_TEMP_TABLES, DEBUGGING, cleanup_tables
import time
from src.create_ii_coordinates_tables.ii_coordinates_primary_generator_structures import main as create_ii_coordinates_tables_query_dataset
from src.create_ii_coordinates_tables.Insert_Representative_Motifs_to_Query_InitialTables import main as insert_representative_motifs_to_dataset_tables
from src.motif_search_primary import main as main_first_step
from src.assign_lowest_rmsd_resi_comb import main as main_second_step
from src.dis_angles_calculations.without_rotations.dis_angles_calc_primary_generator_without_his_rotation import main as main_third_step_without_his_rotation
from src.dis_angles_calculations.with_rotations.dis_angles_calc_primary_generator_with_his_rotation_180deg import main as main_third_step_his_rotation_180deg
from src.scoring_and_compression.add_scores_to_table import final_scoring_and_insertion_to_table  
from src.scoring_and_compression.compress_table_by_proximity import table_compression
from src.refine_results_table import refine_table
from src.create_structure_models_with_predicted_zn.primary_create_structure_models_with_predicted_zn import locate_predicted_zn_within_structures
from src.export_final_table_to_csv_format import export_final_table_to_csv_file
from src.compress_results import compress_unified_results
from src.add_prob.add_prob_to_final_table import add_column_with_probs
from src.db_debugging import debug_print_last_table

def main(list_query_structures_files_paths, boolean_rotamer_examination, path_output,num_cores):
    start_time_prediction = time.time()
    conn = get_db_connection()
    cur = conn.cursor()

    cleanup_tables(cur, conn)
    his_rotation = boolean_rotamer_examination

    start_time_create_tables=time.time()

    create_ii_coordinates_tables_query_dataset(list_query_structures_files_paths,num_cores)
    insert_representative_motifs_to_dataset_tables()

    end_time_create_tables=time.time()
    time_create_tables=end_time_create_tables-start_time_create_tables

    # all_sites_final_table
    cur.execute("""
          CREATE TABLE AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables (
            match_id int,
            pdbid_alphafoldmodel varchar(255),
            chain_resi_1 varchar(255),
            chain_resi_2 varchar(255),
            chain_resi_3 varchar(255),
            chain_resi_4 varchar(255),
            chain_resi_5 varchar(255),
            resi_comb text,
            rmsd_overall float,
            rmsd_close_atom float,
            site_id int,
            distances_list float[],
            dif_angle_base float,
            dif_angle_plane float,
            coordination_anlges float[],
            metalcoord float[]
            )
          """)
    conn.commit()
    print ("mile 1")

    # Execute a SELECT statement to get all the site_id values
    cur.execute("SELECT site_id FROM minimized_training_cluster_information;")
    results = cur.fetchall() # Fetch all the results
    # Extract the site_id values into a list
    site_ids = [result[0] for result in results]
    for site_id in site_ids:
        main_first_step(site_id)
        main_second_step()
        
        if his_rotation:
            main_third_step_his_rotation_180deg(num_cores)
        else:
            main_third_step_without_his_rotation()


        cur = conn.cursor()

        # Define the name of your source table
        source_table = 'af_dataset_final_motif_search_table_v2'
        
        # Query the columns in the source table
        cur.execute(f"SELECT column_name FROM information_schema.columns WHERE table_name = '{source_table}';")
        columns = [row[0] for row in cur.fetchall()]
        # Define the columns in the destination table
        destination_columns = ['match_id', 'pdbid_alphafoldmodel', 'chain_resi_1', 'chain_resi_2', 'chain_resi_3', 'chain_resi_4', 'chain_resi_5', 'resi_comb', 'rmsd_overall', 'rmsd_close_atom', 'site_id','distances_list', 'dif_angle_base', 'dif_angle_plane', 'coordination_anlges','metalcoord']
    
         
        # Build the SELECT statement dynamically
        select_clause = []
        for column in destination_columns:
            if column in columns:
                if '_list' in column or 'angles' in column:  # Assuming only the columns with '_list' and 'angles' in their names are of type array
                    select_clause.append(f"COALESCE({column}, ARRAY[]::double precision[])")
                else:
                    select_clause.append(f"COALESCE({column}, NULL)")
            else:
                if '_list' in column or 'angles' in column:  # Assuming only the columns with '_list' and 'angles' in their names are of type array
                    select_clause.append("ARRAY[]::double precision[]")
                else:
                    select_clause.append("NULL")
        
        select_clause = ", ".join(select_clause)
        
        #Insertion of values to aggregated table
        sql_query = f"""
            INSERT INTO AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables ({', '.join(destination_columns)})
            SELECT
                {select_clause}
            FROM
                {source_table};
        """
        cur.execute(sql_query)
        conn.commit()
        cur.execute("DROP TABLE AF_DATASET_final_motif_search_table_v2")
        conn.commit()

    if not KEEP_TEMP_TABLES:
            cur.execute("DROP TABLE af_dataset_coordinates_table_v2")
            cur.execute("DROP TABLE af_dataset_detailed_coordinates_table_v2")
            cur.execute("DROP TABLE af_dataset_ii_table_v2")
    conn.commit()
    print ("KEEP_TEMP_TABLES=", KEEP_TEMP_TABLES)

    if KEEP_TEMP_TABLES:
        cur.execute("CREATE TABLE dropped_with_erroredrows AS SELECT * FROM AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables WHERE dif_angle_base IS NULL OR rmsd_overall IS NULL;")   #  creates new table with the dropped rows, for scenarios of debugging.
    cur.execute("DELETE FROM AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables WHERE dif_angle_base IS NULL OR rmsd_overall IS NULL;")
    conn.commit()

    final_scoring_and_insertion_to_table(list_query_structures_files_paths)
    if not KEEP_TEMP_TABLES:
        cur.execute("DROP TABLE AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables")
    conn.commit()

    table_compression(conn)
    conn.commit()

    boolean_if_any_predicted_site = add_column_with_probs(conn)
    if not boolean_if_any_predicted_site: return False # no predicted sites
    refine_table()
    conn.commit()

    if not KEEP_TEMP_TABLES:
        cur.execute("DROP TABLE scored_af_dataset_with_aggregated_final_tables")
    conn.commit()

    start_time_create_models= time.time()
    locate_predicted_zn_within_structures(conn,list_query_structures_files_paths, path_output)
    end_time_create_models = time.time()
    

    export_final_table_to_csv_file(path_output)
    compressed_results_path = compress_unified_results('sample_id', his_rotation, path_output)    #TODO: Add option of input sample_id

    if DEBUGGING == True:
        debug_print_last_table(conn,'final_compressed_table_with_scored_binding_sites')

    if not KEEP_TEMP_TABLES:
            cur.execute("DROP TABLE final_compressed_table_with_scored_binding_sites")
    print ("time_create_II_Coordinates_tables",time_create_tables)
    print ("total time for predictions, not including downloading structures and compression: ", start_time_create_models- start_time_prediction)
    print ("total time for create structures: ", end_time_create_models- start_time_create_models)
    conn.commit()
    conn.close()

    return compressed_results_path

if __name__=="__main__":
    from src.settings import QUERY_STRUCTURES_DIR, RESULTS_DIR

    list_of_uniprot_accessions = ["A0A068N621", "A0A0F6AZI6", "A0A292DHH8", "A0A2U3D0N8", "A0A3F2YM30", "A0A5H1ZR49", "G8ZFK7", "O60232", "P0A6G5", "P0DUH5", "P37659", "P38164", "Q03760", "Q08281", "Q2K0Z2", "Q2UFA9", "Q5W0Q7", "Q66K64", "Q68EN5", "Q6CXX6", "Q7MVV4", "Q86T03", "Q8N8R7", "Q8NBJ9", "Q96JC1", "Q9BWG6", "Q9D1N4", "Q9KP27", "Q9M1V3", "Q9NUN7", "Q9NXF7"]


    files = os.listdir(QUERY_STRUCTURES_DIR)

    list_query_structures_files_paths =[]     # need to set at least one to empty
    for file in files:
        list_query_structures_files_paths.append(os.path.join(QUERY_STRUCTURES_DIR, file))

    boolean_rotamer_examination = True
    path_output = RESULTS_DIR

    main(list_query_structures_files_paths, boolean_rotamer_examination, path_output)