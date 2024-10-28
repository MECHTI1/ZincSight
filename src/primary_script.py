#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 18:31:24 2023

@author: mechti
"""
import os
import psycopg2
from src.settings import get_db_connection, QUERY_STRUCTURES_DIR
import time
from src.create_ii_coordinates_tables.ii_coordinates_primary_generator_structures import main as create_ii_coordinates_tables_query_dataset
from src.create_ii_coordinates_tables.Insert_Representative_Motifs_to_Query_InitialTables import main as insert_representative_motifs_to_dataset_tables
from src.motif_search_primary import main as main_first_step
from src.assign_lowest_rmsd_resi_comb import main as main_second_step
from src.dis_angles_calculations.without_rotations.dis_angles_calc_primary_generator_without_his_rotation import main as main_third_step_without_his_rotation
from src.dis_angles_calculations.with_rotations.dis_angles_calc_primary_generator_with_his_rotation_180deg import main as main_third_step_his_rotation_180deg
from src.scoring_and_compression.add_scores_to_table import final_scoring_and_insertion_to_table  
from src.scoring_and_compression.compress_table_by_proximity import table_compression
from src.create_structure_models_with_predicted_zn.primary_create_structure_models_with_predicted_zn import locate_predicted_zn_within_structures
from src.export_final_table_to_csv_format import export_final_table_to_csv_file



def main(list_query_structures_files_paths,boolean_rotamer_examination):
    conn = get_db_connection()
    cur = conn.cursor()
    
    table_creation=True
    His_rotation= boolean_rotamer_examination
    whether_create_structures=True
    
    if table_creation==True:
            start_time_create_tables=time.time()
            
            create_ii_coordinates_tables_query_dataset(list_query_structures_files_paths)

            insert_representative_motifs_to_dataset_tables()
            
            end_time_create_tables=time.time()
            time_create_tables=end_time_create_tables-start_time_create_tables
    
    
         
    start_time_prediction_process_all_sites=time.time()
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
     
    
    # Execute a SELECT statement to get all the site_id values
    cur.execute("SELECT site_id FROM minimized_training_cluster_information;")
    results = cur.fetchall() # Fetch all the results
    
    # Extract the site_id values into a list
    site_ids = [result[0] for result in results]
    print(site_ids)
    
    for site_id in site_ids:
      
        main_first_step(site_id)
        main_second_step()
        
        if His_rotation==True:
            main_third_step_his_rotation_180deg()
              
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
    
    
    
    #Drop table with none type under dif_angle_base
    #This command will not be in the git
    cur.execute("CREATE TABLE dropped_with_erroredrows AS SELECT * FROM AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables WHERE dif_angle_base IS NULL OR rmsd_overall IS NULL;")# and create new table with the dropped rows information which will be called-"dropped_with_erroredrows".
    #This command will be in the git
    cur.execute("DELETE FROM AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables WHERE dif_angle_base IS NULL OR rmsd_overall IS NULL;")# THE DROP COMMAND
    conn.commit()
    
    
    final_scoring_and_insertion_to_table(list_query_structures_files_paths)
        
    # cur.execute("DROP TABLE AF_DATASET_with_metalcoord_test_sites_aggregated_final_tables")
    
    conn.commit()
    table_compression()
    conn.commit()
    
    # cur.execute("DROP TABLE scored_af_dataset_with_aggregated_final_tables")
    #conn.commit()
    conn.commit()
    end_time_prediction_process_all_sites=time.time()
    
    if whether_create_structures== True:
        start_time_create_predcitedmodelstructures = time.time()
        tarred_file_path= locate_predicted_zn_within_structures(conn,list_query_structures_files_paths)
        end_time_create_predcitedmodelstructures = time.time()
    
    
    print ("total time for prediction without create structure models:",end_time_prediction_process_all_sites- start_time_prediction_process_all_sites)    
    if table_creation==True: export_final_table_to_csv_file()
    print ("time_create_II_Coordinates_tables",time_create_tables)
        
    if whether_create_structures== True:
           print ("toatal time for create structures: ", end_time_create_predcitedmodelstructures- start_time_create_predcitedmodelstructures)
    
    cur.close()
    conn.close()
    return (tarred_file_path)


if __name__=="__main__":

   from create_input_argument import primary_return_single_argument_as_paths_list
   
   list_of_uniprot_accessions = ["A0A068N621", "A0A0F6AZI6", "A0A292DHH8", "A0A2U3D0N8", "A0A3F2YM30", "A0A5H1ZR49", "G8ZFK7", "O60232", "P0A6G5", "P0DUH5", "P37659", "P38164", "Q03760", "Q08281", "Q2K0Z2", "Q2UFA9", "Q5W0Q7", "Q66K64", "Q68EN5", "Q6CXX6", "Q7MVV4", "Q86T03", "Q8N8R7", "Q8NBJ9", "Q96JC1", "Q9BWG6", "Q9D1N4", "Q9KP27", "Q9M1V3", "Q9NUN7", "Q9NXF7"]
 
      
   list_query_structures_files_paths=[]
   files = os.listdir(QUERY_STRUCTURES_DIR)
   for file in files:
       list_query_structures_files_paths.append(os.path.join(QUERY_STRUCTURES_DIR, file))
    
   list_query_structures_files_paths =[]     # neeed to set at least one to empty
    


   list_of_paths= primary_return_single_argument_as_paths_list(list_query_structures_files_paths,list_of_uniprot_accessions)
   main(list_of_paths)

            