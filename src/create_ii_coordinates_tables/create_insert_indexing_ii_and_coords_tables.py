#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 17:13:25 2023

@author: mechti
"""

#Create 2 empty tables (if not existed) of IIs and coordinates  
import psycopg2
import psycopg2.extras
from src.settings import get_db_connection
import json

def create_IIs_and_COORDINATES_TABLES():
    conn = get_db_connection()
    cur = conn.cursor()  

    # SQL command to create a II_TABLE_V2 
    II_TABLE_V2_creation_command = """
    CREATE TABLE IF NOT EXISTS AF_DATASET_II_TABLE_V2(
        PDBID_AlphaFoldModel VARCHAR(255),
        AA_pair VARCHAR(255),
        close_atom_dis FLOAT,
        far_atom_dis FLOAT,
        chain_resi_1 VARCHAR(255),
        chain_resi_2 VARCHAR(255)
    );
    """
    
    # SQL command to create a COORDINATES TABLE
    #Real[3](3) -->3 eleements in arrat and 3 number after decimal point
    COORDINATES_TABLE_V2_creation_command = """
    CREATE TABLE IF NOT EXISTS AF_DATASET_COORDINATES_TABLE_V2(
        PDBID_AlphaFoldModel VARCHAR(255),
        chain_resi VARCHAR(255),
        close_atom_coord REAL[], 
        far_atom_coord REAL[],
        B_factor SMALLINT
    );
    """
    
    DETAILED_COORDINATES_TABLE_V2_creation_command = """
    CREATE TABLE IF NOT EXISTS AF_DATASET_DETAILED_COORDINATES_TABLE_V2(
        PDBID_AlphaFoldModel VARCHAR(255),
        chain_resi VARCHAR(255),
        resi_type VARCHAR(10),
        dict_atom_coord json
    );
    """
    
    # Execute the commands
    cur.execute(II_TABLE_V2_creation_command)
    cur.execute(COORDINATES_TABLE_V2_creation_command)
    cur.execute(DETAILED_COORDINATES_TABLE_V2_creation_command)

    # Commit the transaction
    conn.commit()
    # Close the cursor and connection
    cur.close()
    conn.close()

def insert_muliple_rows_from_one_structure(list_of_II_tuples, list_of_coordinates_tuples,all_relevant_atoms_from_residues_list,conn,cur):
    insert_query_II_TABLE_V2 = "INSERT INTO AF_DATASET_II_TABLE_V2 (PDBID_AlphaFoldModel, AA_pair, close_atom_dis, far_atom_dis, chain_resi_1, chain_resi_2) VALUES %s"
    insert_query_COORDINATES_TABLE_V2= "INSERT INTO AF_DATASET_COORDINATES_TABLE_V2 (PDBID_AlphaFoldModel, chain_resi, close_atom_coord, far_atom_coord,B_factor) VALUES %s"
    insert_query_DETAILED_COORDINATES_TABLE_V2= "INSERT INTO AF_DATASET_DETAILED_COORDINATES_TABLE_V2 (PDBID_AlphaFoldModel, chain_resi,resi_type, dict_atom_coord) VALUES %s"
    
    
    all_relevant_atoms_from_residues_list_converted_list = [(pdb, chain, resi, json.dumps(atom_dict)) for pdb, chain, resi, atom_dict in all_relevant_atoms_from_residues_list]
    #print (all_relevant_atoms_from_residues_list_converted_list)


    # Execute the INSERT query
    psycopg2.extras.execute_values(cur, insert_query_II_TABLE_V2, list_of_II_tuples)
    psycopg2.extras.execute_values(cur, insert_query_COORDINATES_TABLE_V2, list_of_coordinates_tuples)
    psycopg2.extras.execute_values(cur, insert_query_DETAILED_COORDINATES_TABLE_V2, all_relevant_atoms_from_residues_list_converted_list)
    
    
    # Commit the transaction
    conn.commit()



#create multiple column indexes. three combinations (for three searches- one in the start of the II combination chain, the second and three in the following ssearches in the chain)
def create_multiple_column_indexes():
    conn = get_db_connection()
    cur = conn.cursor()

    # SQL commands to create indexes
    index_creation_commands = [
        "CREATE INDEX idx_1_v2_II ON AF_DATASET_II_TABLE_V2 (AA_pair, close_atom_dis, far_atom_dis)",
        "CREATE INDEX idx_2_v2_II ON AF_DATASET_II_TABLE_V2 (PDBID_AlphaFoldModel, AA_pair, close_atom_dis, far_atom_dis, chain_resi_1)",
        "CREATE INDEX idx_3_v2_II ON AF_DATASET_II_TABLE_V2 (PDBID_AlphaFoldModel, AA_pair, close_atom_dis, far_atom_dis, chain_resi_2)",
        "CREATE INDEX idx_1_v2_COORD ON AF_DATASET_COORDINATES_TABLE_V2 (PDBID_AlphaFoldModel, chain_resi)",
        "CREATE INDEX idx_1_COORDINATES_v2_COORD ON AF_DATASET_DETAILED_COORDINATES_TABLE_V2 (PDBID_AlphaFoldModel, chain_resi)"
    ]
    
    for command in index_creation_commands:
        try:
            cur.execute(command)
            conn.commit()
        except psycopg2.errors.DuplicateTable:  # The error that PostgreSQL raises when an index already exists
            print(f"Index already exists: {command}")
            conn.rollback()  # Rollback the current transaction to recover from the error

   
        
    cur.close()
    # Commit the transaction
    conn.commit()
    conn.close()















