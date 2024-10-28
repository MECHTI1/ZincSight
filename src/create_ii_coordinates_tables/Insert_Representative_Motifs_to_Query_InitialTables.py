#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:37:53 2023

@author: mechti
"""
import psycopg2
import psycopg2.extras
from src.settings import get_db_connection

def Insert_representative_motif_coordinates_to_User_dataset_coordinates_table(cur,conn):

    
    # SQL command
    Insert_representative_coordinates_table_sql_command = """
    INSERT INTO AF_dataset_coordinates_table_v2
    SELECT * FROM motif_representative_coordinates_table_v2;
    """
    # Execute the SQL command
    cur.execute(Insert_representative_coordinates_table_sql_command)
    # cur.execute("DROP TABLE motif_representative_coordinates_table_v2;")
    
    Insert_representative_detailed_coordinates_table_sql_command = """
    INSERT INTO AF_dataset_detailed_coordinates_table_v2
    SELECT * FROM motif_representative_detailed_coordinates_table_v2;
    """
    #execute the SQL command
    cur.execute(Insert_representative_detailed_coordinates_table_sql_command)
    #cur.execute("DROP TABLE motif_representative_detailed_coordinates_table_v2;")
        
    # Commit the transaction
    conn.commit()
    
    
def Create_indexes_on_II_coordinates_tables(cur,conn):
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
    # commit the transaction
    conn.commit()

def main():
    # Establish a connection to the database
    conn = get_db_connection()
    # Create a new cursor
    cur = conn.cursor()
    
    Insert_representative_motif_coordinates_to_User_dataset_coordinates_table (cur,conn)
    Create_indexes_on_II_coordinates_tables(cur,conn)
        
    # Close the cursor and connection
    cur.close()
    conn.close()


if __name__ == "__main__":
    main()