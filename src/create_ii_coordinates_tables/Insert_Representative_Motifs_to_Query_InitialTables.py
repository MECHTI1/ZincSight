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
    

def main():
    # Establish a connection to the database
    conn = get_db_connection()
    # Create a new cursor
    cur = conn.cursor()
    
    Insert_representative_motif_coordinates_to_User_dataset_coordinates_table (cur,conn)

    # Close the cursor and connection
    cur.close()
    conn.close()


if __name__ == "__main__":
    main()