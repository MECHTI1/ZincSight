#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 21:02:54 2024

@author: mechti
"""

import psycopg2
from src.settings import get_db_connection

# Path to the SQL file (ensure this is the correct path in your Colab environment)
sql_tables_file_path = 'PostgreSQL_4_necessary_tables.sql'
# Read the dump file
with open(dump_file_path, 'r') as f:
    sql = f.read()
try:
    # Read the SQL file
    with open(sql_tables_file_path, 'r') as file:
        sql_commands_create_tables = file.read()

    conn = get_db_connection()
    with conn.cursor() as cursor:
        cursor.execute(sql_commands_create_tables)
        conn.commit()
        print("Tables created successfully!")
    conn.close()

except (IOError, psycopg2.Error, Exception) as e:
    print(f"Error: {e}")
