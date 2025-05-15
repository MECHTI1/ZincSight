#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:46:38 2024

@author: mechti
"""

import os
import psycopg2
from dotenv import load_dotenv



# Set the base directory to the absolute path of the project directory
base_directory =os.path.abspath(os.path.dirname(__file__))
PARENT_DIRECTORY = os.path.abspath(os.path.join(base_directory, os.pardir))

QUERY_STRUCTURES_DIR = os.path.join(PARENT_DIRECTORY, 'Query AlphaFold structures')

RESULTS_DIR= os.path.join(PARENT_DIRECTORY,'results')
TABLES_DIR=os.path.join(RESULTS_DIR,'table')
STRUCTURES_WITH_PREDICTED_ZN= os.path.join(RESULTS_DIR,'structures_with_predicted_zn')

os.makedirs(QUERY_STRUCTURES_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(STRUCTURES_WITH_PREDICTED_ZN, exist_ok=True)

load_dotenv() # Load variables from the .env file into the environment

def get_env_bool(name, default="false"):
    return os.getenv(name, default).strip().lower() == "true"   # returns True only if the string is exactly "true" (case-insensitive).

DEBUGGING = get_env_bool("DEBUGGING_TEMP")
KEEP_TEMP_TABLES = get_env_bool("KEEP_TEMP_TABLES_TEMP")

""" Filtering out predictions that have prob below threshold val stored in EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD var.       
    If FILTER_PROB_AFTER_COMPRESSION is assigned as true: will executed after compression,
 otherwise nothing in the meanwhile""" #TODO:maybe can add option to remove before gor maybe enhance
FILTER_PROB_AFTER_COMPRESSION = get_env_bool("FILTER_PROB_AFTER_COMPRESSION")
MIN_THRESHOLD_PROB = float(os.getenv("EXCLUDE_PREDICTIONS_WITH_PROB_THRESHOLD") or "0.0")   # Retrieve the threshold value and convert it to a float

def get_db_connection():
            load_dotenv()  # Load environment variables from the .env file
            user_id = os.getenv("DB_USER")
            password = os.getenv("PASSWORD")
            host = os.getenv("HOST")
            conn = psycopg2.connect(dbname="zincsight_pipeline_db", user=user_id, host=host, password=password)

            return conn

def cleanup_tables(cur,conn):
    """Keep only required tables, delete others"""
    required_tables = {
        'minimized_training_cluster_information',
        'motif_representative_coordinates_table_v2',
        'motif_representative_detailed_coordinates_table_v2',
        'training_representative_metal_sites_kruskal_v2'
    }
    cur.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'")
    for table in cur.fetchall():
        if table[0] not in required_tables:
            cur.execute(f'DROP TABLE IF EXISTS "{table[0]}" CASCADE')
    conn.commit()

if __name__== "__main__":
    conn = get_db_connection()
