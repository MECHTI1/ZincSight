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

# Load variables from the .env file into the environment
load_dotenv()
# KEEP_TEMP_TABLES = bool(os.getenv('KEEP_TEMP_TABLES'))
KEEP_TEMP_TABLES = True
def get_db_connection():
            load_dotenv()  # Load environment variables from the .env file
            # Retrieve the variables
            user_id = os.getenv("DB_USER")
            password = os.getenv("PASSWORD")
            host = os.getenv("HOST")
            conn = psycopg2.connect(dbname="zincsight_pipeline_db", user=user_id, host=host, password=password)

            return conn



if __name__== "__main__":
    conn = get_db_connection()
