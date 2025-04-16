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

DEBUGGING = get_env_bool("DEBUGGING")
KEEP_TEMP_TABLES = get_env_bool("KEEP_TEMP_TABLES")

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



if __name__== "__main__":
    conn = get_db_connection()
