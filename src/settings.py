#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 08:46:38 2024

@author: mechti
"""

import os
import psycopg2
import random
import string
import subprocess
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


DROP_TEMP_TABLES = bool(os.getenv('DROP_TEMP_TABLES'))


def get_db_connection():
    owner = True
    try:
        if owner == True:
            load_dotenv()  # Load environment variables from the .env file
            user_id = os.getenv('DB_USER')
            password = os.getenv('PASSWORD')
            host = os.getenv('HOST')
            conn = psycopg2.connect(dbname="zincsight_pipeline_db", user=user_id, host=host, password=password)

        else:
            conn = psycopg2.connect(dbname="zincsight_pipeline_db", user="postgres", host="localhost")  # No password

        return conn

    except psycopg2.Error as e:
        print(f"Error connecting to the database: {e}")
        return None


if __name__== "__main__":
    conn = get_db_connection()
