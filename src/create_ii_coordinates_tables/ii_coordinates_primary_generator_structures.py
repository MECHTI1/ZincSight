#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:39:13 2023

@author: mechti
"""

import time
import os
import sys
import multiprocessing
import psycopg2
from tqdm.contrib.concurrent import process_map

from src.settings import get_db_connection
from src.create_ii_coordinates_tables.extract_structure_II_coords import Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure
from src.create_ii_coordinates_tables.create_insert_indexing_ii_and_coords_tables import create_IIs_and_COORDINATES_TABLES, insert_muliple_rows_from_one_structure as Insert, create_multiple_column_indexes as create_indexes

def get_IIs_and_coordinates_lists_of_tuples_and_insert_to_tables(filepath, structure_name):
    IIs_and_coordinates_lists_of_tuples = Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure(
        filepath, structure_name)

    IIs_list_of_tuples = IIs_and_coordinates_lists_of_tuples[0]
    coordinates_lists_of_tuples = IIs_and_coordinates_lists_of_tuples[1]
    all_relevant_atoms_from_residues_list = IIs_and_coordinates_lists_of_tuples[2]

    Insert(IIs_list_of_tuples, coordinates_lists_of_tuples,
            all_relevant_atoms_from_residues_list)


def structures_insert_into_ii_and_coords_tables(filepath):
    base_name = os.path.basename(filepath)  # Get the base name of the file (filename with extension)
    structure_name_without_cif_extension, _ = os.path.splitext(base_name)  # Remove the file extension
    print(structure_name_without_cif_extension)
    print(filepath)
    get_IIs_and_coordinates_lists_of_tuples_and_insert_to_tables(
        filepath, structure_name_without_cif_extension)

def main(list_query_structures_files_paths):
    create_IIs_and_COORDINATES_TABLES()
    process_map(structures_insert_into_ii_and_coords_tables, list_query_structures_files_paths,
                max_workers=os.cpu_count(), chunksize=10)
    create_indexes()


if __name__ == '__main__':
    conn = get_db_connection()
    cur = conn.cursor()
    list_query_structures_files_paths=[]
    directory = '../Query AlphaFold structures'
    files = os.listdir(directory)
    for file in files:
        list_query_structures_files_paths.append(os.path.join(directory, file))
        print (list_query_structures_files_paths)
    main(list_query_structures_files_paths)

