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
from src.settings import get_db_connection
from src.create_ii_coordinates_tables.extract_structure_II_coords import Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure
from src.create_ii_coordinates_tables.create_insert_indexing_ii_and_coords_tables import create_IIs_and_COORDINATES_TABLES, insert_muliple_rows_from_one_structure as Insert, create_multiple_column_indexes as create_indexes, close_connection

def get_IIs_and_coordinates_lists_of_tuples_and_insert_to_tables(filepath, structure_name,conn,cur,counter):
    start_time_create = time.time()
    IIs_and_coordinates_lists_of_tuples = Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure(
        filepath, structure_name,counter)
    end_time_create = time.time()
    total_time_create = end_time_create - start_time_create

    IIs_list_of_tuples = IIs_and_coordinates_lists_of_tuples[0]
    coordinates_lists_of_tuples = IIs_and_coordinates_lists_of_tuples[1]
    all_relevant_atoms_from_residues_list = IIs_and_coordinates_lists_of_tuples[2]

    start_time_insert = time.time()
    Insert(IIs_list_of_tuples, coordinates_lists_of_tuples,
            all_relevant_atoms_from_residues_list,conn,cur)
    end_time_insert = time.time()
    total_time_insert = end_time_insert - start_time_insert

    return total_time_create, total_time_insert


def multithreaded_iteratations_structures_insert_into_ii_and_coords_tables(list_query_structures_files_paths):

    def process_files(list_query_structures_files_paths,result_queue):
      
        # Each process will create its own database connection
        conn = get_db_connection()
        cur = conn.cursor()

        count = 0

        total_time_create = 0
        total_time_insert = 0
        for filepath in list_query_structures_files_paths:
            count += 1
            base_name = os.path.basename(filepath)   # Get the base name of the file (filename with extension)
            structure_name_without_cif_extension, _ = os.path.splitext(base_name)            # Remove the file extension
            print(structure_name_without_cif_extension)
            print(filepath)
            time_create, time_insert = get_IIs_and_coordinates_lists_of_tuples_and_insert_to_tables(
                filepath, structure_name_without_cif_extension, conn,cur, count)
            total_time_create += time_create
            total_time_insert += time_insert
        result_queue.put((total_time_create, total_time_insert))
        
        return total_time_create, total_time_insert

   
   
    result_queue = multiprocessing.Queue()

    mid_point = len(list_query_structures_files_paths) // 2
    
    # Split the file list into two parts
    first_half = list_query_structures_files_paths[:mid_point]
    second_half =list_query_structures_files_paths[mid_point:]
    
    # Create two processes
    process1 = multiprocessing.Process(target=process_files, args=(first_half,result_queue))
    process2 = multiprocessing.Process(target=process_files, args=(second_half,result_queue))
    
    # Start the processes
    process1.start()
    process2.start()
    
    # Wait for the processes to finish
    process1.join()
    process2.join()
    
    total_time_create = 0
    total_time_insert = 0
    while not result_queue.empty():
        time_create, time_insert = result_queue.get()
        total_time_create += time_create
        total_time_insert += time_insert

    return total_time_create, total_time_insert

def main(list_query_structures_files_paths):

    create_IIs_and_COORDINATES_TABLES()
    start_time = time.time()
    total_time_create, total_time_insert = multithreaded_iteratations_structures_insert_into_ii_and_coords_tables(list_query_structures_files_paths)
    
    start_time_create_indexes=time.time()
    create_indexes()
    end_time_create_indexes= time.time()

    close_connection()
    
    end_time = time.time()
    print("Overall time: ", end_time - start_time)
    print("Total time for Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure calls: ", total_time_create)
    print("Total time for Insert calls: ", total_time_insert)
    print("Total time for create indexes: ", end_time_create_indexes-start_time_create_indexes)



if __name__=="__main__":
    conn = get_db_connection()
    cur = conn.cursor()
    list_query_structures_files_paths=[]
    directory = '../Query AlphaFold structures'
    files = os.listdir(directory)
    for file in files:
        list_query_structures_files_paths.append(os.path.join(directory, file))
        print (list_query_structures_files_paths)
    main(list_query_structures_files_paths)

