
"""
Created on Fri Nov 24 11:53:02 2023

@author: mechti
"""


import math
import psycopg2
from src.settings import get_db_connection
import time 

def table_compression():
    def calculate_distance(coord1, coord2):
        """Calculate Euclidean distance between two 3D coordinates."""
        return math.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)
    
    def group_coordinates(coords_score_id):
        """Group tuples containing coordinates, score, and id based on spatial proximity."""
        grouped_coords = []
        processed = set()  # Keep track of processed coordinates
    
        for i, (coord, _, _) in enumerate(coords_score_id):
            if i in processed:
                continue
    
            current_group = [coords_score_id[i]]
            processed.add(i)  #Ifa predcited coord processed alreay and have been put into a group, will not be processed again by the function beacuse already belongs to a group
    
            for j, (other_coord, _, _) in enumerate(coords_score_id):
                if j != i and j not in processed:
                    if all(calculate_distance(other_coord, member_coord[0]) <= 2.6 for member_coord in current_group):          # Check distance with all members in the current group
                        current_group.append(coords_score_id[j])
                        processed.add(j)   #If a predcited coord processed alreay and have been put into a group, will not be processed again by the function beacuse already belongs to a group
    
            grouped_coords.append(current_group)
    
        return grouped_coords
    
    
    
    def select_id_with_lowest_score(groups):
        """Select the row ID with the lowest score from each group in structure (structure have multiple groups)."""
        min_scored_coords_score_id=[]
        whether_was_group_bigger_than_1= False
        for group in groups:

            if len(group)>1:
                whether_was_group_bigger_than_1=True           

            # Filter out None and invalid values from the group
            filtered_group = [item for item in group if item is not None and isinstance(item[1], (int, float))]
    
            if filtered_group:
                # try:
                    # Find the tuple with the lowest score in each group
                    lowest_score_tuple = min(filtered_group, key=lambda x: x[1])
                    # Extract the ID from the tuple
                    min_scored_coords_score_id.append(lowest_score_tuple)
            #     except TypeError as e:
            #         print(f"An error occurred: {e}. Group: {group}")
            #         min_scored_coords_score_id.append(None)
            # # else:
            # #     # Handle the case where the group is empty after filtering
            # #     min_scored_coords_score_id.append(None)
    
        return min_scored_coords_score_id, whether_was_group_bigger_than_1
     
 
    start_time=time.time()

    conn = get_db_connection()
    cur = conn.cursor()
    
    # Step 1: Fetching Data from Database
    cur.execute("SELECT pdbid_alphafoldmodel, metalcoord, score, id FROM scored_af_dataset_with_aggregated_final_tables;")
    data = cur.fetchall()
    
    # Organizing data by pdbid
    pdbid_data = {}
    for pdbid, metalcoord, score, row_id in data:
        if pdbid not in pdbid_data:
            pdbid_data[pdbid] = []
        pdbid_data[pdbid].append((metalcoord, score, row_id))
    
    # Step 2: Clustering and Step 3: Finding Minimum Score in Each Cluster
    list_all_ids_with_min_scores_from_all_structures=[]
    for pdbid, coords_score_id in pdbid_data.items():
       whether_was_group_bigger_than_1= True # Set deafault to True before start the whole process of- 1) grouping to clusters of predicted coords with maximal distance of 2 angstrom between cooords, and 2) select minimal scored coord from each group. check whether can still group with max dis 2 angstrom between predcited coords.
       while whether_was_group_bigger_than_1== True: # did iterations to select the coords Untill no clusters/groups of coord with dis less than 2 angstrom could created
            groups = group_coordinates(coords_score_id) # grouping to clusters of predicted coords with maximal distance of 2 angstrom between cooords
            # print(pdbid, groups)
            coords_score_id,whether_was_group_bigger_than_1= select_id_with_lowest_score(groups)
       
       selected_ids_per_structure= [item[2] for item in coords_score_id]
       list_all_ids_with_min_scores_from_all_structures.extend(selected_ids_per_structure)
    
    # print (list_all_ids_with_min_scores_from_all_structures)
    # print (len(list_all_ids_with_min_scores_from_all_structures))
    if list_all_ids_with_min_scores_from_all_structures:
        print("list_all_ids_with_min_scores_from_all_structures: ", list_all_ids_with_min_scores_from_all_structures)
        query_create_radius_compressed_table = f"""
        CREATE TABLE final_compressed_table_with_scored_binding_sites AS
        SELECT * FROM scored_af_dataset_with_aggregated_final_tables
        WHERE id IN ({','.join(map(str, list_all_ids_with_min_scores_from_all_structures))});
        """
        cur.execute(query_create_radius_compressed_table)
        conn.commit()
    else:
        print("no_predicted_site_in_structure")
    conn.close()
    end_time=time.time()
    print ("time for compression:",end_time- start_time)