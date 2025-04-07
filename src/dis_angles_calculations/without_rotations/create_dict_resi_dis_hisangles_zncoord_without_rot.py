#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 10:40:22 2023

@author: mechti
"""
import numpy as np
import time
from scipy.optimize import least_squares
from src.dis_angles_calculations.without_rotations.histidines_other_calc_given_point_v5_v2 import check_candidate_metal_coord_valid_histidines_angles
from src.dis_angles_calculations.without_rotations.calculate_angles_between_vectors import compute_angles_between_vectors
import traceback

#add- /home/mechti/Documents/PhD project/Create_tables_postgresql/Generate_tables_of_metal_binding_sites/Motif_Search/check_search_motif_algo_results/middle_point_every_shape.py
#add- /home/mechti/Documents/PhD project/Create_tables_postgresql/Generate_tables_of_metal_binding_sites/Motif_Search/check_search_motif_algo_results/histidines_other_calc_given_point_v2.py
aa_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'} 


def calculate_point_at_distance_more_then_3_binding_points(points, distance):
    def equations(vars):
        x, y, z = vars
        return [(x - point["coords"][0])**2 + (y - point["coords"][1])**2 + (z - point["coords"][2])**2 - point["resi_type_dis"]**2 for point in points]

    # Calculate the centroid of the points
    centroid = np.mean([point["coords"]  for point in points], axis=0)


    # Initial guess - at the centroid of the points
    x0 = centroid
    res = least_squares(equations, x0, method='lm')

    return res.x


def calculate_point_at_distance(A, B, C, distance):
    def equations(vars, points, distance):
        x, y, z = vars
        # print ("q",points[0]["coords"][0])
        # print ("w", points[0]["resi_type_dis"])
        return [(x - point["coords"][0])**2 + (y - point["coords"][1])**2 + (z - point["coords"][2])**2 - point["resi_type_dis"]**2 for point in points]

    # Calculate the centroid of the triangle
    centroid = (A ["coords"]+ B["coords"] + C["coords"]) / 3

    # Calculate the normal to the plane of the triangle
    normal = np.cross(B["coords"] - A["coords"], C["coords"] - A["coords"])
    normal = normal / np.linalg.norm(normal)

    # First initial guess - a bit below the centroid of the points
    x0_below = centroid - normal * 0.1  # Adjust this value as needed
    res1 = least_squares(equations, x0_below, args=((A, B, C), distance), method='lm')

    # Second initial guess - a bit above the centroid
    x0_above = centroid + normal * 0.1  # Adjust this value as needed
    res2 = least_squares(equations, x0_above, args=((A, B, C), distance), method='lm')
    # print (res2.x,res1.x) 
    return [res1.x, res2.x]



def Choosing_coordinating_atoms(dict_residues_in_candidate_motif):
    
    lst_mid_points= []
    for resi in dict_residues_in_candidate_motif:
        if dict_residues_in_candidate_motif [resi][0] == 'HIS':
            NE2_coord = dict_residues_in_candidate_motif [resi][1].get ('NE2')
            ND1_coord = dict_residues_in_candidate_motif [resi][1].get ('ND1')
            midpoint =  (NE2_coord + ND1_coord) /2 
            lst_mid_points.append (midpoint)
        if dict_residues_in_candidate_motif [resi][0] == 'GLU':
            OE1_coord = dict_residues_in_candidate_motif [resi][1].get ('OE1')
            OE2_coord = dict_residues_in_candidate_motif [resi][1].get ('OE2')
            midpoint =  (OE1_coord + OE2_coord) /2 
            lst_mid_points.append (midpoint)
        if dict_residues_in_candidate_motif [resi][0] == 'ASP':
            OD1_coord = dict_residues_in_candidate_motif [resi][1].get ('OD1')
            OD2_coord = dict_residues_in_candidate_motif [resi][1].get ('OD2')
            midpoint =  (OD1_coord + OD2_coord) /2 
            lst_mid_points.append (midpoint)
            
        if dict_residues_in_candidate_motif [resi][0] == 'THR':
            OG1_coord = dict_residues_in_candidate_motif [resi][1].get ('OG1')
            CG2_coord = dict_residues_in_candidate_motif [resi][1].get ('CG2')
            midpoint =  (OG1_coord + CG2_coord) /2 
            lst_mid_points.append (midpoint)

        if dict_residues_in_candidate_motif [resi][0] == 'ASN':
            OD1_coord = dict_residues_in_candidate_motif [resi][1].get ('OD1')
            ND2_coord = dict_residues_in_candidate_motif [resi][1].get ('ND2')
            midpoint =  (OD1_coord + ND2_coord) /2 
            lst_mid_points.append (midpoint)

            
        if dict_residues_in_candidate_motif [resi][0] == 'CYS':
            SG_coord = dict_residues_in_candidate_motif [resi][1].get ('SG')
            lst_mid_points.append (SG_coord)
        if dict_residues_in_candidate_motif [resi][0] == 'TYR':
            OH_coord = dict_residues_in_candidate_motif [resi][1].get ('OH')
            lst_mid_points.append (OH_coord)
        if dict_residues_in_candidate_motif [resi][0] == 'MET':
            SD_coord = dict_residues_in_candidate_motif [resi][1].get ('SD')
            lst_mid_points.append (SD_coord)
        if dict_residues_in_candidate_motif [resi][0] == 'SER':
            OG_coord = dict_residues_in_candidate_motif [resi][1].get ('OG')
            lst_mid_points.append (OG_coord)
            




    #print (lst_mid_points)       
    # convert the list of arrays into a 2D numpy array
    stacked_arrays = np.stack(lst_mid_points)
    
    # compute the mean along the first axis (which corresponds to your original lists)
    semi_final_middle_point = np.mean(stacked_arrays, axis=0)
    
    #print("semi_final_middle_point", semi_final_middle_point)  
    
    
    for resi in dict_residues_in_candidate_motif:
      if dict_residues_in_candidate_motif[resi][0] == 'HIS':
          NE2_coord = dict_residues_in_candidate_motif[resi][1].get('NE2')
          ND1_coord = dict_residues_in_candidate_motif[resi][1].get('ND1')
          # CE1_coord = dict_residues_in_candidate_motif[resi][1].get('CE1')
          # CD2_coord = dict_residues_in_candidate_motif[resi][1].get('CD2')
          dist_NE2 = np.linalg.norm(semi_final_middle_point - NE2_coord)
          dist_ND1 = np.linalg.norm(semi_final_middle_point - ND1_coord)
          # dist_CE1 = np.linalg.norm(semi_final_middle_point - CE1_coord)
          # dist_CD2 = np.linalg.norm(semi_final_middle_point -CD2_coord)
          # Creating a dictionary with the distances
          # distances = {'NE2': dist_NE2,
          #              'ND1': dist_ND1,
          #              'CE1': dist_CE1,
          #              'CD2': dist_CD2}
          distances = {'NE2': dist_NE2,
                        'ND1': dist_ND1}
          
          closest_atom = min(distances, key=distances.get)            # Finding the atom with the minimum distance

          dict_residues_in_candidate_motif[resi].append(closest_atom) #The third item in the list is the colsest atom

          #print(f"For {resi}, the closest atom to the semi_final_middle_point is {closest_atom}")
          
      if dict_residues_in_candidate_motif[resi][0] == 'GLU':
          OE1_coord = dict_residues_in_candidate_motif[resi][1].get('OE1')
          OE2_coord = dict_residues_in_candidate_motif[resi][1].get('OE2')
          dist_OE1 = np.linalg.norm(semi_final_middle_point - OE1_coord)
          dist_OE2 = np.linalg.norm(semi_final_middle_point - OE2_coord)
          closest_atom = 'OE2' if dist_OE2 < dist_OE1 else 'OE1'
          dict_residues_in_candidate_motif[resi].append(closest_atom) #The third item in the list is the colsest atom

          #print(f"For {resi}, the closest atom to the semi_final_middle_point is {closest_atom}")

      if dict_residues_in_candidate_motif [resi][0] == 'ASP':
          OD1_coord = dict_residues_in_candidate_motif[resi][1].get('OD1')
          OD2_coord = dict_residues_in_candidate_motif[resi][1].get('OD2')
          dist_OD1 = np.linalg.norm(semi_final_middle_point - OD1_coord)
          dist_OD2 = np.linalg.norm(semi_final_middle_point - OD2_coord)
          closest_atom = 'OD2' if dist_OD2 < dist_OD1 else 'OD1'
          dict_residues_in_candidate_motif[resi].append(closest_atom) #The third item in the list is the colsest atom
          
      if dict_residues_in_candidate_motif [resi][0] == 'THR':
            OG1_coord = dict_residues_in_candidate_motif[resi][1].get('OG1')
            CG2_coord = dict_residues_in_candidate_motif[resi][1].get('CG2')
            dist_OG1 = np.linalg.norm(semi_final_middle_point - OG1_coord)
            dist_CG2 = np.linalg.norm(semi_final_middle_point - CG2_coord)
            closest_atom = 'CG2' if dist_CG2 < dist_OG1 else 'OG1'
            dict_residues_in_candidate_motif[resi].append(closest_atom) #The third item in the list is the colsest atom
          
            #print(f"For {resi}, the closest atom to the semi_final_middle_point is {closest_atom}")

      if dict_residues_in_candidate_motif [resi][0] == 'ASN':
            OD1_coord = dict_residues_in_candidate_motif[resi][1].get('OD1')
            ND2_coord = dict_residues_in_candidate_motif[resi][1].get('ND2')
            dist_OD1 = np.linalg.norm(semi_final_middle_point - OD1_coord)
            dist_ND2 = np.linalg.norm(semi_final_middle_point - ND2_coord)
            closest_atom = 'ND2' if dist_ND2 < dist_OD1 else 'OD1'
            dict_residues_in_candidate_motif[resi].append(closest_atom) #The third item in the list is the colsest atom          
          


      if dict_residues_in_candidate_motif [resi][0] == 'CYS':
          dict_residues_in_candidate_motif[resi].append('SG')
      if dict_residues_in_candidate_motif [resi][0] == 'TYR':
          dict_residues_in_candidate_motif[resi].append('OH')
      if dict_residues_in_candidate_motif [resi][0] == 'MET':
          dict_residues_in_candidate_motif[resi].append('SD')
      if dict_residues_in_candidate_motif [resi][0] == 'SER':
          dict_residues_in_candidate_motif[resi].append('OG')
          

                    
        
          
     
          
                  


    return (dict_residues_in_candidate_motif)   # {'A_171': ['HIS', {'NE2': array([-6.783,  1.6  ,  7.082]), 'ND1': array([-7.078,  3.727,  7.531])}, 'NE2'], 'A_272': ['HIS', {'NE2': array([-8.973,  0.214,  5.061]), 'ND1': array([11.106,  0.169,  5.592])}, 'NE2']}


def calculate_total_difference(candidate_point_angles_stast):
    return sum(candidate_point_angles_stast['candidate_point_angles_stast'].values())




# dict_residues_in_candidate_motif = {
#     'A_166': ['HIS', {
#         'NE2': np.array([-0.363, -0.177, -8.548]),
#         'ND1': np.array([-2.323, -1.159, -8.446]),
#         'CE1': np.array([-1.067, -1.159, -7.970]),
#         'CD2': np.array([-1.238, 0.490, -9.395]),
#         'CG': np.array([-2.455, -0.115, -9.329])
#     }],
#     'A_180': ['HIS', {
#         'NE2': np.array([2.408, 0.017, -10.478]),
#         'ND1': np.array([2.248, 0.467, -12.599]),
#         'CE1': np.array([1.574, 0.081, -11.507]),
#         'CD2': np.array([3.655, 0.440, -10.910]),
#         'CG': np.array([3.556, 0.709, -12.244])
#     }],
#     'A_184': ['GLU', {
#         'OE1': np.array([2.087, 2.208, -8.224]),
#         'OE2': np.array([4.216, 1.473, -8.249])
#     }]
# }


# dict_residues_in_candidate_motif = {
#     'A_107': ['HIS', {
#         'CG': np.array([-2.161, -5.437, -5.642]),
#         'CD2': np.array([-2.261, -6.764, -5.319]),
#         'ND1': np.array([-1.572, -4.773, -4.580]),
#         'CE1': np.array([-1.357, -5.683, -3.629]),
#         'NE2': np.array([-1.763, -6.892, -4.043])
#     }],
#     'A_156': ['HIS', {
#         'CG': np.array([-2.907, -2.639, -0.216]),
#         'CD2': np.array([-2.856, -3.071, -1.509]),
#         'ND1': np.array([-1.677, -2.122, 0.154]),
#         'CE1': np.array([-0.894, -2.275, -0.926]),
#         'NE2': np.array([-1.574, -2.838, -1.938])
#     }],
#     'A_160': ['ASP', {
#         'CG': np.array([-3.262, -1.250, -4.964]),
#         'OD1': np.array([-2.186, -0.617, -4.942]),
#         'OD2': np.array([-3.499, -2.364, -4.459])
#     }],
#     'A_93': ['ASP', {
#         'CG': np.array([1.796, -4.155, -2.054]),
#         'OD1': np.array([1.340, -4.605, -3.121]),
#         'OD2': np.array([1.903, -2.956, -1.756])
#     }]
# }





# dict_residues_in_candidate_motif = {
#     'A_55': ['ASP', {
#         'OD1': np.array([39.092, 19.909, 7.029]),
#         'OD2': np.array([41.131, 19.639, 6.336])
#     }],
#     'A_69': ['HIS', {
#         'CG': np.array([36.103, 22.677, 8.102]),
#         'CD2': np.array([35.042, 21.875, 7.997]),
#         'ND1': np.array([37.137, 21.804, 8.229]),
#         'CE1': np.array([36.760, 20.545, 8.200]),
#         'NE2': np.array([35.451, 20.570, 8.060])
#     }],
#     'A_122': ['ASP', {
#         'OD1': np.array([41.260, 24.248, 9.765]),
#         'OD2': np.array([39.156, 23.415, 9.817])
#     }],
#     'A_118': ['HIS', {
#         'CG': np.array([39.775, 19.425, 11.866]),
#         'CD2': np.array([39.236, 20.460, 11.168]),
#         'ND1': np.array([40.582, 18.720, 11.012]),
#         'CE1': np.array([40.552, 19.280, 9.808]),
#         'NE2': np.array([39.755, 20.299, 9.919])
#     }]



# dict_residues_in_candidate_motif = {
#     'A_125': ['CYS', {'SG': [-22.835, 8.634, 6.667]}], 
#     'A_128': ['CYS', {'SG': [-22.924, 4.868, 7.574]}], 
#     'A_90': ['CYS', {'SG': [-25.438, 6.124, 5.579]}], 
#     'A_93': ['CYS', {'SG': [-25.521, 7.357, 9.181]}]
# }




        
def caclulate_dis_CoordinationAngles_HISangles_stats (dict_residues_in_candidate_motif):
 
  try:

    # Sorting keys based on whether the value starts with 'CYS'
    sorted_keys = sorted(dict_residues_in_candidate_motif.keys(), key=lambda x: dict_residues_in_candidate_motif[x][0] != 'CYS')

    # Reconstructing dictionary with sorted keys
    dict_residues_in_candidate_motif = {key: dict_residues_in_candidate_motif[key] for key in sorted_keys}

    # print(dict_residues_in_candidate_motif)
    # start_time = time.time()

    # Convert the numeric data to numpy arrays
    for key, value in dict_residues_in_candidate_motif.items():
        for subkey, subvalue in value[1].items():
            dict_residues_in_candidate_motif[key][1][subkey] = np.array(subvalue)
    
    Choosing_coordinating_atoms(dict_residues_in_candidate_motif) #change the dict_residues_in_candidate_motif to have also the binding atoms
    #print ("dict_residues_in_candidate_motif",dict_residues_in_candidate_motif)
    
    lst_resi_binding_atoms_coordinates=[]
    for resi in dict_residues_in_candidate_motif:
        binding_atom= dict_residues_in_candidate_motif[resi][2]
        coords_binding_atom= dict_residues_in_candidate_motif[resi][1][binding_atom]
        lst_resi_binding_atoms_coordinates.append({
    "resi_type_dis": 2.32 if dict_residues_in_candidate_motif[resi][0] == 'CYS' else 2.15, 
    "coords": coords_binding_atom
})
    # print (lst_resi_binding_atoms_coordinates)
    
     # Distance from each vertex
    total_dis=0
    for resi in dict_residues_in_candidate_motif:
         if dict_residues_in_candidate_motif [resi][0] == 'CYS':
             total_dis+=2.32
         else:
             total_dis+=2.15
             
    distance= total_dis/ (len(dict_residues_in_candidate_motif))
    #print (distance)
    # tolerance_distance = 0.5  # Set a suitable 
     
    if len (dict_residues_in_candidate_motif) <4:
        A = lst_resi_binding_atoms_coordinates[0]
        B = lst_resi_binding_atoms_coordinates[1]
        C = lst_resi_binding_atoms_coordinates[2]
        
        # print ("a,b,c", A,B,C)
        mid_points = calculate_point_at_distance(A, B, C, distance)
        
      # Boolean_correctness=True
     
        tolerance_distinct=0.01
        
        # for i, point in enumerate(mid_points, start=1):
        #     #print(f"E{i}: ", point)
        #     distances = [np.linalg.norm(point - A["coords"]), np.linalg.norm(point - B["coords"]), np.linalg.norm(point - C["coords"])]
        #     #print("Distances from E{i} to A, B, C: ", distances)
        # # if all(abs(d - distance) <= tolerance_distance for d in distances):
        #     #print(f"All distances from E{i} are within the tolerance of the desired distance.")
        # else:
            #print(f"At least one distance from E{i} is not within the tolerance of the desired distance. maybe delete the match")
           # Boolean_correctness=False
            # if all(abs(d - distance) <= tolerance_distance for d in distances):
            #     #print(f"All distances from E{i} are within the tolerance of the desired distance.")
            # # else:
            #     for i in range(len(lst_resi_binding_atoms_coordinates)):
            #         for j in range(i+1, len(lst_resi_binding_atoms_coordinates)):
            #             distance = np.linalg.norm(lst_resi_binding_atoms_coordinates[i]-lst_resi_binding_atoms_coordinates[j])
            #             if distance > 6:
                            #print(f"Distance between point {i+1} and point {j+1} is greater than 6 Angstroms: {distance} Angstroms")
                        # else:
                            #print(f"Distance between point {i+1} and point {j+1} is not greater than 6 Angstroms: {distance} Angstroms")
        
        # Check if the found points are distinct
        # if np.allclose(mid_points[0], mid_points[1], atol=tolerance_distinct):
            #print("The found points are not distinct.Which means maybe planar")
            #Boolean_correctness=False
        # else:
            #print("The found points are distinct")
     
    
        dict_point_1_candidate_point_angles_stast= check_candidate_metal_coord_valid_histidines_angles (mid_points[0], dict_residues_in_candidate_motif)
        dict_point_2_candidate_point_angles_stast= check_candidate_metal_coord_valid_histidines_angles (mid_points[1], dict_residues_in_candidate_motif)
        
        # calculation SUM dif angles_stast
        total_difference_point_1 = calculate_total_difference(dict_point_1_candidate_point_angles_stast)
        total_difference_point_2 = calculate_total_difference(dict_point_2_candidate_point_angles_stast)      

        # comparison and output
        if total_difference_point_1 < total_difference_point_2:
            #print("Point 1 is better with a total difference of", total_difference_point_1)
            dict_most_reasonable_candidate_angle_stats= dict_point_1_candidate_point_angles_stast
            mid_p= mid_points[0]
        else:
            #print("Point 2 is better with a total difference of", total_difference_point_2)
            dict_most_reasonable_candidate_angle_stats= dict_point_2_candidate_point_angles_stast
            mid_p= mid_points[1]
    
        distances_list = [float(np.linalg.norm(mid_p - A["coords"])), float(np.linalg.norm(mid_p - B["coords"])), float(np.linalg.norm(mid_p - C["coords"]))]

    
    
    else:
        
        
        
        mid_p = calculate_point_at_distance_more_then_3_binding_points(lst_resi_binding_atoms_coordinates, distance)
        coords_list = [d["coords"] for d in lst_resi_binding_atoms_coordinates]
        # START_ANGLE_TIME=time.time()
        #print("mid_p: ", mid_p)
        distances_list=[]
        for i, point_i in enumerate(coords_list, start=1):
            distance_to_point = float(np.linalg.norm(mid_p - point_i))
            distances_list.append (distance_to_point)
            #print(f"Distance from the point to A{i}: ", distance_to_point)
            # if not (distance - tolerance_distance <= distance_to_point <= distance + tolerance_distance):
                #print(f"The distance from point A{i} to the calculated point is outside the tolerance.")
                
        dict_most_reasonable_candidate_angle_stats= check_candidate_metal_coord_valid_histidines_angles (mid_p, dict_residues_in_candidate_motif)
    
    
    
    # Extracting only the coordinates from each dictionary
    coords_list = [d["coords"] for d in lst_resi_binding_atoms_coordinates]
    # START_ANGLE_TIME=time.time()
    angles_between_vectors= compute_angles_between_vectors(coords_list, mid_p)
    # END_ANGLE_TIME=time.time()
    # print(END_ANGLE_TIME-START_ANGLE_TIME)
    #print (distances_list)
    #print (dict_most_reasonable_candidate_angle_stats)
    #print (angles_between_vectors)
            
    # end_time = time.time()
    
    #print("Execution time:", end_time - start_time, "seconds")
    
    list_of_aa = [aa_dict[residue_info[0]] for residue_info in dict_residues_in_candidate_motif.values()]
    sorted_list_of_aa = sorted(list_of_aa)     # Sort the list
    

    # Convert numpy array to a Python list
    list_format_metalcoord = mid_p.tolist()

    
    return ({"sorted_list_of_aa": sorted_list_of_aa,"distances_list":distances_list, "HIS_anlges_stats":dict_most_reasonable_candidate_angle_stats, "Coordination_anlges": angles_between_vectors, "metalcoord":list_format_metalcoord})
    #assess points by angles 
    #calculate_mid_point

  except Exception as e:
    # print("An error occurred:")
    print(e)
    traceback.print_exc()  # This will print the traceback including the line number

    # print (dict_residues_in_candidate_motif)
    return ("error")

# dict_residues_in_candidate_motif = {
#     'A_125': ['CYS', {'SG': np.array([-22.835, 8.634, 6.667])}], 
#     'A_128': ['CYS', {'SG': np.array([-22.924, 4.868, 7.574])}], 
#     'A_90': ['CYS', {'SG': np.array([-25.438, 6.124, 5.579])}], 
#     'A_93': ['CYS', {'SG': np.array([-25.521, 7.357, 9.181])}]
# }

# dict_residues_in_candidate_motif = {
#     'A_125': ['CYS', {'SG': [-22.835, 8.634, 6.667]}], 
#     'A_128': ['CYS', {'SG': [-22.924, 4.868, 7.574]}], 
#     'A_90': ['CYS', {'SG': [-25.438, 6.124, 5.579]}], 
#     'A_93': ['CYS', {'SG': [-25.521, 7.357, 9.181]}]
# }


# dict_residues_in_candidate_motif = {
#     'A_123': ['HIS', {
#         'CG': np.array([5.787, -7.635, -3.778]),
#         'CD2': np.array([6.796, -8.343, -4.379]),
#         'ND1': np.array([4.767, -7.543, -4.694]),
#         'CE1': np.array([5.148, -8.172, -5.816]),
#         'NE2': np.array([6.359, -8.738, -5.651])}],  
#     'A_128': ['CYS', {'SG': [2.247,-7.357,-2.372]}], 
#     'A_90': ['CYS', {'SG': [1.407, -8.539, -5.768]}]
# }
# print (caclulate_dis_CoordinationAngles_HISangles_stats(dict_residues_in_candidate_motif))
########################################################################################################################3
# start_time=time.time()

# dict_residues_in_candidate_motif = {
#     'A_125': ['CYS', {'SG': [-22.835, 8.634, 6.667]}], 
#     'A_128': ['CYS', {'SG': [-22.924, 4.868, 7.574]}], 
#     'A_90': ['CYS', {'SG': [-25.438, 6.124, 5.579]}], 
#     'A_93': ['CYS', {'SG': [-25.521, 7.357, 9.181]}]
# }

# dict_residues_in_candidate_motif = {
#     'A_125': ['CYS', {'SG': [-22.835, 8.634, 6.667]}], 
#     'A_128': ['CYS', {'SG': [-22.924, 4.868, 7.574]}], 
#     'A_90': ['CYS', {'SG': [-25.438, 6.124, 5.579]}],
# }
# dict_residues_in_candidate_motif = {
#     'A_123': ['HIS', {
#         'CG': np.array([5.787, -7.635, -3.778]),
#         'CD2': np.array([6.796, -8.343, -4.379]),
#         'ND1': np.array([4.767, -7.543, -4.694]),
#         'CE1': np.array([5.148, -8.172, -5.816]),
#         'NE2': np.array([6.359, -8.738, -5.651])}],  
#     'A_128': ['CYS', {'SG': [2.247,-7.357,-2.372]}], 
#     'A_90': ['CYS', {'SG': [1.407, -8.539, -5.768]}]
# }

# dict_residues_in_candidate_motif = {
#     'A_55': ['ASP', {
#         'OD1': np.array([39.092, 19.909, 7.029]),
#         'OD2': np.array([41.131, 19.639, 6.336])
#     }],
#     'A_69': ['HIS', {
#         'CG': np.array([36.103, 22.677, 8.102]),
#         'CD2': np.array([35.042, 21.875, 7.997]),
#         'ND1': np.array([37.137, 21.804, 8.229]),
#         'CE1': np.array([36.760, 20.545, 8.200]),
#         'NE2': np.array([35.451, 20.570, 8.060])
#     }],
#     'A_122': ['ASP', {
#         'OD1': np.array([41.260, 24.248, 9.765]),
#         'OD2': np.array([39.156, 23.415, 9.817])
#     }],
#     'A_118': ['HIS', {
#         'CG': np.array([39.775, 19.425, 11.866]),
#         'CD2': np.array([39.236, 20.460, 11.168]),
#         'ND1': np.array([40.582, 18.720, 11.012]),
#         'CE1': np.array([40.552, 19.280, 9.808]),
#         'NE2': np.array([39.755, 20.299, 9.919])
#     }]}
# print (caclulate_dis_CoordinationAngles_HISangles_stats(dict_residues_in_candidate_motif))
# end_time=time.time()
# print (end_time-start_time)
