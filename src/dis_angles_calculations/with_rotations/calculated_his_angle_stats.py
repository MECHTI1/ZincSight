#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 13:02:33 2023

@author: mechti
"""

import numpy as np
import time
#ce1 instead ne2, cd2 instead nd1
def prepare_A_B_C(residues_dict):
    residues = []
    for res in residues_dict.values():
        if res[-1] == 'NE2':
            residues.append({
                "A": res[1]['CE1'],
                "B": res[1]['NE2'],
                "C": res[1]['CD2']
            })
        elif res[-1] == 'ND1':
            residues.append({
                "A": res[1]['CE1'],
                "B": res[1]['ND1'],
                "C": res[1]['CG']
            })
        elif res[-1] == 'CD2':
             residues.append({
                 "A": res[1]['NE2'],
                 "B": res[1]['CD2'],
                 "C": res[1]['CG']   
                }) 
        elif res[-1] == 'CE1':
                 residues.append({
                "A": res[1]['NE2'],
                "B": res[1]['CE1'],
                "C": res[1]['ND1']   
                }) 
                 
    return residues


def calculate_angles(mid_P, A, B, C):
    def angle_between(v1, v2):
        """Returns the angle in degrees between vectors 'v1' and 'v2'"""
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))  

    AC_vector = A - C
    D = (A + C) / 2
    DP_vector = mid_P - D
    angle_DP_AC = angle_between(DP_vector, -AC_vector)

    normal = np.cross(B - A, C - A)
    normal /= np.linalg.norm(normal)
    angle_DP_plane = 90 - angle_between(DP_vector, normal)

    return angle_DP_AC, angle_DP_plane



def calculate_difference_stats(values, target=0):
    deviations = np.array(values) - target
    
    mse = np.mean(deviations ** 2)
    rmse = np.sqrt(mse)
    # rmse=mse

    return rmse

def calculations_differences(candidate_point_angles):

    dif_angle_base = calculate_difference_stats (candidate_point_angles["angle_to_base"], target=90)
    dif_angle_plane = calculate_difference_stats (candidate_point_angles["angle_to_plane"], target=0)
    
    return {"candidate_point_angles_stast": {"dif_angle_base": dif_angle_base, "dif_angle_plane": dif_angle_plane}}
    

def check_candidate_metal_coord_valid_histidines_angles(mid_P, dict_residues_in_candidate_motif):

    residues = prepare_A_B_C(dict_residues_in_candidate_motif)
    
    angle_to_HIS_plane=[]
    angle_to_middle_point_base_vector=[]
    # calculate angles for each residue
    for i, residue in enumerate(residues):
        angle_DP_AC, angle_DP_plane = calculate_angles(mid_P, **residue)
        
        angle_to_middle_point_base_vector.append (angle_DP_AC)
        angle_to_HIS_plane.append (angle_DP_plane)

        # print(f"For residue {i + 1}:")
        # print(f"The angle at D between AC and the vector to the midpoint is {angle_DP_AC} degrees")
        # print(f"The angle between the vector to the midpoint and the plane formed by ABC is {angle_DP_plane} degrees")
        # print()
    
    # print (residues)
    # print ("angle_to_middle_point_base_vector", angle_to_middle_point_base_vector)
    # print ("angle_to_HIS_plane", angle_to_HIS_plane)
    # Return dictionary
    candidate_point_dict_angles = {
        "angle_to_base": angle_to_middle_point_base_vector,
        "angle_to_plane": angle_to_HIS_plane
    }
    
    
    calculated_his_angle_stats = calculations_differences(candidate_point_dict_angles)
    return calculated_his_angle_stats


# dict_residues_in_candidate_motif = {
#     'A_119': ['HIS', 
#               {'NE2': np.array([-0.891,  4.314, -1.945]), 
#                'ND1': np.array([-0.066,  2.318, -2.398]), 
#                'CE1': np.array([-0.679,  3.399, -2.897]), 
#                'CD2': np.array([-0.368,  3.785, -0.772]), 
#                'CG': np.array([ 0.145,  2.548, -1.061])}, 'ND1'], 
#     'A_194': ['HIS', 
#               {'NE2': np.array([ 1.705, -0.313, -3.149]), 
#                'ND1': np.array([ 2.964, -1.36 , -1.68 ]), 
#                'CE1': np.array([ 1.708, -1.126, -2.087]), 
#                'CD2': np.array([ 3.045, -0.004, -3.422]), 
#                'CG': np.array([ 3.816, -0.669, -2.504])}, 'NE2'], 
#     'A_196': ['HIS', 
#               {'NE2': np.array([ 1.288,  2.155, -5.458]), 
#                'ND1': np.array([ 1.927,  3.357, -7.185]), 
#                'CE1': np.array([ 1.357,  2.211, -6.797]), 
#                'CD2': np.array([ 1.806,  3.37 , -4.992]), 
#                'CG': np.array([ 2.2  ,  4.107, -6.073])}, 'NE2']
# }


# #mid_P = np.array([ 1.67754912,  1.7677847,  -3.43110193])
# mid_P = np.array([ 0.30961529,  1.00821601, -3.99592961])

# print (check_candidate_metal_coord_valid_histidines_angles (mid_P, dict_residues_in_candidate_motif))


