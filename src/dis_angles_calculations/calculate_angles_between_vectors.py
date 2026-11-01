#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 13:00:03 2023

@author: mechti
"""

import numpy as np
from itertools import combinations
import time

def compute_angles_between_vectors(lst_resi_binding_atoms_coordinates, mid_p):
    # Create vectors from each point to the midpoint
    vectors = [point - mid_p for point in lst_resi_binding_atoms_coordinates]

    # Initialize a list to hold all the calculated angles
    angles = []

    # Iterate over all possible combinations of two vectors
    for vector_1, vector_2 in combinations(vectors, 2):
        # Calculate the dot product of the two vectors
        dot_product = np.dot(vector_1, vector_2)

        # Calculate the lengths (norms) of the vectors
        norm_1 = np.linalg.norm(vector_1)
        norm_2 = np.linalg.norm(vector_2)

        # Compute the cosine of the angle between the vectors
        cosine_angle = dot_product / (norm_1 * norm_2)

        # Convert the cosine to an angle in degrees
        angle = np.arccos(cosine_angle) * (180.0 / np.pi)

        # Append the calculated angle to the list of angles
        angles.append(angle)
        
    # Return the list of all calculated angles
    return angles


# lst_resi_binding_atoms_coordinates = [np.array([-0.363, -0.177, -8.548]), np.array([2.408, 0.017, -10.478]), np.array([2.087, 2.208, -8.224])]
# mid_p = np.array([1.70909773, 0.16120613, -8.50297083])

# start_time = time.time()

# print(compute_angles_between_vectors(lst_resi_binding_atoms_coordinates, mid_p))

# end_time = time.time()

# print(f"Execution time: {end_time - start_time} seconds")
