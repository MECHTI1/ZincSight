#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 15:06:52 2024

@author: mechti
"""

import numpy as np
from itertools import product
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation as R
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import time



def rotate_point(point, axis, theta):
    rotation = R.from_rotvec(theta * axis / np.linalg.norm(axis))
    return rotation.apply(point)

def rotate_histidine_both_directions(residue_info,rotation_angle):
    # print (residue_info)
    atom_positions = residue_info[1]
    nd1 = atom_positions['NE2']
    ce1 = atom_positions['CE1']
    cg = atom_positions['CG']
    
    midpoint = (nd1 + ce1) / 2.0
    rotation_axis = cg - midpoint
    rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize the axis

    results = {
        'original': atom_positions.copy(),
        '15_deg_clockwise': {},
        '15_deg_counterclockwise': {}
    }
    
    theta_clockwise = np.radians(rotation_angle)
    for atom, position in atom_positions.items():
        results['15_deg_clockwise'][atom] = rotate_point(position - midpoint, rotation_axis, theta_clockwise) + midpoint

    theta_counterclockwise = -theta_clockwise
    for atom, position in atom_positions.items():
        results['15_deg_counterclockwise'][atom] = rotate_point(position - midpoint, rotation_axis, theta_counterclockwise) + midpoint
    
    return results


def generate_rotation_combinations(input_dict,rotation_angle):
    # Filter HIS residues and generate rotations
    his_residues = {key: value for key, value in input_dict.items() if value[0] == 'HIS'}
    
    all_rotations = {residue: rotate_histidine_both_directions(info,rotation_angle) for residue, info in his_residues.items()}
    
    # Generate all combinations of rotations
    keys, rotation_options = zip(*all_rotations.items())
    combinations = [dict(zip(keys, combo)) for combo in product(*rotation_options)]
    
    # Compile detailed combinations
    detailed_combinations = []
    for combo in combinations:
        detailed_combo = input_dict.copy()  # Copy the original input_dict
        for residue, rotation in combo.items():
            # Retain the residue name and update the atom positions
            detailed_combo[residue] = [input_dict[residue][0], all_rotations[residue][rotation]]
        detailed_combinations.append(detailed_combo)
    # print ("iterating_rounds")
    return detailed_combinations








if __name__=='__main__':
    
    def plot_atoms_with_equal_aspect(atom_positions, title="Atom positions"):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plotting atoms
        for atom, position in atom_positions.items():
            ax.scatter(position[0], position, position[2], label=atom)

        # Connecting each atom with lines
        atom_keys = list(atom_positions.keys())
        for i in range(len(atom_keys)):
            for j in range(i + 1, len(atom_keys)):
                atom1 = atom_keys[i]
                atom2 = atom_keys[j]
                ax.plot([atom_positions[atom1][0], atom_positions[atom2][0]], 
                        [atom_positions[atom1][1], atom_positions[atom2][1]], 
                        [atom_positions[atom1][2], atom_positions[atom2][2]], 'gray')

        # Setting the aspect ratio to be equal
        max_range = np.array([np.ptp([atom[coord] for atom in atom_positions.values()]) for coord in range(3)]).max() / 2.0
        mid_x = np.mean([atom[0] for atom in atom_positions.values()])
        mid_y = np.mean([atom[1] for atom in atom_positions.values()])
        mid_z = np.mean([atom[2] for atom in atom_positions.values()])

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        elev = 30  # Set the elevation angle
        azim = 45  # Set the azimuthal angle
        # Setting the view angle
        ax.view_init(elev=elev, azim=azim)
        ax.legend()
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')
        plt.title(title)
        plt.show()
    
    start_calculations_time=time.time()
    input_dict = {
            'A_69': ['HIS', {
            'CG': np.array([36.103, 22.677, 8.102]),
            'CD2': np.array([35.042, 21.875, 7.997]),
            'ND1': np.array([37.137, 21.804, 8.229]),
            'CE1': np.array([36.760, 20.545, 8.200]),
            'NE2': np.array([35.451, 20.570, 8.060])
        }],
        'A_122': ['ASP', {
            'OD1': np.array([41.260, 24.248, 9.765]),
            'OD2': np.array([39.156, 23.415, 9.817])
        }],
        'A_118': ['HIS', {
            'CG': np.array([39.775, 19.425, 11.866]),
            'CD2': np.array([39.236, 20.460, 11.168]),
            'ND1': np.array([40.582, 18.720, 11.012]),
            'CE1': np.array([40.552, 19.280, 9.808]),
            'NE2': np.array([39.755, 20.299, 9.919])
        }]} 
    rotation_angle=10
    combinations = generate_rotation_combinations(input_dict,rotation_angle)
    # for combo in combinations:
    
    # print (combinations)
    end_calculations_time=time.time()
    print (combinations)
    start_plot_time=time.time()
    # Test the function with the atom positions
    for i in range(9):
        plot_atoms_with_equal_aspect(
        combinations[i]['A_69'][1])

    plot_atoms_with_equal_aspect({
    'CG': np.array([39.775, 19.425, 11.866]),
    'CD2': np.array([39.236, 20.460, 11.168]),
    'ND1': np.array([40.582, 18.720, 11.012]),
    'CE1': np.array([40.552, 19.280, 9.808]),
    'NE2': np.array([39.755, 20.299, 9.919])
    })
    end_plot_time=time.time()
    
    print ("calculations_time: " , end_calculations_time-start_calculations_time)
    print ("plot_time: ", end_plot_time-start_plot_time)
