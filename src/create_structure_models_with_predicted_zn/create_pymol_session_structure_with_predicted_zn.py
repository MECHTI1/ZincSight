#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:14:53 2024

@author: mechti
"""

# Import the PyMOL module
from pymol import cmd
from pathlib import Path
import os
##TODO:DETELEAFTER: from src.settings import STRUCTURES_WITH_PREDICTED_ZN

# Custom color dictionary based on score
import json

# Load the dictionary from the JSON file
with open('src/create_structure_models_with_predicted_zn/score_precision_dict.json', 'r') as file:
    threshold_precision_dict = json.load(file)


def get_precision(score):
    for score_thresold in sorted(threshold_precision_dict.keys(), key=float):
        if score <= float(score_thresold):
            precision = float(threshold_precision_dict[score_thresold])
            return precision
    return float(
        min(threshold_precision_dict.values()))  # Default to higher than highest threshold the precision is lowest


def create_pymol_session_structure_with_predicted_zn(path_structure, list_of_predicted_ZN, path_output):
    # Initialize PyMOL
    cmd.reinitialize()

    structure_name = Path(path_structure).stem
    cmd.load(path_structure, structure_name)  # Load an existing CIF file

    cmd.set_color("light_blue", [0.85, 0.85, 1])  # Define a custom soft blue-purple color for the protein
    cmd.select("protein_only", "polymer.protein and not resn HOH")  # Select the protein, excluding cofactors and water
    cmd.color("light_blue", "protein_only")  # Color the protein with the custom color

    if list_of_predicted_ZN is not None:
        # Create a selection for all predicted zinc atoms
        cmd.select("predicted_zincs", "none")

        for index, score_coordinates in enumerate(list_of_predicted_ZN):
            score, coordinates = score_coordinates[0], score_coordinates[1]

            # Add a zinc atom at the specified position
            [x, y, z] = coordinates
            selection_name = f"pZN_{index}_{round(score, 2)}"
            cmd.pseudoatom(object=selection_name, elem="Zn", pos=[x, y, z], name="Zn")

            # Add this zinc to the overall selection
            cmd.select("predicted_zincs", f"predicted_zincs or {selection_name}")

            # Change the representation to sphere
            cmd.show("spheres", selection_name)
            cmd.set("sphere_scale", 0.7, selection_name)

            # Set b-factor to the threshold value
            precision = get_precision(score)
            cmd.alter(selection_name, f"b={precision}")

            # Select and show as sticks and color the protein residues within 3 Ångström of the zinc atom
            nearby_residues_selection = f"nearby_residues_{index}"
            cmd.select(nearby_residues_selection,
                       f"byres ((name CG,CD2,ND1,CE1,NE2,OE1,OE2,SG,OD1,OD2,SD,OG1,CG2,OH,OG,ND2) within 3.0 of {selection_name})")
            cmd.set_color("dark_purple", [0.6, 0.4, 0.8])
            cmd.show("sticks", nearby_residues_selection)
            cmd.color("dark_purple", nearby_residues_selection)
            cmd.util.cnc(nearby_residues_selection)

        # Apply the color spectrum to all predicted zinc atoms
        thresholds = list(map(float, threshold_precision_dict.values()))

        # Create a spectrum using these custom colors
        cmd.spectrum("b", "red_red_orange_yellow_green", "predicted_zincs", minimum=min(thresholds),
                     maximum=max(thresholds))


    # Create the subdirectory for structures with predicted Zn
    path_structures_with_predicted_zn = os.path.join(path_output, "structures_with_predicted_zn")
    os.makedirs(path_structures_with_predicted_zn, exist_ok=True)    # Create a directory, do nothing if it already exists
    ##TODO:DETELEAFTER: os.makedirs(STRUCTURES_WITH_PREDICTED_ZN, exist_ok=True)
    cmd.delete("protein_only")  # delete selection

    if list_of_predicted_ZN is not None:
        with_or_no = "with"
    else:
        with_or_no = "no"

    # Save the entire session
    full_path_pymol_session=os.path.join(path_structures_with_predicted_zn, f"{structure_name}_{with_or_no}_predicted_ZN.pse")
    ##TODO:DETELEAFTER: full_path_pymol_session = os.path.join(STRUCTURES_WITH_PREDICTED_ZN,f"{structure_name}_{with_or_no}_predicted_ZN.pse")


    cmd.save(full_path_pymol_session)

    print(f"Session saved as '{full_path_pymol_session}'.")


if __name__ == "__main__":
    from src.settings import RESULTS_DIR
    path_output=RESULTS_DIR

    list_of_predicted_ZN = [(17, [-20.313, 5.506, 50.392]), (0.8, [-18.313, 5.506, 50.392]),
                            (22, [-23.313, 5.506, 50.392])]

    #TODO: Need to add path to structure
    create_pymol_session_structure_with_predicted_zn(path_structure, list_of_predicted_ZN, path_output)