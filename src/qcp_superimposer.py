#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:35:09 2023

@author: mechti
"""
#There us maybe option to wirte the same coordiantion twice but then the superimposition will not be equal to cg and metal binfing atom!!!!!!!!!!!!!!!!!!!!!!!
from Bio.PDB.QCPSuperimposer import QCPSuperimposer 
import numpy as np
import time

time_start=time.time()


def calculate_rmsd(coords1, coords2):
  try: 
    qcp = QCPSuperimposer()    # Initialize QCPSuperimposer object

    qcp.set(coords1[0::2], coords2[0::2]) # Set the coordinates to be superimposed
    qcp.run() # Superimpose the coordinate sets
    rmsd_close_atom = qcp.get_rms() # Calculate the RMSD between the first atoms in each set
    
    # qcp.set(coords1[1::2], coords2[1::2]) # Set the coordinates to be superimposed
    # qcp.run() # Superimpose the coordinate sets
    # rmsd_far_atom = qcp.get_rms() # Calculate the RMSD between all other atoms
    
    qcp.set(coords1, coords2) # Set the coordinates to be superimposed
    qcp.run() # Superimpose the coordinate sets
    rmsd_overall = qcp.get_rms()# Calculate the RMSD between the first atoms in each set
       #    weighted_rmsd= (2*(rmsd_close_atom )+(rmsd_far_atom))/3
   
    return (rmsd_overall, rmsd_close_atom)
  except Exception as e:
      return ("error","error")



def main():
    time_start = time.time()

    # Set two sets of 4 coordinates each
    coords1 = np.array([
        [158.164, 155.680, 190.176],  # Close atom
        [157.405, 157.353, 191.378],  # Far atom
        [156.138, 162.207, 195.344],  # Close atom
        [157.260, 160.785, 194.121]   # Far atom
    ])

    coords2 = np.array([
        [24.912,   6.002,    6.883],  # Close atom
        [24.907,   4.228,    5.610],  # Far atom
        [25.307,   0.852,    2.713],  # Close atom
        [26.133,   2.689,    3.578]   # Far atom
    ])

    rmsd_overall, rmsd_close_atom = calculate_rmsd(coords1, coords2)
    
        
    print("rmsd_close_atom:", rmsd_close_atom)
    print("rmsd_overall:", rmsd_overall)
    
    # Print the results
    time_end = time.time()
    overall_time = time_end - time_start
    print("overall_time", overall_time)

if __name__ == "__main__":
    main()