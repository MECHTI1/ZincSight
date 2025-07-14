# -*- coding: utf-8 -*-
"""
Spyder Editor

Created on Tue May 23 18:03:13 2023

This is a temporary script file.
"""
from Bio.PDB import *
from itertools import combinations
import numpy as np
from functools import lru_cache
import time
import math
from src.settings import DEBUGGING

# Declare the function signature
aa_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'} 
metal_binding_residues= ['CYS', 'HIS', 'ASP', 'GLU', 'MET', 'TYR', 'SER', 'THR','ASN']
metal_have_to_be_binding_residues= ['CYS', 'HIS', 'ASP', 'GLU']




def create_all_relevant_atoms_from_residues_list(structure, PDBID):
    # Define the atoms to be included for each residue
    atom_dict = {
        'HIS': ['CG', 'CD2', 'ND1', 'CE1', 'NE2', 'CB'],
        'GLU': ['OE1', 'OE2', 'CD', 'CG'],
        'CYS': ['SG', 'CB'],
        'ASP': ['OD1', 'OD2', 'CG', 'CB'],
        'MET': ['SD', 'CG'],
        'THR': ['OG1', 'CG2', 'CA'],
        'TYR': ['OH', 'CZ'],
        'SER': ['OG', 'CB'],
        'ASN': ['ND2', 'OD1', 'CB']
    }

    # Initialize the list of residues
    residue_list = []

    for model in structure:
        for chain in model:
            for residue in chain:
                residue_name = residue.get_resname()
    
                # If the residue is in the atom_dict, create a new entry for it
                if residue_name in atom_dict:
                    residue_id = str(chain.get_id()) + "_" + str(residue.get_id()[1])
                    atoms = atom_dict[residue_name]
                    atom_coords = {}
    
                    # Add the coordinates for the specified atoms
                    for atom in atoms:
                        if atom in residue:
                            # Get coordinates and round them to 4 decimal places
                            coords = residue[atom].get_coord().tolist()
                            atom_coords[atom] = [round(coord, 4) for coord in coords]

                    # Add the residue to the list as a tuple
                    residue_list.append((PDBID, residue_id, residue_name, atom_coords))

    return residue_list


@lru_cache(maxsize=None)
def check_b_factor(residue,return_BinaryForAlphaFold_or_value):
    c_alpha_atom_b_factor= residue["CA"].get_bfactor()
    if return_BinaryForAlphaFold_or_value== "Binary":
        if c_alpha_atom_b_factor>=0:
            return True
        else:
            return None
    else:                        # if the return_BinaryForAlphaFold_or_value== "Value"
        return round (c_alpha_atom_b_factor)

@lru_cache(maxsize=None)
def get_residue_id(residue):
   residue_full_id=residue.get_full_id()
   Identifier_format_resi=(str(residue_full_id[2])+("_")+str(residue_full_id[3][1]))
   return Identifier_format_resi

@lru_cache(maxsize=None)
def get_atom_close_coordinates(residue):
    atom_dict = {
        'HIS': ('NE2', 'ND1'),  # AVG (NE2, ND1)
        'GLU': ('OE1', 'OE2'),
        'CYS': 'SG',
        'THR': 'OG1',
        'ASP': ('OD1','OD2'),
        'MET': 'SD',
        'TYR': 'OH',
        'SER': 'OG',
        'ASN': 'OD1'
    }

    try:
        residue_name = residue.get_resname()
        atom_key = atom_dict[residue_name]
        if isinstance(atom_key, tuple):  # If it's a tuple we need to average
            atom_coords = [residue[atom].get_coord() for atom in atom_key]
            return np.mean(atom_coords, axis=0)
        else:  # If it's a single atom
            return residue[atom_key].get_coord()
    except:
        return None


@lru_cache(maxsize=None)
def get_far_atom_coordinates(residue):
    atom_dict = {
        'HIS': 'CG',
        'GLU': 'CG',
        'CYS': 'CB',
        'THR': 'CA',
        'ASP': 'CB',
        'MET': 'CG',
        'TYR': 'CZ',
        'SER': 'CB',
        'ASN': 'CB'
    }

    try:
        residue_name = residue.get_resname()
        atom_key = atom_dict[residue_name]
        return residue[atom_key].get_coord()
    except:
        return None          


def calc_distance(coord1, coord2):
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    return round(math.sqrt(dx*dx + dy*dy + dz*dz),1)



        
def primary_residue_coordinates_table(PDB_ID, resi_relevant):
    primary_residue_coordinates_list = []
    for residue in resi_relevant:
        res_id_suitable_format = get_residue_id(residue)
        atom_close_coordinates = get_atom_close_coordinates(residue)
        atom_far_coordinates = get_far_atom_coordinates(residue)

        atom_close_rounded_coordinates = np.round(atom_close_coordinates, 4).tolist()
        atom_far_rounded_coordinates = np.round(atom_far_coordinates, 4).tolist()
        b_factor=check_b_factor(residue,"Value")
        primary_residue_coordinates_list.append((PDB_ID, res_id_suitable_format, atom_close_rounded_coordinates, atom_far_rounded_coordinates, b_factor))
    return primary_residue_coordinates_list


            
def create_inverted_index_from_a_pair(pair,PDBID):
#    try:
        residue1,residue2=pair[0], pair[1]
        
        
        #get_name
        residue1_full_name,residue2_full_name= residue1.get_resname(),residue2.get_resname()
        
        if residue1_full_name not in metal_have_to_be_binding_residues and residue2_full_name not in metal_have_to_be_binding_residues:
            return
        
        residue1_short_name,residue2_short_name= aa_dict[residue1_full_name],aa_dict[residue2_full_name] # The seq1 module prodice short name

        #sort/order residues by name alphabet for faster search in future
        if residue1_short_name> residue2_short_name:
            residue1,residue2=residue2,residue1
            residue1_short_name,residue2_short_name=residue2_short_name,residue1_short_name
         
        atom_close_coordinates_1 = get_atom_close_coordinates(residue1)
        atom_close_coordinates_2 = get_atom_close_coordinates(residue2) 
        

        # Call the function
        Distance_close_atom_pair = calc_distance(atom_close_coordinates_1, atom_close_coordinates_2)
        
        
        if Distance_close_atom_pair>7: return
        
        atom_far_coordinates_1 = get_far_atom_coordinates(residue1)
        atom_far_coordinates_2 = get_far_atom_coordinates(residue2) 
        
        Distance_far_atom_pair = calc_distance(atom_far_coordinates_1, atom_far_coordinates_2)


        res_id1_suitable_format = get_residue_id(residue1)
        res_id2_suitable_format = get_residue_id(residue2)
        
        Header_name= str(residue1_short_name)+str(residue2_short_name)
       
        inverted_str_to_return=(PDBID,Header_name, Distance_close_atom_pair,Distance_far_atom_pair, res_id1_suitable_format, res_id2_suitable_format)
       
        return inverted_str_to_return
  
    #except:return

def Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure(PDBpath,PDBID):
        TIME_START=time.time()

        primary_list_inverted_index_one_protein=[]                
        
        parser = MMCIFParser() 

        # structure = parser.get_structure("PHA-L", PDBID +".cif") # required download of file before to computer
        structure = parser.get_structure(PDBID, PDBpath)

        # Detect NMR and handle only first model
        method = structure.header.get('structure_method', '') or structure.header.get('_exptl.method', '')
        methods = method.upper().split(';')
        is_nmr_structure = any('NMR' in m for m in methods)

        if DEBUGGING:
            print(f"for structureID={PDBID} , is_nmr_structure = {is_nmr_structure}")           #  Print whether NMR only if debugging is enabled

        if is_nmr_structure:
            models = list(structure.get_models())
            model_count = len(models)
            if model_count > 1:
                structure = [models[0]]  # Use only the first model


        all_relevant_atoms_from_residues_list = create_all_relevant_atoms_from_residues_list(structure,PDBID)

        if isinstance(structure, list): # In case it was NMR structure with multiple models
            all_resi_in_structure = structure[0].get_residues()
        else:
            all_resi_in_structure = structure.get_residues()

        resi_relevant_AA_type = [residue for residue in all_resi_in_structure if residue.get_resname() in metal_binding_residues]
        
        resi_relevant = []
        if "AF-" in PDBID or "MGYP" in PDBID:
                for residue in resi_relevant_AA_type:
                    if check_b_factor(residue,"Binary") is not None:
                       resi_relevant.append(residue)

        else:
            for residue in resi_relevant_AA_type:
                if get_atom_close_coordinates(residue) is not None and get_far_atom_coordinates(residue) is not None and residue.has_id("CA"):  # append residue to resi_relevant list only if have all the coordinates
                       resi_relevant.append(residue)

        
        primary_list_coordinates_one_protein= primary_residue_coordinates_table(PDBID,resi_relevant)
        iter_objects= combinations(resi_relevant, 2)
        
        for pair in iter_objects:
            pair_inverted_index= create_inverted_index_from_a_pair (pair,PDBID)
            if pair_inverted_index!= None: #Distance_pair> 6

                primary_list_inverted_index_one_protein.append (pair_inverted_index)

        TIME_FINISH=time.time()
        print ("overall_time", TIME_FINISH-TIME_START  )
        return (primary_list_inverted_index_one_protein, primary_list_coordinates_one_protein,all_relevant_atoms_from_residues_list)
                
if __name__=="__main__":
    Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure("/media/mechti/Data/download_structures_api_rcsb/structures_downloaded_rcsb_api/1a0q.cif","1a0q",1)

    #Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure("/media/mechti/Data/download_structures_api_rcsb_v2/Downloaded_structures_until_23.5.23/4wxx.cif","4wxx")
    #Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure("/home/mechti/Downloads/AF-A0A1E1XE96-F1-model_v4.cif","AF-A0A1E1XE96-F1-model_v4")
    #Create_IIs_Coordinates_ListOfTuples_of_ProteinStructure("/home/mechti/Downloads/8dqf.cif","8dqf")
    #4hk6
