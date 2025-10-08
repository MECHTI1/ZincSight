#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 00:30:56 2024

@author: mechti
"""
import os
import time
import requests
from concurrent.futures import ThreadPoolExecutor
# from src.settings import QUERY_STRUCTURES_DIR

# Generalized download function with file format handling
def download_files(url_template, ids, directory, file_extension):
    start_time=time.time()

    def download_file(file_id):
        url = url_template.format(file_id, file_extension)
        response = requests.get(url)
        if response.status_code == 200:
            filename = os.path.join(directory, f'{file_id}.{file_extension}')
            with open(filename, 'wb') as file:
                file.write(response.content)
            print(f'Downloaded {filename}')
        else:
            print(f'Failed to download {file_id}.{file_extension}')

    with ThreadPoolExecutor(max_workers=3) as executor:
        executor.map(download_file, ids)

    print("Overall Running Time:", time.time() - start_time)


# Specific download functions for AF, PDB (in CIF), and esm (in PDB)
def download_structures_af(uniprot_accessions, directory):
    af_ids = ["AF-" + accession + "-2-F1-model_v6" for accession in uniprot_accessions]
    url_template = 'https://alphafold.ebi.ac.uk/files/{}.{}'

    download_files(url_template, af_ids, directory, file_extension="cif")


def download_structures_pdb(pdb_ids, directory):
    url_template = 'https://files.rcsb.org/download/{}.{}'
    download_files(url_template, pdb_ids, directory, file_extension="cif")


def download_structures_esm(esm_ids, directory):
    url_template = 'https://api.esmatlas.com/fetchPredictedStructure/{}.{}'  # Replace with actual URL for esm PDBs
    download_files(url_template, esm_ids, directory, file_extension="pdb")


def main(uniprot_accessions, pdb_ids, esm_ids, path_query_structures):
    # AlphaFold and PDB downloads use CIF format
    download_structures_af(uniprot_accessions, path_query_structures)
    download_structures_pdb(pdb_ids, path_query_structures)

    # esm downloads use PDB format
    download_structures_esm(esm_ids, path_query_structures)


if __name__ == "__main__":
    uniprot_accessions = ['A0A2K5XT84', 'G3QSU8', 'A5Z1T7']
    pdb_ids = ['1CRN', '4HHB', '5XNL']
    esm_ids = ['MGYP002537940442', 'MGYP001215146166', 'MGYP001823580159']
    main(uniprot_accessions, pdb_ids, esm_ids)
