#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main_add_templates.py

Orchestrate the three steps:
 1) Export CSVs from SQL dump
 2) Append the new motif
 3) Generate a full SQL dump (COPY format)

Works both as a script (CLI) and importable in notebooks.

Usage (CLI):
    python main_add_templates.py \
      --sql-dump origin_sql_file.sql \
      --cif-dir ../mmcif_with_ZN \
      --pdbid 5ROB \
      --residues B_8,B_9,B_19

Defaults:
  - csv output dir: editable_template_files_sql
  - schema SQL: origin_sql_file.sql
  - output full dump: dump.sql

Importable:
    from main_add_templates import main, create_csv_tables, add_template, write_dump
    main()
"""
import argparse
import os
from pathlib import Path

from src.inject_new_motif_into_templates.step1_create_csv_tables_from_the_originsql import parse_and_export
from src.inject_new_motif_into_templates.step2_add_summary_motif import main as add_template
from src.inject_new_motif_into_templates.step3_csv_to_sql import write_full_dump
import requests


def download_pdb_file(pdbid):
    response = requests.get(f'https://files.rcsb.org/download/{pdbid}.cif')
    response.raise_for_status()
    download_dir = Path(__file__).parent.parent.parent / 'cif_templates'
    download_dir.mkdir(exist_ok=True)
    pdb_path = download_dir / f'{pdbid}.cif'
    with open(pdb_path, 'wb') as file:
        file.write(response.content)
    print(f'Downloaded {pdb_path}')
    return pdb_path


def main(
    pdbid: str = '5ROB',
    residues: str = 'B_8,B_9,B_19'
) -> None:

    pdb_path=download_pdb_file(pdbid)

    # Step1: extract CSVs
    parse_and_export()

    # Step2: append new motif
    add_template(pdb_path,pdbid,residues)
    # Step3: generate SQL dump
    write_full_dump()
