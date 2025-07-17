#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
csvs_to_sql.py

Generate a PostgreSQL-compatible SQL dump from four CSV tables. The script
includes a hardcoded SQL preamble and schema, followed by COPY data blocks
generated from the CSV files.

Interface:
    write_full_dump(csv_dir: str, output_path: str) -> None

CLI (defaults assume you’re in the project root):
    python csvs_to_sql.py

    -- Uses `editable_template_files_sql/` for the CSVs,
    -- Emits `dump.sql` in the root.

To override:
    python csvs_to_sql.py \
      --csv-dir path/to/csvs \
      --output out.sql
"""
import os
import csv
import argparse
import sys
from pathlib import Path

# Hardcoded SQL preamble and table schema
SQL_PREAMBLE_AND_SCHEMA = """--
-- PostgreSQL database dump
--

-- Dumped from database version 15.7 (Ubuntu 15.7-1.pgdg22.04+1)
-- Dumped by pg_dump version 16.3 (Ubuntu 16.3-1.pgdg22.04+1)

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_table_access_method = heap;

SET SEARCH_PATH = public;

DROP TABLE IF EXISTS minimized_training_cluster_information;
DROP TABLE IF EXISTS motif_representative_coordinates_table_v2;
DROP TABLE IF EXISTS motif_representative_detailed_coordinates_table_v2;
DROP TABLE IF EXISTS training_representative_metal_sites_kruskal_v2;


CREATE TABLE minimized_training_cluster_information (
    site_id integer,
    pdbid character varying(255),
    resi_comb character varying(255),
    num_of_residues integer,
    cluster_size integer
);

CREATE TABLE motif_representative_coordinates_table_v2 (
    pdbid character varying(255),
    chain_resi character varying(255),
    close_atom_coord real[],
    far_atom_coord real[],
    b_factor smallint
);

CREATE TABLE motif_representative_detailed_coordinates_table_v2 (
    pdbid character varying(255),
    chain_resi character varying(255),
    resi_type character varying(10),
    dict_atom_coord json
);

CREATE TABLE training_representative_metal_sites_kruskal_v2 (
    site_id integer,
    ii_order_id integer,
    pdbid character varying(255),
    aa_pair character varying(255),
    close_atom_dis double precision,
    far_atom_dis double precision,
    chain_resi_1 character varying(255),
    chain_resi_2 character varying(255)
);
"""

def generate_copy_blocks(csv_dir: str) -> str:
    """Generate PostgreSQL COPY blocks from CSV files in a directory."""
    tables = [
        'minimized_training_cluster_information',
        'training_representative_metal_sites_kruskal_v2',
        'motif_representative_coordinates_table_v2',
        'motif_representative_detailed_coordinates_table_v2'
    ]
    parts = []
    for tbl in tables:
        csv_path = Path(__file__).parent /"editable_template_files_sql"/f"{tbl}.csv"
        with open(csv_path, newline='', encoding='utf-8') as f:
            reader = csv.reader(f)
            headers = next(reader)
            rows = list(reader)
        parts.append(f"COPY {tbl} ({', '.join(headers)}) FROM stdin;")
        for row in rows:
            parts.append('\t'.join(row))
        parts.append('\\.')
        parts.append('')
    return '\n'.join(parts)


def write_full_dump() :
    csv_path ="editable_template_files_sql"
    output_path=Path(__file__).parent.parent/"setup_pg_db_with_tables"/"PostgreSQL_4_necessary_tables.sql"
    """Writes the complete SQL dump file."""
    # Get hardcoded preamble and CREATE statements
    header = SQL_PREAMBLE_AND_SCHEMA
    # Generate COPY data blocks
    csv_dir="editable_template_files_sql"
    copy_data = generate_copy_blocks(csv_dir)
    # Combine and write the final dump file
    dump = header + '\n\n' + copy_data + '\n-- PostgreSQL database dump complete\n'
    with open(output_path, 'w', encoding='utf-8') as out:
        out.write(dump)
    print(f"Full SQL dump written to {output_path}")
