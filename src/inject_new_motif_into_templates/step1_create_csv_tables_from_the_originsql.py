#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
export_sql_dump_to_csv.py

Parse a PostgreSQL dump SQL file containing COPY ... FROM stdin operations
and export each table's data into a separate CSV file with headers.
If run without arguments, prompts interactively for paths.
"""

import argparse
import os
import csv
import re
from pathlib import Path

def parse_and_export():
    output_dir= Path(__file__).parent /"editable_template_files_sql"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    sql_path= Path(__file__).parent / "origin_sql_file.sql"
    copy_re = re.compile(
        r"^COPY\s+([\w\.\"]+)\s*\(([^)]+)\)\s+FROM\s+stdin;", re.IGNORECASE
    )
    current = None
    columns = []
    writer = None
    csv_file = None

    with open(sql_path, 'r', encoding='utf-8') as f:
        for line in f:
            if current is None:
                m = copy_re.match(line)
                if m:
                    full_table = m.group(1)
                    table_name = full_table.split('.')[-1].strip('"')
                    columns = [c.strip() for c in m.group(2).split(',')]
                    csv_path = os.path.join(output_dir, f"{table_name}.csv")
                    csv_file = open(csv_path, 'w', newline='', encoding='utf-8')
                    writer = csv.writer(csv_file)
                    writer.writerow(columns)
                    current = (table_name, columns)
                continue

            if line.strip() == '\\.':
                csv_file.close()
                current = None
                writer = None
                csv_file = None
                continue

            row = line.rstrip('\n').split('\t')
            writer.writerow(row)

    print(f"Exported tables to CSV files in '{output_dir}'")