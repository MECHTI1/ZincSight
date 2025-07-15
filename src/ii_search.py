

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:57:55 2023

@author: mechti
"""

import psycopg2
from src.settings import get_db_connection
def create_motif_search_table_temp_1(conn):

    cur = conn.cursor()

    # SQL command to create Primary_metal_binding_site_table
    MOTIF_SEARCH_TABLE_TEMP_creation_command = """
     CREATE TABLE IF NOT EXISTS AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_1(
        ROW_ID SERIAL PRIMARY KEY,
        PDBID_AlphaFoldModel VARCHAR(255),
        chain_resi_1 VARCHAR(255),
        chain_resi_2 VARCHAR(255),
        chain_resi_3 VARCHAR(255),
        chain_resi_4 VARCHAR(255),
        chain_resi_5 VARCHAR(255),
        chain_resi_6 VARCHAR(255),
        chain_resi_7 VARCHAR(255),
        chain_resi_8 VARCHAR(255),
        chain_resi_9 VARCHAR(255),
        chain_resi_10 VARCHAR(255)
     );
     """
    cur.execute(MOTIF_SEARCH_TABLE_TEMP_creation_command)
    conn.commit()     # Commit the transaction


def fill_search_table_with_template_metal_binding_site(conn, site_id):
    cursor = conn.cursor()

    # Get the pdbid for the given site_id
    cursor.execute(
        "SELECT pdbid FROM training_representative_metal_sites_kruskal_v2 WHERE site_id = %s LIMIT 1", (site_id,))
    pdbid_result = cursor.fetchone()
    print ("pdbid_result",pdbid_result)
    pdbid_of_site_id = pdbid_result[0]

    # Initialize a new row in search_table_temp with the pdbid
    cursor.execute(
        "INSERT INTO AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_1 (PDBID_AlphaFoldModel) VALUES (%s)", (pdbid_of_site_id,))
    conn.commit()

    # Loop over the rows for the given site_id, ordered by ii_order_id
    cursor.execute(
        "SELECT * FROM training_representative_metal_sites_kruskal_v2 WHERE site_id = %s ORDER BY ii_order_id", (site_id,))
    rows = cursor.fetchall()
    print(rows)
    columns_names_in_representative_metal_sites_kruskal = [
        'site_id', 'ii_order_id', 'pdbid', 'aa_pair', 'close_atom_dis', 'far_atom_dis', 'chain_resi_1', 'chain_resi_2']
    for row in rows:
        # Process both residues for each row
        for i in range(1, 3):  # Iterate from 1 to 2
            v_residue_column = 'chain_resi_' + str(i)
            v_residue = row[columns_names_in_representative_metal_sites_kruskal.index(
                v_residue_column)]
            if v_residue is not None:
                # Update the first null residue field in the row if the residue doesn't already exist in the row
                for j in range(1, 11):
                    v_residue_column = 'chain_resi_' + str(j)
                    cursor.execute("UPDATE AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_1 SET {} = %s WHERE {} IS NULL".format(
                        v_residue_column, v_residue_column), (v_residue, ))
                    conn.commit()
                    cursor.execute("SELECT {} FROM AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_1 WHERE PDBID_AlphaFoldModel = %s".format(
                        v_residue_column), (pdbid_of_site_id,))
                    result = cursor.fetchone()

                    if result[0] == v_residue:
                        break


def create_first_round_insert_query():
    insert_sql_query_first_round = """
    INSERT INTO AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_1 (PDBID_AlphaFoldModel, chain_resi_1, chain_resi_2)
    SELECT t2.pdbid_alphafoldmodel,t2.chain_resi_1, t2.chain_resi_2            
    FROM AF_DATASET_ii_table_v2 AS t2
    INNER JOIN training_representative_metal_sites_kruskal_v2 AS t1
    ON t2.aa_pair = t1.aa_pair AND t1.close_atom_dis BETWEEN 0 AND 7
    WHERE t1.site_id = %s AND t1.ii_order_id = 1
    UNION ALL
    SELECT t2.pdbid_alphafoldmodel,t2.chain_resi_2, t2.chain_resi_1            
    FROM AF_DATASET_ii_table_v2 AS t2
    INNER JOIN training_representative_metal_sites_kruskal_v2 AS t1
    ON t2.aa_pair = t1.aa_pair AND t1.close_atom_dis BETWEEN 0 AND 7
    WHERE t1.site_id = %s AND t1.ii_order_id = 1 AND SUBSTRING(t1.aa_pair, 1, 1) = SUBSTRING(t1.aa_pair, 2, 1);
    ;
    """
    return insert_sql_query_first_round


def find_common_resi_kruskal_II_column_header(conn, site_id, ii_order_id):
    cur = conn.cursor()

    # Write your query
    query = """
        SELECT 
            chain_resi
        FROM
            (
            SELECT chain_resi_1 AS chain_resi, ii_order_id FROM training_representative_metal_sites_kruskal_v2 WHERE site_id = %s
            UNION ALL
            SELECT chain_resi_2 AS chain_resi, ii_order_id FROM training_representative_metal_sites_kruskal_v2 WHERE site_id = %s
            ) AS AllRows
        GROUP BY 
            chain_resi
        HAVING 
            %s = ANY(ARRAY_AGG(ii_order_id)) AND MIN(ii_order_id) < %s;
    """
    
    # Execute the query
    cur.execute(query, (site_id, site_id, ii_order_id, ii_order_id))
   # Fetch the results
    common_resi_metal_binding_site_kruskal = cur.fetchone()[0]
    print("common_resi_metal_binding_site_kruskal",
          common_resi_metal_binding_site_kruskal)

    # find in what column header does the common residue ,of the Kruskal II,found in the MOTIF_SEARCH_TABLE_TEMP
    query = "SELECT * FROM AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version} LIMIT 1".format(
        previous_table_version=ii_order_id-1)
    cur.execute(query)
    first_row = cur.fetchone()
    print ("first_row", first_row)
    column_names = [desc[0] for desc in cur.description]
    column_header_with_value = None
    for column, value_in_column in zip(column_names, first_row):
        if value_in_column == common_resi_metal_binding_site_kruskal:
            column_header_with_value = column
            break

    print("column_header_with_value", column_header_with_value)

    return column_header_with_value


def retrieve_AApair_CloseAtomDis_FarAtomDis_Of_kruskal_II(conn, site_id, ii_order_id):
    cur = conn.cursor()

    # Insert the results of the select query into the MOTIF_SEARCH_TABLE_TEMP table
    Retrieve_AApair_CloseAtomDis_FarAtomDis_Of_II_with_defined_IIorderId = """
    SELECT aa_pair,close_atom_dis, far_atom_dis
    FROM training_representative_metal_sites_kruskal_v2
    WHERE site_id = %s AND ii_order_id = %s;
    """
    cur.execute(Retrieve_AApair_CloseAtomDis_FarAtomDis_Of_II_with_defined_IIorderId, ((
        site_id, ii_order_id)))

    list_with_AApair_CloseAtomDis_FarAtomDis_Of_kruskal_II = cur.fetchone()
    return (list_with_AApair_CloseAtomDis_FarAtomDis_Of_kruskal_II)


def assign_dis_ranges_for_II_search(close_atom_dis, far_atom_dis, tolerance_value_close_atom, tolerance_value_far_atom):
    lower_limit_close_atom_dis = close_atom_dis-tolerance_value_close_atom
    upper_limit_close_atom_dis = close_atom_dis+tolerance_value_close_atom
    lower_limit_far_atom_dis = far_atom_dis-tolerance_value_far_atom
    upper_limit_far_atom_dis = far_atom_dis+tolerance_value_far_atom

    return (lower_limit_close_atom_dis, upper_limit_close_atom_dis, lower_limit_far_atom_dis, upper_limit_far_atom_dis)


def create_string_with_all_columns_with_new_appended_resi(resi_num_header_skip, ii_table_chain_resi_index, ii_order_id):
    string_to_add_to_command = "AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.ROW_ID, AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.pdbid_alphafoldmodel".format(
        previous_table_version=ii_order_id-1)
    for i in range(1, 11):
        if i != resi_num_header_skip:
            string_to_add_to_command += ", AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.chain_resi_{i}".format(
                previous_table_version=ii_order_id-1, i=i)
        else:
            string_to_add_to_command += ", AF_DATASET_ii_table_v2.chain_resi_{} AS chain_resi_{}".format(
                ii_table_chain_resi_index, i)

    return (string_to_add_to_command)


def generate_condition_string_that_verify_that_inserted_resi_not_existed_in_row(ii_order_id, II_table_resi_num):

    string_to_add_to_command = ""
    for i in range(1, ii_order_id+1):
        string_to_add_to_command += " AND AF_DATASET_ii_table_v2.chain_resi_{II_table_resi_num} != AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.chain_resi_{i}".format(
            II_table_resi_num=II_table_resi_num, previous_table_version=ii_order_id-1, i=i)
        # MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.chain_resi_{i}".format(previous_table_version=ii_order_id-1,i=i
    return (string_to_add_to_command)


def return_sql_query_single_round_motif_search_not_first_round(string_1, string_2, common_resi_header, ii_order_id, string_condition_1_verify_resi_new_for_row, string_condition_2_verify_resi_new_for_row):

    query = """
    CREATE TABLE IF NOT EXISTS AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{current_table_version} AS 
    SELECT * FROM (
        (
            SELECT * FROM AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version} LIMIT 1
        )
        UNION
        (
            SELECT {string_1} 
            FROM AF_DATASET_ii_table_v2
            INNER JOIN AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}
            ON (
                AF_DATASET_ii_table_v2.pdbid_alphafoldmodel= AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.pdbid_alphafoldmodel
                AND AF_DATASET_ii_table_v2.chain_resi_2 = AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.{common_resi_header}
                AND AF_DATASET_ii_table_v2.aa_pair = %s
                AND AF_DATASET_ii_table_v2.close_atom_dis BETWEEN 0 AND 7
                AND AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.ROW_ID>1
                {string_condition_1_verify_resi_new_for_row}
                )   
        )
        UNION
        (
            SELECT {string_2} 
            FROM AF_DATASET_ii_table_v2
            INNER JOIN AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}
            ON (
                AF_DATASET_ii_table_v2.pdbid_alphafoldmodel= AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.pdbid_alphafoldmodel
                AND AF_DATASET_ii_table_v2.chain_resi_1 = AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.{common_resi_header}
                AND AF_DATASET_ii_table_v2.aa_pair = %s
                AND AF_DATASET_ii_table_v2.close_atom_dis BETWEEN 0 AND 7
                AND AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_table_version}.ROW_ID>1
                {string_condition_2_verify_resi_new_for_row}
            )
        )
    ) AS subquery
    ORDER BY ROW_ID;
    """.format(current_table_version=ii_order_id, previous_table_version=ii_order_id-1, string_1=string_1, string_2=string_2,  common_resi_header=common_resi_header, string_condition_1_verify_resi_new_for_row=string_condition_1_verify_resi_new_for_row, string_condition_2_verify_resi_new_for_row=string_condition_2_verify_resi_new_for_row)
    
    print (query)
    return query


def create_FINAL_MOTIF_SEARCH_TABLE(conn,ii_order_id):
    cur = conn.cursor()

    # Specify the number of chain_resi columns you want
    num_chain_resi =ii_order_id+1
    
    # Create the CREATE TABLE query dynamically
    sql_command_create_empty_FINAL_MOTIF_SEARCH_TABLE= "CREATE TABLE IF NOT EXISTS AF_DATASET_FINAL_MOTIF_SEARCH_TABLE_V2 (match_id SERIAL PRIMARY KEY, PDBID_AlphaFoldModel VARCHAR(255)"
    
    
    for i in range(1, num_chain_resi + 1):
        sql_command_create_empty_FINAL_MOTIF_SEARCH_TABLE += f", chain_resi_{i} VARCHAR(255)"
    
    sql_command_create_empty_FINAL_MOTIF_SEARCH_TABLE += ");"
    

     # Define the SQL command
    cur.execute( sql_command_create_empty_FINAL_MOTIF_SEARCH_TABLE)
    
    string_for_insert_command = "pdbid_alphafoldmodel"
    for i in range(1, ii_order_id + 2):
        string_for_insert_command += f", chain_resi_{i}"
   

    sql_command_insert_to_FINAL_MOTIF_SEARCH_TABLE = """
       INSERT INTO AF_DATASET_FINAL_MOTIF_SEARCH_TABLE_V2 ({string_for_insert_command})
    SELECT {string_for_insert_command}
    FROM (
        (SELECT 1 AS sort_order, {string_for_insert_command} 
        FROM AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{last_table_version} LIMIT 1)
        UNION 
        (SELECT 2 AS sort_order, {string_for_insert_command}
        FROM AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{last_table_version}
        WHERE ROW_ID > 1)
    ) AS alias_subquery
    ORDER BY sort_order, pdbid_alphafoldmodel;
   ;
""".format(string_for_insert_command=string_for_insert_command, last_table_version= ii_order_id)

    print(sql_command_insert_to_FINAL_MOTIF_SEARCH_TABLE)
  # Execute the SQL command
    cur.execute(sql_command_insert_to_FINAL_MOTIF_SEARCH_TABLE)
    conn.commit()    # Commit the transaction


def iterative_motif_search(conn, site_id, Boolean_whether_delete_all_motif_search_table_temp_after_iteration):
    cur = conn.cursor()

    fill_search_table_with_template_metal_binding_site(conn, site_id)

    tolerance_value_close_atom = 1.5
    tolerance_value_far_atom = 3

    # First round Insert the results of the select query into the MOTIF_SEARCH_TABLE_TEMP table
    first_round_insert_query = create_first_round_insert_query()

    # Execute the query with variables
    cur.execute(first_round_insert_query, (site_id,site_id))
    conn.commit()    # Commit the transaction

    # search with following II AND with specific residues that are correlated to the soeicifc resi in the metal temPlate motifs. cant have been written before.
    # iterate over single site_id to fins structure motif match
    cur.execute(
        "SELECT COUNT(*) FROM training_representative_metal_sites_kruskal_v2 WHERE site_id = %s", (site_id,))
    number_of_ii_order_ids_in_metal_site = cur.fetchone()[0]
    range_ii_order_ids = range(2, number_of_ii_order_ids_in_metal_site+1)

    for ii_order_id in range_ii_order_ids:

        # call function that returns common reside column header with one of previous kruskal II
        common_resi_header = find_common_resi_kruskal_II_column_header(
            conn, site_id, ii_order_id)

        # retrieve aa_pair,close_atom_dis, far_atom_dis Of II with defined IIorderId
        aa_pair, close_atom_dis, far_atom_dis = retrieve_AApair_CloseAtomDis_FarAtomDis_Of_kruskal_II(
            conn, site_id, ii_order_id)

        lower_limit_close_atom_dis, upper_limit_close_atom_dis, lower_limit_far_atom_dis, upper_limit_far_atom_dis = assign_dis_ranges_for_II_search(
            close_atom_dis, far_atom_dis, tolerance_value_close_atom, tolerance_value_far_atom)

        resi_num_header_skip = ii_order_id+1
        string_1 = create_string_with_all_columns_with_new_appended_resi(
            resi_num_header_skip, 1, ii_order_id)
        string_2 = create_string_with_all_columns_with_new_appended_resi(
            resi_num_header_skip, 2, ii_order_id)

        string_condition_1_verify_resi_new_for_row = generate_condition_string_that_verify_that_inserted_resi_not_existed_in_row(
            ii_order_id, 1)
        string_condition_2_verify_resi_new_for_row = generate_condition_string_that_verify_that_inserted_resi_not_existed_in_row(
            ii_order_id, 2)

        query_single_round_motif_search_not_first_round = return_sql_query_single_round_motif_search_not_first_round(
            string_1,  string_2, common_resi_header, ii_order_id, string_condition_1_verify_resi_new_for_row, string_condition_2_verify_resi_new_for_row)

        cur.execute(query_single_round_motif_search_not_first_round, (aa_pair, aa_pair))
        # Commit the transaction
        conn.commit()
        
        if Boolean_whether_delete_all_motif_search_table_temp_after_iteration:
            query_delete_previous_table = "DROP TABLE AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{previous_temp_table_version}".format(previous_temp_table_version=ii_order_id-1)
            cur.execute(query_delete_previous_table)
            conn.commit()
        
            
    create_FINAL_MOTIF_SEARCH_TABLE(conn, ii_order_id)
    
    if Boolean_whether_delete_all_motif_search_table_temp_after_iteration:
        query_delete_previous_table = "DROP TABLE AF_DATASET_MOTIF_SEARCH_TABLE_TEMP_{last_temp_table_version}".format(last_temp_table_version= ii_order_id)
        print ("last_command_delete_table", query_delete_previous_table)
        cur.execute(query_delete_previous_table)
        conn.commit()
        

def main(site_id,Boolean_whether_delete_all_motif_search_table_temp_after_iteration):
    conn = get_db_connection()
    create_motif_search_table_temp_1(conn)
    iterative_motif_search(conn, site_id, Boolean_whether_delete_all_motif_search_table_temp_after_iteration)
    conn.close()


if __name__ == '__main__':
    site_id = 2
    Boolean_whether_delete_all_motif_search_table_temp_after_iteration=False
    main(site_id,Boolean_whether_delete_all_motif_search_table_temp_after_iteration)
