# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 14:18:02 2023

@author: mechti
"""


from src.settings import get_db_connection
# create table of results


def create_checking_table_and_insert_matches(conn):
    cur = conn.cursor()
    cur.execute(
        """ALTER TABLE AF_DATASET_final_motif_search_table_V2 
      ADD COLUMN rmsd_overall FLOAT, 
      ADD COLUMN rmsd_close_atom FLOAT;
      """
    )

    cur.execute(
        """UPDATE AF_DATASET_final_motif_search_table_V2
           SET rmsd_overall = AF_DATASET_rmsd_of_matches_table_V2.rmsd_overall,
           rmsd_close_atom = AF_DATASET_rmsd_of_matches_table_V2.rmsd_close_atom
           FROM AF_DATASET_rmsd_of_matches_table_V2
           WHERE AF_DATASET_final_motif_search_table_V2.match_id = AF_DATASET_rmsd_of_matches_table_V2.match_id;
        """
    )

    conn.commit()

def remove_duplicates_leave_only_lowest_rmsd(conn, num_residues):
   
    # create a new cursor
   cur = conn.cursor()
   cur.execute("ALTER TABLE AF_DATASET_rmsd_of_matches_table_V2 DROP CONSTRAINT AF_DATASET_rmsd_of_matches_table_V2_match_id_fkey;")
   
   # Building the chain_resi_n sequence
   residues_sequence = ", ".join([f"chain_resi_{i+1}" for i in range(num_residues)])
   
   # SQL statement
   sql = f"""
   WITH sorted_table AS (
       SELECT 
           pdbid_alphafoldmodel, 
           (SELECT STRING_AGG(residue, '_' ORDER BY residue) FROM unnest(ARRAY[{residues_sequence}]) as residue) as sorted_residues,
           MIN(rmsd_close_atom) as min_rmsd
       FROM 
           AF_DATASET_final_motif_search_table_V2
       GROUP BY 
           pdbid_alphafoldmodel, sorted_residues
   ),
   duplicates AS (
       SELECT 
           t.match_id
       FROM 
           AF_DATASET_final_motif_search_table_V2 t
       JOIN 
           sorted_table s ON t.pdbid_alphafoldmodel = s.pdbid_alphafoldmodel
           AND (SELECT STRING_AGG(residue, '_' ORDER BY residue) FROM unnest(ARRAY[{residues_sequence}]) as residue) = s.sorted_residues
           AND t.rmsd_close_atom > s.min_rmsd
   )
   DELETE FROM 
       AF_DATASET_final_motif_search_table_V2
   WHERE 
       match_id IN (SELECT match_id FROM duplicates)
   """
   # execute the SQL command
   cur.execute(sql)

   # commit the changes
   conn.commit()
     
      
     
def get_number_of_resi_metal_binding_site(conn):        
     
    cur = conn.cursor()
 
    # get number of residues in metal binding site
    cur.execute("""
    SELECT *
    FROM AF_DATASET_coordinates_of_matches_table_V2    
    WHERE match_id = 1
    """)
    # Fetch all the results
    rows_original_metal_binding_site = cur.fetchall()
    num_residues_in_metal_binding_site= len (rows_original_metal_binding_site)  # get number of residues in metal binding site, thats in order to know how many rows to select in a batch (dont want to fetch partially matching sites)
    return num_residues_in_metal_binding_site
    

def main():
    conn = get_db_connection()
    
    
    cur = conn.cursor()
    
    
    num_residues= get_number_of_resi_metal_binding_site(conn)
    
    
    create_checking_table_and_insert_matches (conn)

    Boolean_check = True
    
    if Boolean_check == True:
    
    #     add_site_id_to_similiar_places(conn)
        remove_duplicates_leave_only_lowest_rmsd(conn, num_residues)
        
    Boolean_delete_tables_rmsd_of_matches_table_coordinates_of_matches_table_matchid_to_siteid_table_after_running= True
    if Boolean_delete_tables_rmsd_of_matches_table_coordinates_of_matches_table_matchid_to_siteid_table_after_running:
        cur.execute ("DROP TABLE AF_DATASET_Coordinates_of_Matches_Table_V2, AF_DATASET_rmsd_of_matches_table_V2")
        conn.commit()
    
    conn.close()


if __name__ == "__main__":
    main()
