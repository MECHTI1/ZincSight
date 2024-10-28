import psycopg2
from src.settings import get_db_connection
import csv
import os 
from src.settings import TABLES_DIR

def export_final_table_to_csv_file():
    conn = get_db_connection()

    # Query to select specific columns
    query = """
    SELECT 
        pdbid_alphafoldmodel, 
        chain_resi_1, 
        chain_resi_2, 
        chain_resi_3, 
        chain_resi_4, 
        chain_resi_5, 
        resi_comb, 
        rmsd_overall, 
        rmsd_close_atom, 
        distances_list, 
        dif_angle_base, 
        dif_angle_plane, 
        metalcoord, 
        score
    FROM final_compressed_table_with_scored_binding_sites
    """
    
    
    
    # Create a cursor
    cur = conn.cursor()
    
    # Execute the query
    cur.execute(query)
    
    # Fetch all rows from the executed query
    rows = cur.fetchall()
    
    # Get column names from the cursor
    column_names = [desc[0] for desc in cur.description]
    
    csv_file_path = os.path.join(TABLES_DIR,'table_predicted_zn_binding_sites.csv')     # Path to the CSV file (current directory)
    
    # Write the rows to a CSV file
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write the column names as the header
        writer.writerow(column_names)
        # Write the data rows
        writer.writerows(rows)
    
    # Close the cursor and connection
    cur.close()
    conn.close()
    
    print(f"Data exported successfully to {csv_file_path}")
