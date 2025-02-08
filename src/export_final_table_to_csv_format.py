from src.settings import get_db_connection
import csv
import os

def export_final_table_to_csv_file(path_output):
    conn = get_db_connection()

    # Query to select specific columns
    query = """
    SELECT 
        structure_id, 
        chain_resi_1, 
        chain_resi_2, 
        chain_resi_3, 
        chain_resi_4, 
        resi_comb_ordered, 
        template_rmsd, 
        template_id,
        angles_to_mn_rms,
        angles_to_planes_rms,
        dis_rmsd,
        predicted_ion_pos,
        score,
        prob
    FROM final_compressed_table_with_scored_binding_sites
    ORDER BY prob DESC
    """
    
    
    
    # Create a cursor
    cur = conn.cursor()
    
    # Execute the query
    cur.execute(query)
    
    # Fetch all rows from the executed query
    rows = cur.fetchall()
    
    # Get column names from the cursor
    column_names = [desc[0] for desc in cur.description]

    table_dir = os.path.join(path_output, 'table')  # Path to 'table' directory
    os.makedirs(table_dir, exist_ok=True)   # Create 'table' directory if it doesn't exist
    csv_file_path = os.path.join(table_dir, 'table_predicted_zn_binding_sites.csv')     # Full path for saving the final CSV file 'table_predicted_zn_binding_sites.csv'

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

if __name__ == '__main__':
    from src.settings import RESULTS_DIR
    path_output = RESULTS_DIR
    export_final_table_to_csv_file(path_output)