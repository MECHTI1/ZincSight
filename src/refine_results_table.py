from src.settings import get_db_connection

def generate_binding_resi_update():
    # Define the residue mappings with their template IDs
    mappings = {
        "{C,C,C,C}": ["{C,C,C,C}", 1],
        "{D,D,H,H}": ["{D,H,H,D}", 2],
        "{E,H,H}": ["{E,H,H}", 3],
        "{C,C,H}": ["{C,C,H}", 4],
        "{H,H,H}": ["{H,H,H}", 5],
        "{D,E,E,H}": ["{D,E,H,E}", 6],
        "{E,E,E,H}": ["{E,E,H,E}", 7],
        "{D,D,H}": ["{D,H,D}", 8],
        "{D,E,H}": ["{D,H,E}", 9],
        "{C,C,C,H}": ["{C,C,C,H}", 10],
        "{D,C,H,H}": ["{C,D,H,H}", 11],
        "{D,C,C,H}": ["{C,D,C,H}", 12],
        "{D,D,E,H}": ["{E,H,D,D}", 13],
        "{C,C,C}": ["{C,C,C}", 14],
        "{D,H,H}": ["{D,H,H}", 15],
        "{E,H,H,Y}": ["{E,H,Y,H}", 16],
        "{C,H,H}": ["{C,H,H}", 17],
        "{C,E,H}": ["{C,E,H}", 18],
        "{C,C,H,H}": ["{C,C,H,H}", 19],
        "{H,H,H,Y}": ["{H,Y,H,H}", 20],
        "{E,E,H,H}": ["{E,H,E,H}", 21],
        "{D,H,H,H}": ["{D,H,H,H}", 22],
        "{N,D,D}": ["{D,N,D}", 23],
        "{D,D,E}": ["{D,D,E}", 24],
        "{E,E,H}": ["{E,E,H}", 25],
        "{D,C,C,C}": ["{C,D,C,C}", 26],
        "{D,C,C}": ["{C,C,D}", 27],
        "{C,C,E,H}": ["{E,H,C,C}", 28],
        "{D,E,H,H}": ["{D,H,E,H}", 29],
        "{E,E,E}": ["{E,E,E}", 30],
        "{D,C,H}": ["{C,D,H}", 31],
        "{E,H,H,H}": ["{H,H,H,E}", 32],
        "{D,D,D}": ["{D,D,D}", 33],
        "{D,E,E}": ["{D,E,E}", 34],
        "{C,H,H,H}": ["{C,H,H,H}", 35],
        "{C,C,E}": ["{C,E,C}", 36],
        "{N,D,H,H}": ["{D,N,H,H}", 37],
        "{H,H,H,H}": ["{H,H,H,H}", 38],
        "{D,D,H,S}": ["{D,H,D,S}", 39],
        "{D,D,E,E}": ["{D,E,E,D}", 40],
        "{D,D,H,T}": ["{E,T,D,H}", 41],
        "{C,C,C,E}": ["{C,E,C,C}", 42]
    }

    # Generate the single update command for both resi_comb and template_id
    case_statements_resi = [
        f"WHEN resi_comb = '{k}' THEN '{v[0]}'"
        for k, v in mappings.items()
    ]

    case_statements_template = [
        f"WHEN resi_comb = '{k}' THEN {v[1]}"
        for k, v in mappings.items()
    ]

    update_sql = f"""
    UPDATE final_compressed_table_with_scored_binding_sites
    SET 
        resi_comb = CASE 
            {' '.join(case_statements_resi)}
            ELSE resi_comb
        END,
        template_id = CASE
            {' '.join(case_statements_template)}
            ELSE template_id
        END;
    """

    return update_sql

def refine_table():
    conn = get_db_connection()
    cur = conn.cursor()

    # Rename columns if they exist in the table
    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN pdbid_alphafoldmodel TO structure_id
    """)

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN dif_angle_base TO angles_to_mn_rms
    """)

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN dif_angle_plane TO angles_to_planes_rms
    """)

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN metalcoord TO predicted_ion_pos
    """)

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN coordination_anlges TO coordination_angles
    """)

    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN rmsd_close_atom TO template_rmsd
    """)

    conn.commit()

    # Add template_id column if it doesn't exist
    cur.execute("ALTER TABLE final_compressed_table_with_scored_binding_sites ADD COLUMN IF NOT EXISTS template_id INT;")
    conn.commit()

    # Get and execute the update SQL
    sql_update = generate_binding_resi_update()
    cur.execute(sql_update)
    conn.commit()

    # Rename resi_comb to resi_comb_ordered
    cur.execute("""
        ALTER TABLE final_compressed_table_with_scored_binding_sites
        RENAME COLUMN resi_comb TO resi_comb_ordered
    """)
    conn.commit()

    cur.close()
    conn.close()

if __name__ == "__main__":
    refine_table()