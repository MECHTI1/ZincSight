def debug_print_last_table(conn, table_name):
    cur = conn.cursor()
    cur.execute(f"SELECT * FROM {table_name};")
    rows = cur.fetchall()
    colnames = [desc[0] for desc in cur.description]
    print(" | ".join(colnames))
    print("-" * 50)
    for row in rows:
        print(" | ".join(str(val) for val in row))
