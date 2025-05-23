psql -U postgres -tc "SELECT 1 FROM pg_database WHERE datname = 'zincsight_pipeline_db';" | grep -q 1 || psql -U postgres -c "CREATE DATABASE zincsight_pipeline_db;"
psql -U postgres -c "GRANT ALL PRIVILEGES ON DATABASE zincsight_pipeline_db TO postgres;"
psql -U postgres -c "ALTER USER postgres WITH PASSWORD 'postgres';"
psql -U postgres -d zincsight_pipeline_db -P pager=off -f "src/setup_pg_db_with_tables/PostgreSQL_4_necessary_tables.sql"
