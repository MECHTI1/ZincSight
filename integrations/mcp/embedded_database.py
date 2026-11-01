"""Optional, job-local PGlite adapter; not a general PostgreSQL replacement.

The MCP worker admits one structure and owns one disposable database. Autocommit
avoids nested connections waiting on a transaction in PGlite's single session.
Any failed job is discarded entirely; databases are never shared across jobs.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import shutil
import subprocess
import time


@contextmanager
def embedded_database(root: Path, launcher: Path):
    node = shutil.which('node')
    if node is None:
        raise RuntimeError('The pglite backend requires Node.js 22 or newer.')
    ready = root / 'database-ready.json'
    dump = root / 'src/setup_pg_db_with_tables/PostgreSQL_4_necessary_tables.sql'
    with (root / 'database.log').open('w') as log:
        process = subprocess.Popen([node, str(launcher), str(dump), str(ready)],
                                   stdout=log, stderr=log)
        try:
            deadline = time.monotonic() + 30
            while not ready.exists():
                if process.poll() is not None:
                    raise RuntimeError('Embedded PostgreSQL failed to start; check the optional npm dependencies.')
                if time.monotonic() >= deadline:
                    raise RuntimeError('Embedded PostgreSQL startup exceeded 30 seconds.')
                time.sleep(.05)
            metadata = json.loads(ready.read_text())

            # Install before importing the prediction modules, which import this
            # function by name. The standard backend and source files stay intact.
            import psycopg2
            from src import settings
            original = settings.get_db_connection

            def connect():
                conn = psycopg2.connect(metadata['uri'], connect_timeout=5)
                conn.autocommit = True
                return conn

            settings.get_db_connection = connect
            try:
                yield metadata
            finally:
                settings.get_db_connection = original
        finally:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
