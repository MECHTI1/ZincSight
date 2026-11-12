# Tamarind mounts the submitted repository at /app at runtime. Dependencies
# and the pinned prediction engine must therefore live outside /app.
FROM node:22-bookworm-slim AS node_runtime
FROM python:3.12-slim-bookworm

COPY --from=node_runtime /usr/local/bin/node /usr/local/bin/node
COPY --from=node_runtime /usr/local/lib/node_modules /usr/local/lib/node_modules
RUN ln -s /usr/local/lib/node_modules/npm/bin/npm-cli.js /usr/local/bin/npm \
    && apt-get update \
    && apt-get install -y --no-install-recommends build-essential ca-certificates curl libgl1 libglu1-mesa \
    && rm -rf /var/lib/apt/lists/*

# Bump revision and archive digest together, then rerun the offline checks.
ENV ZINCSIGHT_ENGINE_REVISION=8d663023e8e0d4b9b6fe850307652beaff3a0946
RUN curl --fail --location --retry 3 \
      "https://codeload.github.com/MECHTI1/ZincSight/tar.gz/${ZINCSIGHT_ENGINE_REVISION}" \
      -o /tmp/zincsight.tar.gz \
    && echo "7b54c485f7ba9afda9eb7ed2f11a979cfca5e4fcf40420c274184117e8078026  /tmp/zincsight.tar.gz" | sha256sum -c - \
    && mkdir -p /opt/zincsight \
    && tar -xzf /tmp/zincsight.tar.gz --strip-components=1 -C /opt/zincsight \
    && printf '%s\n' "$ZINCSIGHT_ENGINE_REVISION" > /opt/zincsight/SOURCE_REVISION \
    && rm /tmp/zincsight.tar.gz

ENV CC=gcc \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHON_DOTENV_DISABLED=1 \
    MPLBACKEND=Agg \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1
RUN python -m pip install --no-cache-dir \
      -r /opt/zincsight/plugins/zincsight-preview/skills/zincsight-preview/scripts/requirements-lock.txt \
    && npm ci --prefix /opt/zincsight/integrations/mcp/pglite --ignore-scripts --no-audit --no-fund \
    && python -c "import pymol, psycopg2, Bio, scipy, sklearn" \
    && python -m pip check

# Avoid reliance on loopback TCP in Tamarind's isolated runtime. Keep the pinned
# database and SQL intact; communicate through the job's private Unix socket.
RUN python - <<'PY'
from pathlib import Path
path = Path('/opt/zincsight/integrations/mcp/pglite/server.mjs')
source = path.read_text()
changes = {
    "import { readFile, writeFile, rename } from 'node:fs/promises';":
        "import { readFile, writeFile, rename } from 'node:fs/promises';\nimport { dirname, join } from 'node:path';",
    "const server = new PGLiteSocketServer({db, host: '127.0.0.1', port: 0, maxConnections: 64});":
        "const socketDir = dirname(readyPath);\nconst server = new PGLiteSocketServer({db, path: join(socketDir, '.s.PGSQL.5432'), maxConnections: 64});",
    "const uri = `postgresql://postgres@${server.getServerConn()}/postgres?sslmode=disable`;":
        "const uri = `postgresql://postgres@/postgres?host=${encodeURIComponent(socketDir)}&port=5432&sslmode=disable`;",
}
for old, new in changes.items():
    if source.count(old) != 1:
        raise RuntimeError('The pinned PGlite launcher changed; review the Unix-socket patch.')
    source = source.replace(old, new)
path.write_text(source)
PY

WORKDIR /app
CMD ["bash", "/app/run.sh"]
