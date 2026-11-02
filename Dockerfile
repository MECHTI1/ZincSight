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
ENV ZINCSIGHT_ENGINE_REVISION=7bf80c27fd7f0e8ce7b6b100e305a7bc18620a03
RUN curl --fail --location --retry 3 \
      "https://codeload.github.com/MECHTI1/ZincSight/tar.gz/${ZINCSIGHT_ENGINE_REVISION}" \
      -o /tmp/zincsight.tar.gz \
    && echo "703a51cc5bcb1bed28ac59510d97149b10293a61f7919b1703ed6458392df061  /tmp/zincsight.tar.gz" | sha256sum -c - \
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

WORKDIR /app
CMD ["bash", "/app/run.sh"]
