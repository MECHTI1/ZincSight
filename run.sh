#!/usr/bin/env bash
set -euo pipefail
cd /app
exec python /app/integrations/tamarind/predict.py
