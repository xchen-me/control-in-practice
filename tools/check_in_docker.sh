#!/bin/sh
# Author: Xu Chen, University of Washington
# Run every demo, starter and test inside the course image.
# Docker Desktop cannot bind-mount cloud-synced folders (Google Drive, OneDrive, iCloud): mounting one can
# hang the Docker daemon. So we copy this folder to a local mirror first and mount the mirror.
#   usage: tools/check_in_docker.sh [filter]     e.g. tools/check_in_docker.sh 05_
set -e
SRC="$(cd "$(dirname "$0")/.." && pwd)"
MIRROR="${CIP_MIRROR:-$HOME/.cache/control-in-practice-mirror}"
mkdir -p "$MIRROR"
rsync -a --delete --exclude out --exclude __pycache__ --exclude .ipynb_checkpoints "$SRC/" "$MIRROR/"
docker build -q -t control-in-practice:2026 "$MIRROR" >/dev/null
docker run --rm -v "$MIRROR":/course control-in-practice:2026 \
    sh -c "python run_all.py $1 && pytest -q -p no:cacheprovider tests"
echo "figures: $MIRROR/out"
