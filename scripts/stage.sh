#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source project.conf
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export CUDA_VISIBLE_DEVICES=""
"$PYTHON" -m pipeline.cli "$@"
