#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == --help || "${1:-}" == -h || $# == 0 ]]; then
    echo 'Usage: bash scripts/stage.sh STAGE [CLI options]; use configure, sources, structures, humanize, ground_truth, sequences, diversity, annotations, states, benchmark, pilot, design, discovery, freeze, heldout, sensitivity, figures, report, verify or smoke'; exit 0
fi
source project.conf
export OMP_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export CUDA_VISIBLE_DEVICES=""
"$PYTHON" -m pipeline.cli "$@"
