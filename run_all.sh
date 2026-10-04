#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mode=core
from_stage=
configure_args=()
while (($#)); do
    case "$1" in
        --help|-h) echo 'Usage: bash run_all.sh [--mode smoke|core|full] [--from STAGE] [--threads N] [--ram GB] [--disk GB] [--seed N] [--parallel] [--yes]'; exit 0 ;;
        --mode) mode="${2:?Missing mode}"; shift 2 ;;
        --from) from_stage="${2:?Missing stage}"; shift 2 ;;
        --threads|--ram|--disk|--seed) configure_args+=("$1" "${2:?Missing value}"); shift 2 ;;
        --yes|--parallel) configure_args+=("$1"); shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$mode" in smoke|core|full) ;; *) echo "Invalid mode: $mode" >&2; exit 2 ;; esac
if [[ ! -f project.conf ]] || ((${#configure_args[@]})); then
    bash scripts/00_configure.sh "${configure_args[@]}"
fi
source project.conf
if [[ "$mode" == smoke ]]; then
    bash scripts/stage.sh smoke
    exit 0
fi
# Core and full use the pilot's largest approved compact CPU budget. The mode
# label cannot override the measured resource ceiling or an existing freeze.
mapfile -t stages < <("$PYTHON" -c 'from pipeline.cli import STAGES; print("\n".join(STAGES[:-1]))' | tr -d '\r')
if [[ -n "$from_stage" ]] && [[ ! " ${stages[*]} " == *" $from_stage "* ]]; then
    echo "Unknown resume stage: $from_stage" >&2; exit 2
fi
start=false
for stage in "${stages[@]}"; do
    if [[ -z "$from_stage" || "$stage" == "$from_stage" ]]; then start=true; fi
    if [[ "$start" == true ]]; then bash scripts/stage.sh "$stage"; fi
done
