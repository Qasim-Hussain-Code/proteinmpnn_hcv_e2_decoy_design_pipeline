#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
    echo 'Usage: bash scripts/finish.sh; resumes discovery through the completed report'; exit 0
fi
source project.conf
for stage in discovery freeze heldout sensitivity figures report; do
    bash scripts/stage.sh "$stage"
done
