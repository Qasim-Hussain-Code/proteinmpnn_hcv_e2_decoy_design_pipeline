#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source project.conf
for stage in discovery freeze heldout sensitivity figures report; do
    bash scripts/stage.sh "$stage"
done
