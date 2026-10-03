#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source project.conf
bash scripts/stage.sh all
