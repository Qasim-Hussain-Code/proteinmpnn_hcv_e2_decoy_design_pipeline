#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source project.conf
"$PYTHON" -m pytest -q
if command -v shellcheck >/dev/null 2>&1; then
    shellcheck scripts/*.sh run_all.sh
else
    echo "shellcheck unavailable; recorded in verification rather than silently passed"
fi
"$PYTHON" -m pipeline.cli verify
