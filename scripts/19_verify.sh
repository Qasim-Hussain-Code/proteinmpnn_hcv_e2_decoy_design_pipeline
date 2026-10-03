#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
    echo 'Usage: bash scripts/19_verify.sh; runs Python tests, Bash static analysis and scientific provenance checks'; exit 0
fi
source project.conf
"$PYTHON" -m pytest -q
if command -v shellcheck >/dev/null 2>&1; then
    shellcheck scripts/*.sh run_all.sh
elif [[ -f vendor/shellcheck/shellcheck.exe ]]; then
    vendor/shellcheck/shellcheck.exe -x scripts/*.sh run_all.sh
else
    echo "shellcheck unavailable; recorded in verification rather than silently passed"
fi
"$PYTHON" -m pipeline.cli verify
