#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ -x .venv/Scripts/python.exe ]]; then
    python_bin=.venv/Scripts/python.exe
elif [[ -x .venv/bin/python ]]; then
    python_bin=.venv/bin/python
else
    python_bin="${PYTHON:-python3}"
fi
"$python_bin" -m pipeline.cli configure "$@"
