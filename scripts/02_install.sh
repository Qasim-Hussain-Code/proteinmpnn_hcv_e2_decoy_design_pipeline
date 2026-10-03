#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
source project.conf
mkdir -p vendor .cache
export PIP_CACHE_DIR="$PWD/.cache/pip"
export TMPDIR="$PWD/.cache/tmp"
mkdir -p "$TMPDIR"
"$PYTHON" -m pip install --no-cache-dir -r requirements.lock.txt --extra-index-url https://download.pytorch.org/whl/cpu
if [[ ! -d vendor/ProteinMPNN/.git ]]; then
    git clone --filter=blob:none --no-checkout https://github.com/dauparas/ProteinMPNN.git vendor/ProteinMPNN
    git -C vendor/ProteinMPNN sparse-checkout init --no-cone
    git -C vendor/ProteinMPNN sparse-checkout set protein_mpnn_run.py protein_mpnn_utils.py LICENSE README.md vanilla_model_weights/v_48_020.pt soluble_model_weights/v_48_020.pt
fi
mpnn_commit=$("$PYTHON" -c 'import json; print(json.load(open("config/upstream_lock.json"))["ProteinMPNN"])')
evo_commit=$("$PYTHON" -c 'import json; print(json.load(open("config/upstream_lock.json"))["EvoEF2"])')
git -C vendor/ProteinMPNN fetch origin "$mpnn_commit"
git -C vendor/ProteinMPNN checkout --detach "$mpnn_commit"
if [[ ! -d vendor/EvoEF2/.git ]]; then
    git clone https://github.com/tommyhuangthu/EvoEF2.git vendor/EvoEF2
fi
git -C vendor/EvoEF2 fetch origin "$evo_commit"
git -C vendor/EvoEF2 checkout --detach "$evo_commit"
if [[ "${OS:-}" == Windows_NT ]]; then
    g++ -O3 -static -o vendor/EvoEF2/EvoEF2_local.exe vendor/EvoEF2/src/*.cpp
else
    g++ -O3 -o vendor/EvoEF2/EvoEF2 vendor/EvoEF2/src/*.cpp
fi
"$PYTHON" -c 'import torch; assert torch.version.cuda is None, "CUDA build refused"'
