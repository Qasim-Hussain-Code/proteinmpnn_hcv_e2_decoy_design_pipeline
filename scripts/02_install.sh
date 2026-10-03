#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
    echo 'Usage: bash scripts/02_install.sh; installs locked CPU dependencies and compiles pinned EvoEF2'; exit 0
fi
source project.conf
mkdir -p vendor .cache
export PIP_CACHE_DIR="$PWD/.cache/pip"
export TMPDIR="$PWD/.cache/tmp"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"
"$PYTHON" -c 'from pipeline.common import check_resources; check_resources(2_000_000_000,1,"installation size")'
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
if "$PYTHON" -c 'from pipeline.common import ROOT,read_tsv,sha256; p=ROOT/"results/software_manifest.tsv"; rows=read_tsv(p) if p.exists() else []; matches=[r for r in rows if r["tool"]=="EvoEF2"]; assert matches and all((ROOT/r["model_or_weight_file"]).exists() and sha256(ROOT/r["model_or_weight_file"])==r["weight_hash"] for r in matches)' 2>/dev/null; then
    echo 'Verified cached EvoEF2 executable hash; compilation retained'
elif [[ "${OS:-}" == Windows_NT ]]; then
    g++ -O3 -static -o vendor/EvoEF2/EvoEF2_local.exe vendor/EvoEF2/src/*.cpp
else
    g++ -O3 -o vendor/EvoEF2/EvoEF2 vendor/EvoEF2/src/*.cpp
fi
"$PYTHON" -c 'import torch; assert torch.version.cuda is None, "CUDA build refused"'
