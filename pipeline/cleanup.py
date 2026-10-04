"""Remove only audited generated PDB files after completed frozen evaluation."""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import re
from .common import ROOT, MeasuredStage, now, read_tsv, sha256, write_json, write_tsv
from .selection import require_freeze
from .structures import read_structure, sequence

MANIFEST = ROOT/'results/disposable_model_manifest.tsv'
ALLOWED = ('receptor_models', 'evo', 'scoring')

def disposable_models():
    frozen = require_freeze()
    generated = {r['candidate_id']: r for r in read_tsv(ROOT/'results/generated_sequences.tsv.gz')}
    discovery = read_tsv(ROOT/'results/discovery_scores.tsv.gz')
    heldout = read_tsv(ROOT/'results/heldout_scores.tsv')
    states = {r['haplotype_id'] for r in discovery}
    completed = defaultdict(set)
    for row in discovery:
        if row['score_status'] != 'passed':
            raise RuntimeError('Cleanup refused: incomplete discovery scores')
        completed[row['candidate_id']].add(row['haplotype_id'])
    if any(completed[c] != states for c in completed):
        raise RuntimeError('Cleanup refused: missing discovery states')
    held_states = {r['haplotype_id'] for r in heldout}
    for candidate in ['wild_type', *frozen['candidate_ids']]:
        rows = [r for r in heldout if r['candidate_id'] == candidate]
        if {r['haplotype_id'] for r in rows if r['score_status'] == 'passed'} != held_states:
            raise RuntimeError('Cleanup refused: incomplete held-out evaluation')
        if candidate != 'wild_type' and not (ROOT/f'data/processed/frozen_candidates/{candidate}.pdb').is_file():
            raise RuntimeError('Cleanup refused: frozen model missing')
    scored_sequences = {generated[c]['sequence'] for c in completed if c != 'wild_type'}
    rows = []
    for folder in ALLOWED:
        parent = ROOT/'data/work'/folder
        for path in sorted(parent.rglob('*.pdb')) if parent.exists() else []:
            resolved = path.resolve()
            if path.is_symlink() or not resolved.is_relative_to(parent.resolve()):
                raise RuntimeError(f'Cleanup refused: path escapes known model directory: {path}')
            retained = ''
            if folder == 'receptor_models':
                candidate = path.stem
                if candidate == 'wild_type':
                    retained = 'data/processed/humanized_repaired.pdb'
                elif candidate in generated:
                    actual = sequence(read_structure(path)[0]['R'])
                    if actual != generated[candidate]['sequence'] or actual not in scored_sequences:
                        raise RuntimeError(f'Cleanup refused: candidate identity or score coverage mismatch: {candidate}')
                    if candidate in frozen['candidate_ids']:
                        retained = f'data/processed/frozen_candidates/{candidate}.pdb'
                else:
                    continue
                reason = 'Sequence verified; complete discovery scores retained; frozen selections copied exactly'
            elif path.name == 'score_complex.pdb':
                reason = 'Disposable recombination input overwritten for each scored state; compact state files retained'
            elif re.fullmatch(r'(decoy_\d{4}|wild_type)_monomer\.pdb', path.name):
                candidate = path.stem.removesuffix('_monomer')
                if candidate != 'wild_type' and candidate not in generated:
                    continue
                reason = 'Disposable chain extraction for completed monomer scoring; biophysics table retained'
            else:
                continue
            digest = sha256(path)
            if retained and digest != sha256(ROOT/retained):
                raise RuntimeError(f'Cleanup refused: retained model differs: {path}')
            rows.append(dict(path=path.relative_to(ROOT).as_posix(), size_bytes=path.stat().st_size,
                             sha256=digest, retained_model=retained, justification=reason))
    return rows

def cleanup():
    require_freeze()
    if not MANIFEST.exists():
        raise RuntimeError('Run python -m pipeline.cleanup --plan and review the model manifest first')
    rows = read_tsv(MANIFEST)
    removed = 0
    removed_bytes = 0
    for row in rows:
        path = ROOT/row['path']
        resolved = path.resolve()
        if path.suffix != '.pdb' or path.is_symlink() or not any(
                resolved.is_relative_to((ROOT/'data/work'/folder).resolve()) for folder in ALLOWED):
            raise RuntimeError(f'Cleanup refused: manifest path outside audited model directories: {path}')
        if not path.exists():
            continue
        if sha256(path) != row['sha256']:
            raise RuntimeError(f'Cleanup refused: model changed after audit: {path}')
        if row['retained_model'] and sha256(ROOT/row['retained_model']) != row['sha256']:
            raise RuntimeError('Cleanup refused: retained model changed after audit')
        path.unlink()
        removed += 1
        removed_bytes += int(row['size_bytes'])
    previous = json.loads((ROOT/'results/cleanup_summary.json').read_text()) if (ROOT/'results/cleanup_summary.json').exists() else {}
    write_json(ROOT/'results/cleanup_summary.json', dict(timestamp=now(),
        manifest_sha256=sha256(MANIFEST), audited_files=len(rows),
        deleted_files_this_run=removed, deleted_bytes_this_run=removed_bytes,
        deleted_files_total=previous.get('deleted_files_total', 0)+removed,
        deleted_bytes_total=previous.get('deleted_bytes_total', 0)+removed_bytes,
        retained='Frozen models, controls, compact states, checkpoints, original FASTA/NPZ, source snapshots and logs',
        method='Nonrecursive unlink of individually hashed generated PDB files; no directory deletion'))
    print(f'Removed {removed} audited disposable PDB files ({removed_bytes} bytes).', flush=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', action='store_true', help='Validate outputs and write a reviewable per-file manifest without deleting files')
    parser.add_argument('--apply', action='store_true', help='Delete only unchanged files in the reviewed manifest')
    args = parser.parse_args()
    if args.plan == args.apply:
        parser.error('Choose exactly one of --plan or --apply')
    with MeasuredStage('cleanup_plan' if args.plan else 'cleanup', 'python -m pipeline.cleanup'):
        if args.plan:
            rows = disposable_models()
            if rows:
                write_tsv(MANIFEST, rows)
            elif not MANIFEST.exists():
                raise RuntimeError('No verified disposable model files found; no manifest created')
            print(f'Audited {len(rows)} disposable PDB files ({sum(r["size_bytes"] for r in rows)} bytes). No deletion performed.')
        else:
            cleanup()

if __name__ == '__main__':
    main()
