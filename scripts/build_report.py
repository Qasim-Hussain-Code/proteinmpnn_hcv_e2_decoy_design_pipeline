"""Render the integrated report from completed, unchanged analysis tables."""
from pathlib import Path
from string import Template
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def snapshot():
    rows = read_tsv(ROOT/'followup/original_snapshot.tsv')
    assert any(r['original_path'] == 'results/candidate_freeze.tsv' for r in rows)
    for row in rows:
        actual = hashlib.sha256((ROOT/row['path']).read_bytes()).hexdigest()
        if actual != row['sha256']:
            raise RuntimeError('Original scientific record changed: '+row['path'])


def build():
    snapshot()
    sources = [
        'results/run_summary.json', 'results/benchmark_summary.json',
        'results/diversity_summary.json', 'results/heldout_summary.json',
        'results/model_comparison_effects.tsv', 'followup/ground_truth.tsv',
        'followup/benchmark_summary.tsv', 'followup/model_audit.tsv',
        'followup/control_contacts.tsv', 'followup/coverage_summary.json',
        'followup/balanced_panel_summary.json',
    ]
    load = lambda p: json.loads((ROOT / p).read_text(encoding='utf-8'))
    run = load(sources[0]); benchmark = load(sources[1])
    cohort = load(sources[2]); heldout = load(sources[3])
    coverage = load(sources[9]); balanced = load(sources[10])
    truth = read_tsv(ROOT / sources[5])
    summaries = read_tsv(ROOT / sources[6])
    models = read_tsv(ROOT / sources[7])
    fixed = {r['model']: r for r in summaries
             if r['mode'] == 'fixed_e2' and float(r['tolerance']) == .5}
    repacked = {r['model']: r for r in summaries
                if r['mode'] == 'repacked' and float(r['tolerance']) == .5}
    effects = {(r['design_arm'], r['metric']): r
               for r in read_tsv(ROOT / sources[4])}
    contacts = {int(r['cd81_position']): r for r in read_tsv(ROOT / sources[8])
                if r['model'] == 'archive'}
    included = [r for r in truth if r['included'] == 'true']
    metrics = dict(
        original_correct=benchmark['concordant_directional'],
        control_total=benchmark['directional_binding_controls'],
        baseline_correct=benchmark['majority_direction_baseline_correct'],
        benchmark_assay_rows=sum(r['followup_role'] != 'external_source_challenge' for r in included),
        included_rows=len(included), model_count=len(models),
        best_bound_correct=int(fixed['bound_ae']['original_single_correct']),
        second_bound_correct=int(fixed['bound_bh']['original_single_correct']),
        external_missed=int(fixed['bound_ae']['external_reduced_total']) - int(fixed['bound_ae']['external_reduced_correct']),
        generated=run['generated'], unique_passed=run['unique_passed'],
        frozen=heldout['frozen'], heldout_states=heldout['states'],
        heldout_accessions=heldout['sequences'], shared_candidates=heldout['shared_candidates'],
        strategy_candidates=heldout['single_state_candidates'],
        effect=heldout['worst_delta_effect_escape_minus_single'],
        ci_low=heldout['ci95'][0], ci_high=heldout['ci95'][1],
        rmsd=run['humanization_rmsd'], original_retrieved=cohort['retrieved'],
        original_included=cohort['included'], original_development=cohort['discovery'],
        original_states=cohort['discovery_states'],
        original_positions=len(cohort['interface_positions']),
        original_coverage=coverage['archival_24_state_coverage'],
        retrieved=coverage['retrieved'], eligible=coverage['included'],
        excluded=coverage['excluded'], unique_e2=coverage['unique_e2'],
        development=coverage['development'], future_reserve=coverage['future_reserve'],
        expanded_positions=len(coverage['expanded_interface_positions']),
        global_states=coverage['states_for_90_percent'], global_coverage=coverage['achieved_coverage'],
        balanced_states=balanced['states'], balanced_coverage=balanced['coverage'],
        balanced_minimum=balanced['minimum_genotype_coverage'],
        unknown_genotype=coverage['unknown_genotype'],
        modeled_positions=len(coverage['expanded_interface_positions']) - len(coverage['unmodeled_interface_positions']),
        min_contact163=float(contacts[163]['min_distance']),
        min_contact196=float(contacts[196]['min_distance']),
        wt_repulsion=float(fixed['archive']['wt_vdwrep']),
        peak_rss_gb=run['peak_rss_bytes']/1e9, peak_project_gb=run['peak_project_bytes']/1e9,
    )
    assert not benchmark['benchmark_passed']
    assert metrics['included_rows'] == 34 and metrics['benchmark_assay_rows'] == 28
    substitutions = {k: str(v) for k, v in metrics.items()}
    for k in ['effect', 'ci_low', 'ci_high', 'rmsd', 'min_contact163', 'min_contact196', 'wt_repulsion']:
        substitutions[k] = f'{metrics[k]:.2f}'
    for k in ['original_coverage', 'global_coverage', 'balanced_coverage', 'balanced_minimum']:
        substitutions[k] = f'{metrics[k]:.2%}'
    for k in ['peak_rss_gb', 'peak_project_gb']:
        substitutions[k] = f'{metrics[k]:.3f}'
    for arm in ['A', 'B']:
        r = effects[arm, 'mean_hydrophobic_sasa_fraction']
        substitutions['exposure_' + arm] = f"{float(r['soluble_minus_standard']):.4f}"
        substitutions['exposure_ci_' + arm] = f"[{float(r['ci95_low']):.4f}, {float(r['ci95_high']):.4f}]"
    substitutions['model_rows'] = '\n'.join(
        f"| `{r['model']}` | {r['clashes_after']} | {float(fixed[r['model']]['wt_vdwrep']):.2f} | "
        f"{repacked[r['model']]['original_single_correct']}/11 | "
        f"{fixed[r['model']]['original_single_correct']}/11 | "
        f"{fixed[r['model']]['external_reduced_correct']}/4 |" for r in models)
    template = ROOT / 'scripts/report_template.md'
    text = Template(template.read_text(encoding='utf-8')).substitute(substitutions)
    assert 'figures/figure_1.png' not in text and 'figures/figure_4.png' not in text
    (ROOT / 'README.md').write_text(text, encoding='utf-8', newline='\n')
    assets = {'docs_report.md': 'docs/report.md', 'supplements.md': 'followup/README.md',
              'figures_index.md': 'figures/README.md'}
    for name, target in assets.items():
        (ROOT/target).write_bytes((ROOT/'scripts/report_assets'/name).read_bytes())
    evidence = dict(metrics=metrics, source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                                  for p in sources},
                    template_sha256=hashlib.sha256(template.read_bytes()).hexdigest(),
                    presentation_asset_sha256={target: hashlib.sha256((ROOT/target).read_bytes()).hexdigest()
                                               for target in assets.values()},
                    report_sha256=hashlib.sha256(text.encode()).hexdigest(),
                    interpretation='Integrated reporting; scientific inputs, scores and candidate freeze unchanged')
    (ROOT/'publication/unified_report_provenance.json').write_text(
        json.dumps(evidence, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    # Only documentation and archive-location records are refreshed.
    manifest = ROOT / 'followup/artifact_manifest.tsv'
    rows = read_tsv(manifest)
    for r in rows:
        if r['path'] in {'followup/README.md', 'followup/original_snapshot.tsv'}:
            p = ROOT/r['path']; r['sha256'] = hashlib.sha256(p.read_bytes()).hexdigest(); r['bytes'] = p.stat().st_size
    write_tsv(manifest, rows)
    snapshot()
    print('Integrated report rendered; original scientific freeze and archive verified.')


if __name__ == '__main__':
    build()
