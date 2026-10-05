# Supplementary tables and reproduction

The [integrated study report](../README.md) contains the scientific analysis and interpretation. This directory holds the structural diagnostics, corrected evidence and expanded coverage artifacts used in that report. Its historical `followup` name is retained so that recorded file paths and scripts continue to resolve.

## Evidence and structural diagnostics

| Artifact | Contents |
|---|---|
| [Corrected ground truth](ground_truth.tsv) | Assay-specific labels, construct identities, exclusions and source provenance |
| [Curation audit](curation_audit.tsv) | Correction of the soluble double-mutant attribution |
| [Locked diagnostic protocol](protocol.json) | Structural comparisons and prespecified choices |
| [Model audit](model_audit.tsv) | Template identity, fit atoms, missing residues, contacts and clashes |
| [Benchmark summary](benchmark_summary.tsv) | Control recovery by model, preparation mode and declared tolerance |
| [Mutation scores](benchmark.tsv) | Mutant-minus-WT score changes and predictions |
| [Energy terms](energy_terms.tsv) | Individual EvoEF2 contributions |
| [Control contacts](control_contacts.tsv) | Distances relevant to interpreting mutation effects |
| [WT preparation sensitivity](wt_repair_sensitivity.tsv) | Matched additional WT repair checks |

## Natural-sequence coverage

| Artifact | Contents |
|---|---|
| [Sequence QC](sequence_manifest.tsv) | Included and excluded accessions, with reasons |
| [Retrieval records](downloads.tsv) | Queries, identifiers, download dates and hashes |
| [Interface positions](interface_positions.tsv) | Expanded position definition and structural mapping |
| [Coverage summary](coverage_summary.json) | Cohort counts, limitations and reserve definition |
| [Balanced development panel](balanced_development_panel.tsv) | Final development haplotypes and selection provenance |
| [Balanced panel summary](balanced_panel_summary.json) | Overall and minimum genotype coverage |
| [Per-genotype coverage](balanced_genotype_coverage.tsv) | Denominators and represented accessions |
| [State manifest](state_manifest.tsv) | Expected coordinate hashes, mutations and missing positions |

## Scientific reproduction

Use the pinned CPU environment, verified original source snapshots and software described in the [software manifest](../results/software_manifest.tsv). Exact snapshot hashes are required to reproduce this completed run; a fresh live query may return different records. The published checkout retains its scientific freeze and refuses to silently replace it. A separately initiated analysis requires its own configuration, source audit and freeze.

The original design/evaluation workflow uses Git Bash:

```bash
bash scripts/00_configure.sh --threads 2 --ram 14 --disk 13 --seed 20261004 --yes
bash scripts/02_install.sh
bash run_all.sh --mode core --parallel
bash scripts/19_verify.sh
```

With the original source/software caches available, the diagnostic and coverage stages run from the project root:

```powershell
.\.venv\Scripts\python.exe -m pipeline.followup prepare
.\.venv\Scripts\python.exe -m pipeline.followup models
.\.venv\Scripts\python.exe -m pipeline.followup sequences
.\.venv\Scripts\python.exe -m pipeline.followup coverage
.\.venv\Scripts\python.exe -m pipeline.followup benchmark
.\.venv\Scripts\python.exe -m pipeline.followup states
.\.venv\Scripts\python.exe -m pipeline.followup controls
.\.venv\Scripts\python.exe -m pipeline.followup report
.\.venv\Scripts\python.exe -m pipeline.followup verify
```

The historical scientific report generators retain their original output conventions. After running them, restore the current integrated report with `python scripts/build_report.py`. The command reads the completed result tables; it does not generate receptor sequences, score candidates or change the freeze.

Online retrieval requires network access. EvoEF2 uses separate process work directories, short input filenames and resolved executable paths. Pinned source and executable hashes remain required. The [code manifest](code_manifest.tsv), [protocol](protocol.json) and recorded control-addition timestamps distinguish scientific choices from subsequent technical corrections.

## Published files and verification

Final result tables, necessary WT/mutation controls, original frozen candidates and report figures are tracked. Expanded viral state caches, fixed-E2 coordinate duplicates and raw EvoEF2 logs are regenerated locally and ignored. [The reproduction manifest](reproduction_manifest.tsv) retains their expected hashes. [The artifact manifest](artifact_manifest.tsv) describes the retained files.

The offline fixture and tests need no raw viral downloads or model weights. Full scientific verification also needs the ignored source/software and reconstructed structure caches. Use `python scripts/verify_publication.py --check-only` to audit the retained publication tree, or `python scripts/reverify_science.py` for the cached full scientific audit. See [publication checks](../publication/README.md) and the [recorded scientific verification](verification.json) for the scope and measured results.

The [original snapshot](original_snapshot.tsv) still resolves every original scientific record to its exact bytes, including archived documentation. Earlier reports are retained in `archive/`; they are historical records, while the root README is the current report. Initial diagnostic stages lacked continuous peak-resource telemetry, and no unmeasured peak is claimed.
