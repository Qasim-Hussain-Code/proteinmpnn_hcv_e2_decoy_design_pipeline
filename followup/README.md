# Supplementary tables and reproduction

The [study report](../README.md) presents the analysis and interpretation. This directory contains its experimental evidence, structural diagnostics and natural-sequence coverage data.

## Evidence and structural diagnostics

| Artifact | Contents |
|---|---|
| [Experimental ground truth](ground_truth.tsv) | Assay-specific labels, construct identities, exclusions and source provenance |
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

Use the pinned CPU environment, verified source snapshots and software described in the [software manifest](../results/software_manifest.tsv). Exact snapshot hashes are required to reproduce the study; a fresh live query may return different records. Candidate freeze hashes bind the study's sequence selection and inputs. An independent analysis requires its own configuration, source audit and freeze.

The design and evaluation workflow uses Git Bash:

```bash
bash scripts/00_configure.sh --threads 2 --ram 14 --disk 13 --seed 20261004 --yes
bash scripts/02_install.sh
bash run_all.sh --mode core --parallel
bash scripts/19_verify.sh
```

With the source and software caches available, the structural and coverage stages run from the project root:

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

Run `python scripts/build_report.py` after the analysis stages to render the study report and supplementary indexes from the completed result tables.

Online retrieval requires network access. EvoEF2 uses separate process work directories, short input filenames and resolved executable paths. The [code manifest](code_manifest.tsv) records the required source hashes, and the [protocol](protocol.json) specifies the structural comparisons and their scope.

## Published files and verification

Final result tables, necessary WT/mutation controls, frozen candidates and report figures are tracked. Viral state caches, fixed-E2 coordinate duplicates and raw EvoEF2 logs are regenerated locally and ignored. [The reproduction manifest](reproduction_manifest.tsv) retains their expected hashes. [The artifact manifest](artifact_manifest.tsv) describes the retained files.

The offline fixture and tests need no raw viral downloads or model weights. Full scientific verification also needs the ignored source/software and reconstructed structure caches. Use `python scripts/verify_publication.py --check-only` to audit the retained publication tree, or `python scripts/reverify_science.py` for the cached full scientific audit. See [publication checks](../publication/README.md) and the [recorded scientific verification](verification.json) for the scope and measured results.

The [source snapshot](original_snapshot.tsv) records file hashes for the design and evaluation inputs and outputs. Complete diagnostic peak-resource measurements are unavailable; the report's resource values describe only the metered workflow.
