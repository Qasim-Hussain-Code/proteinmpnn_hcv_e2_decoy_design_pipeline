# Verification and reproducibility

The [study report](../README.md) presents experimental-control benchmarking, held-out candidate evaluation, structural sensitivity and natural-sequence coverage. [Report provenance](unified_report_provenance.json) records its source hashes and displayed metrics. Run `python scripts/build_report.py` to render the report from the completed result tables.

## Verification scope

The offline analytical fixture and all 43 pytest tests pass in a clean clone using the documented CPU environment. Bash syntax and ShellCheck pass for all six shell scripts. The publication audit checks retained files, source and artifact hashes, scientific table consistency, commit format, tracking policy and file-size limits.

Scientific verification includes 87 design/evaluation checks and 2,075 applicable structural and coverage checks against the published artifacts and locally retained caches. The complete local cache audit covers 4,550 checks. These checks establish computational consistency and provenance; the experimental scoring benchmark fails and does not establish binding improvement.

[Design and evaluation traceability](../results/readme_traceability.tsv), [structural and coverage traceability](readme_traceability.tsv), [report provenance](unified_report_provenance.json) and [data-license review](data_licenses.tsv) connect the reported values and artifacts to their sources. Verification records identify their tested commits and scope.

## Running checks

```powershell
python -m pytest -q
python scripts/verify_publication.py --check-only
python scripts/build_report.py
```

The offline fixture and publication audit require no raw viral downloads or model weights. The full scientific audit additionally requires the pinned source/software and reconstructed structure caches:

```powershell
python scripts/reverify_science.py
```

Use the locked CPU environment and the [stage commands](../followup/README.md#scientific-reproduction). Exact downloaded snapshot hashes are required for reproducing the recorded scientific run. Contemporary queries can return different records.

## Files and resources

The repository tracks final scientific tables, necessary controls, frozen candidates, report assets and compact provenance records. Bulk downloads, virtual environments, weights, intermediate coordinate populations and state caches are excluded. All tracked files are below the 50 MB limit. Regenerated structures are checked against their expected hashes.

Git history uses the existing author identity, starts with `.gitignore`, and uses two- or three-word underscore commit messages. Third-party software and weights are downloaded separately under their recorded usage terms. Original implementation code is MIT licensed.

The design/evaluation workflow has measured peak-resource records. Complete peak measurements are unavailable for the structural diagnostics, so measured resource values are scoped to the metered workflow. EvoEF2 uses short input filenames and resolved executable paths; its source and executable hashes are pinned.
