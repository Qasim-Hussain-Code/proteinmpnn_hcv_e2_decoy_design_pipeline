# Publication checks

The original scientific verification and all 4,550 follow-up integrity checks were repeated before preparing this publication branch. Their audit tables and timings are retained here. The score benchmark remains failed; publication does not establish binding improvement.

The last research checkpoint's five-word message and regenerable cache tracking did not meet the requested Git rules. This independent publication branch starts at the compliant original history, retains the completed scientific tables and necessary controls, and uses two- or three-word underscore commit messages. The original local research branch is preserved and is not uploaded.

The full scientific workflow uses separately downloaded source snapshots, a locked CPU environment and pinned external tools. Large downloads, model weights, virtual environments and intermediate coordinate populations are excluded. The published supplementary tables retain expected hashes for locally regenerated structures. Historical source bytes, including startup records that a new clone regenerates, are preserved in the archive.

Run `python scripts/verify_publication.py` to inspect the published tree, history, hashes and scientific table consistency. Run the documented offline smoke fixture and pytest suite in a clean clone. Full scientific verification additionally requires the ignored local source/software and structure caches; regenerate them with the documented stages. Clean-clone testing is performed on Windows with Git Bash and the existing locked CPU environment.

See [data-license review](data_licenses.tsv), [follow-up numeric traceability](readme_traceability.tsv), [original numeric traceability](../results/readme_traceability.tsv), and the generated publication check report. Initial follow-up stages lacked continuous peak-resource telemetry; those historical peaks cannot be reconstructed and are not claimed. The repeat audits retain observed timing and provenance, not invented historical peaks.

The publication branch passed a fresh clone check: the offline smoke fixture, all 43 pytest tests, Bash syntax and ShellCheck for all six shell scripts, and 820 publication audit checks. The clone contained no raw downloads, third-party vendor directory, virtual environment or expanded structural-state cache. [The recorded check](clean_clone_verification.json) identifies its exact source commit and code hashes. The retained publication tree also passed all 87 original scientific checks and 2,075 applicable follow-up checks using the separately retained local caches.

Reverification exposed a long executable-alias path that EvoEF2 could not handle. The follow-up launcher now resolves that alias before execution, and a regression test checks the behavior. This runtime correction does not change the frozen experiment or scoring interpretation. Use `python scripts/verify_publication.py --check-only` for a read-only repeat audit.
