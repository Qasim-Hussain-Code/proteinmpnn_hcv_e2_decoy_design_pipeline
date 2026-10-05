# Follow-up: scoring remains unreliable; broader coverage is prepared

The requested scoring, structural, and viral-coverage follow-up is complete. The original experiment is preserved byte for byte by [the snapshot](original_snapshot.tsv). No original candidate is re-ranked, and no new designed receptor sequence is generated. These are exploratory diagnostic comparisons with already consulted literature labels.

## Findings

The original model still recovers **5/11** distinct directional binding controls. A human-sequence model on the first bound tamarin backbone recovers **8/11**, but the second experimental bound pose recovers **3/11**. The first model misses all **four** additional binding-loss challenge mutations from [Drummer 2005](https://pubmed.ncbi.nlm.nih.gov/15670777/). Always predicting reduced binding recovers **10/11** original directional labels. Changing the model or the score threshold after seeing these labels cannot establish validation. The failed original experiment remains failed.

The original WT contains a **1.09 A** T163-to-E2 contact after repair and a **93.15** inter-chain van der Waals repulsion term. In the follow-up fixed-E2 T163A comparison, its score delta is **-6.16**, with a repulsive-term delta of **-9.35**. This is consistent with clash removal contributing to the apparent favorable score; it does not prove that the experimentally enhanced interaction has the same mechanism. D196 is **7.70 A** from E2 in the archived model, beyond EvoEF2's 6 A interaction cutoff. Folding, dimerization, receptor dynamics, glycans, membrane context, and assay expression effects are not represented by this interface score.

Sidechain preparation is an additional confound. Repacking both proteins versus recombining every mutant receptor with the exact WT E2 coordinates changes several model/control classifications. WT receives one repair in the initial model, whereas a mutant receives BuildMutant and another repair. [An additional WT repair control](wt_repair_sensitivity.tsv) measures this asymmetry explicitly; it is not a new validation gate. Model-specific shifts and all score terms are retained. Scores are model units, not measured affinities or calibrated free energies.

## Evidence correction

The original curation incorrectly assigned two soluble GST-LEL rows to single F186L and single E188K. [The primary Methods and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC111874/) describe the **F186L+E188K double mutant** for the soluble construct. The follow-up excludes those two misassigned rows and adds the double-mutant observation. Valid full-length cell-surface single-mutant observations are retained. The corrected original evidence has **28 included assay rows**; adding six source-external challenge observations gives **34**. Correcting this error does not change the original 11 distinct single-mutation directional labels or their 5/11 result. See [the correction audit](curation_audit.tsv) and [the corrected table](ground_truth.tsv). The abstract/Results discrepancy for D196E is retained as assay-context uncertainty. The original files remain an archival record, not the corrected evidence table.

The six new observations are explicitly abstract-level evidence. K124T and V146E retained binding, without establishing WT-equivalent affinity; they are not mislabeled neutral or counted in directional accuracy. The other four report binding loss in recombinant LEL. Full-length and soluble phenotypes are not pooled. Mutations affecting disulfides and dimerization challenge the score's scope.

## Structural comparison

Eight prespecified models use archived 3X0E, both 3X0E chains, both [1G8Q](https://www.rcsb.org/structure/1G8Q) chains, [5TCX](https://www.rcsb.org/structure/5TCX), and both experimental [7MWX](https://www.rcsb.org/structure/7MWX) binding pairs. All human entities and mapped native identities are verified. Missing residues are recorded and are never filled with invented coordinates. The core-fit sensitivity excludes the variable head by the ranges specified in [the locked protocol](protocol.json). Bound-template models use observed tamarin backbone coordinates with five human substitutions and are **not experimentally determined human complexes**. Their missing receptor termini are documented.

| Model | Clashes <2 A | WT repulsive term | Controls, repacked | Controls, fixed E2 | Additional loss controls, fixed E2 |
|---|---:|---:|---:|---:|---:|
| archive | 2 | 93.15 | 5/11 | 5/11 | 1/4 |
| x0e_b | 41 | 379.70 | 1/11 | 1/11 | 0/4 |
| x0e_core | 36 | 441.95 | 1/11 | 4/11 | 1/4 |
| g8q_a | 3 | 112.28 | 1/11 | 5/11 | 1/4 |
| g8q_b | 39 | 368.22 | 1/11 | 7/11 | 4/4 |
| tcx_a | 27 | 378.24 | 1/11 | 6/11 | 0/4 |
| bound_ae | 0 | 10.27 | 8/11 | 8/11 | 0/4 |
| bound_bh | 0 | 13.89 | 3/11 | 3/11 | 0/4 |

All native WT models retain the two native disulfide geometries. Alternative human rigid fits can worsen clashes markedly; lower CA RMSD alone is not a reliable quality criterion. Zero clashes in a humanized bound-template model do not establish its accuracy or its affinity. No best model is selected using these benchmark outcomes.

![Structural scoring sensitivity](figures/benchmark.png)

## Viral coverage

The bounded retrieval expanded from 309 to **878 unique accessions**, with **532 eligible**, **346 excluded**, and **444 distinct mapped E2 sequences**. The interface definition expands from 20 to **37 H77 positions**, combining the archived human geometry, both native bound interfaces, and the five already curated functional positions. Queries request at most 120 default-order records per stratum and explicitly inspect patent references; the ignored original patent query filter is not reused. Query strings, UIDs, dates, source bytes, and hashes are retained. This is a convenience sample, not a worldwide prevalence estimate.

The original 24-state panel covered **39.44%** of its development accessions; reaching 90% on that same archived cohort requires **118 states** under the same selection rule. The expanded 37-position development panel requires **284 states**, covering **90.11%** of 526 development accessions. Its genotype 1 coverage is only 72.13%. A documented secondary extension to **306 states** raises overall coverage to **94.30%** and the minimum coverage within each included recorded genotype to **90.16%**. All of these coordinate models are built and hashed. The panel extension uses sequence counts only, with no candidate or assay score selection.

Coverage is conditional on structural eligibility. **Genotype 8 remains outside this coordinate panel**: its four metadata-confirmed records contain interface insertions that a fixed-backbone sidechain-substitution model cannot represent. Those accessions and reasons are retained; no claim of genotype 8 coverage is made. There are **22 unknown-genotype eligible accessions**. Three expanded interface positions, **415-417**, lack coordinates in the archived E2 chain; every state records this missing structural coverage. Thus sequence haplotype coverage at 37 sites and coordinate coverage at the 34 observed sites are distinct. The prepared states are partial structural proxies, not complete models of all natural E2 differences.

The excluded-record counts are: {'excessive_ambiguous_amino_acids': 214, 'expanded_interface_missing': 3, 'no_annotated_CDS_translation': 60, 'patent_record': 56, 'severe_alignment_failure': 4, 'unmodelable_interface_insertion': 9}. Strict mapping and indel exclusions can bias the retained cohort, so broad protective efficacy cannot be inferred from these percentages.

A future reserve contains **6 newly retrieved accessions** in novel connected components, separated from all previously retrieved accessions and all old interface haplotypes within Hamming distance one. These reserve sequences are not modeled, scored, or used to select the panel. The old held-out data are already exposed and remain development data here. A six-accession reserve is too small for a convincing broad-coverage claim; an independent future cohort remains necessary.

![Accession coverage curve](figures/coverage.png)

![Coverage by recorded genotype](figures/genotype_coverage.png)

## Reproduction and artifacts

Use the existing locked CPU environment and pinned EvoEF2 executable. No new dependency stack or GPU predictor is installed. From the project root, run:

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
.\.venv\Scripts\python.exe -m pytest -q
```

Online retrieval needs network access. Original raw inputs and the pinned software are prerequisites, as in the original experiment; source caches are outside Git. Retrieval hashes enforce the downloaded snapshot. EvoEF2 receives short local filenames to avoid its upstream path buffer limitation. Separate process work directories prevent simultaneous state and receptor jobs from overwriting inputs. Repair/mutation caches are retained. The initial code hash, final code manifest, and control-addition timestamp distinguish prespecified choices from subsequent technical fixes and coverage extension. This repository adds a supplement; it does not rewrite the original manuscript.

The primary tables are [benchmark results](benchmark.tsv), [energy terms](energy_terms.tsv), [model audits](model_audit.tsv), [sequence QC](sequence_manifest.tsv), [coverage summary](coverage_summary.json), [balanced panel](balanced_development_panel.tsv), and [state hashes](state_manifest.tsv). Five figures are supplied as PNG and PDF. Verification checks archive immutability, model sequence identities, score logs, state identities, coordinate coverage limitations, leakage boundaries, and per-genotype target coverage. See [verification](verification.json) for the measured result. Resource measurements from the original experiment remain untouched; follow-up stage timings printed during execution and the three-state pilot are separate. Earlier follow-up stages did not continuously meter RSS, and no unmeasured peak is claimed.

## Next work

The scoring problem requires a method that represents receptor conformational change and fold integrity, with matched WT/mutant preparation and glycan/complex context, followed by validation on independently curated controls. Receptor-only fold/stability diagnostics and an ensemble constrained by experimental structures are the next computational steps. A change in this score or the score threshold is insufficient. Collecting an independent, geographically broader cohort and modeling interface indels are necessary before genotype 8 can enter a structural breadth claim. Future receptor ranking should use a separately frozen protocol and a fresh independent test set. Experimental binding and soluble-protein behavior would then determine whether any therapeutic-decoy hypothesis is supported.


## Published files and local regeneration

The GitHub tree retains final result tables, figures, experimental WT models, frozen original candidates, and the mutation controls needed for reproducibility. The 306 expanded viral states, fixed-E2 duplicate models, and raw EvoEF2 execution logs are rebuildable local outputs and are excluded from Git. Their expected hashes are retained in [the reproduction manifest](reproduction_manifest.tsv). The complete scientific audit was rerun against these local outputs before preparing the publication tree.

The existing `prepare`, `models`, `sequences`, `coverage`, `benchmark`, `states`, and `controls` stages regenerate those ignored outputs. Full `followup verify` requires the regenerated outputs and the original source/software caches. A clone's offline fixture and Python tests require no viral downloads or model weights; use the original locked CPU environment, or install the documented dependencies. Publication verification separately checks the retained files, history, tracking policy, and data consistency.

Immutable copies of the original README and host-dependent startup/resource records are preserved under `archive/original/`. The original scientific freeze files and inputs retain their exact bytes. Subsequent smoke-test execution may regenerate the live host-dependent records without changing the archived evidence.
