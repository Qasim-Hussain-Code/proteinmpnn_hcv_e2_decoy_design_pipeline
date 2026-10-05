# Benchmarking CD81 sequence design across HCV E2 diversity

This study combines constrained sequence design, experimental-control benchmarking, structural diagnostics and natural-sequence coverage in one analysis. It asks whether a human CD81 scaffold can support a credible computational receptor-decoy comparison across naturally observed HCV E2 variants.

**Scoring remains unreliable.** EvoEF2 recovered ${original_correct}/${control_total} directional controls in the frozen experiment, compared with ${baseline_correct}/${control_total} from always predicting reduced binding. Across ${model_count} structural models, the strongest bound-template comparison recovered ${best_bound_correct}/${control_total}, but missed all ${external_missed} additional binding-loss controls. Broader sampling improved coverage of the eligible development accessions to ${balanced_coverage}; it did not validate the scorer or demonstrate experimental binding for any generated sequence.

## Study design

HCV E2 interacts with the large extracellular loop of CD81. Here, ProteinMPNN proposes amino-acid sequences on an experimental human CD81 backbone. This is constrained de novo sequence design; it does not generate a new fold. The designed molecule is the receptor scaffold. E2 remains fixed during sequence generation.

The human complex is a modeled hybrid. [7MWX](https://www.rcsb.org/structure/7MWX) supplies the experimentally bound E2 orientation and **tamarin**, rather than human, CD81. The human scaffold from [3X0E](https://www.rcsb.org/structure/3X0E) is aligned onto that orientation, with a matched C-alpha RMSD of ${rmsd} Å. Missing coordinates are left missing, glycans and other excluded heteroatoms are audited, and side chains are repaired without docking or large backbone minimization. The modeled protein-only interface does not represent complete-virion accessibility.

Standard and soluble ProteinMPNN use matched seeds, temperatures and generation budgets. Arm A preserves the geometric interface and recognition constraints; Arm B permits the remaining reliable contact positions to vary. Both retain native disulfide cysteines. Single-state selection favors the reference score; escape-aware selection favors the worst candidate-minus-WT score across discovery states, then its median. Sequence probability, physical scores, charge and exposure are recorded separately.

WT means the unmodified human CD81 sequence. A paired score delta subtracts the WT score under the same viral-state geometry; positive values favor WT within this model. A haplotype is the combination of amino acids at the selected interface positions, and a structural state is its partial coordinate model. These distinctions keep sequence coverage, model performance and experimental evidence separate.

The analysis comprises experimental-control benchmarking, held-out candidate evaluation, structural sensitivity and natural-sequence coverage. Candidate sequences and selection criteria are frozen before held-out evaluation. Structural comparisons use consulted experimental labels and are exploratory; they do not constitute independent blinded validation. Coverage panels are constructed from sequence frequencies without candidate scores. The [design protocol](config/design.yml), [structural protocol](followup/protocol.json), [candidate freeze](results/candidate_freeze_manifest.json) and [methods](docs/methods_notes.md) define these analysis roles.

## Experimental evidence and scoring

The [experimental evidence table](followup/ground_truth.tsv) contains **${included_rows} included assay observations**: ${benchmark_assay_rows} in the benchmark evidence group and six in the source-external challenge group. It distinguishes soluble-protein binding, cell-surface binding and viral entry; heterogeneous assay values are not pooled or converted into invented affinities.

[Higginbottom's Methods and Results](https://pmc.ncbi.nlm.nih.gov/articles/PMC111874/) describe the soluble GST-LEL **F186L+E188K double-mutant construct**. This observation is treated as a double-mutant assay. Cell-surface single-mutant observations are assigned to their reported constructs and assay contexts. The directional benchmark comprises ${control_total} distinct single-mutation labels, of which the reference model recovers ${original_correct}.

The additional observations come from [Drummer 2005](https://pubmed.ncbi.nlm.nih.gov/15670777/). Four report binding loss; two report retained binding without establishing WT-equivalent affinity. The latter are not labeled neutral or counted in directional accuracy. These consulted labels form an exploratory challenge, not an independent blinded validation set.

Eight structural models assess alternative human receptor conformations, core alignment and both observed 7MWX receptor-binding pairs. Models based on a bound tamarin backbone carry the human sequence, but are not experimentally determined human complexes. Both preparation modes use the declared 0.5 score-unit tolerance: repacking both proteins, or keeping every mutant's E2 coordinates identical to its model's WT E2.

| Model | Clashes below 2 Å | WT repulsive term | Benchmark controls, repacked | Benchmark controls, fixed E2 | Challenge loss controls, fixed E2 |
|---|---:|---:|---:|---:|---:|
${model_rows}

The first bound-template pose recovers ${best_bound_correct}/${control_total} controls; the second recovers ${second_bound_correct}/${control_total}. The first misses all four added binding-loss controls. Model-dependent results and missed challenge controls prevent treating the best-looking outcome as validation. No model or threshold was selected to improve the reported benchmark.

![Directional-control recovery across structural models](followup/figures/benchmark.png)

*Structural sensitivity of the scoring benchmark. Results depend on receptor geometry and preparation; the diagnostic comparisons do not replace the failed frozen validation gate. [Underlying summaries](followup/benchmark_summary.tsv).*

The frozen WT has a ${min_contact163} Å T163-to-E2 contact and a ${wt_repulsion} inter-chain repulsive term. This suggests that clash relief can contribute to a favorable mutation score. It does not establish the experimental mechanism. D196 lies ${min_contact196} Å from E2, outside the scorer's 6 Å interaction cutoff. Receptor folding, oligomerization, dynamics, expression and membrane context can affect experimental outcomes without appearing in this interface score. [Contacts](followup/control_contacts.tsv), [energy terms](followup/energy_terms.tsv) and [matched WT repair controls](followup/wt_repair_sensitivity.tsv) make these limitations inspectable.

## Frozen sequence-design comparison

ProteinMPNN generated ${generated} sequences; ${unique_passed} unique sequences passed the sequence constraints. The frozen shortlist contains ${frozen} distinct computational candidates: ${strategy_candidates} selections per strategy, with ${shared_candidates} shared selections. Every untested sequence remains a **computational candidate**.

Evaluation used ${heldout_states} structural states from ${heldout_accessions} held-out accessions. The escape-aware-minus-single-state difference in mean candidate worst paired score was **+${effect} EvoEF2 units**, with paired state-bootstrap interval **[${ci_low}, ${ci_high}]**. Lower scores are favored within this model, so escape-aware selection did not improve this held-out measure. Correlated haplotype components and the small test set limit the interval's interpretation. The failed experimental benchmark prevents interpreting either strategy's scores as experimental binding improvement.

![Frozen held-out strategy comparison](figures/figure_6.png)

*The candidate sets and held-out states were fixed before evaluation. This is a comparison of model outputs under a failed experimental scoring gate. [Strategy statistics](results/heldout_strategy_comparison.tsv), [bootstrap summary](results/heldout_summary.json), [candidate table](results/candidate_summary.tsv).*

The soluble model reduced the hydrophobic SASA fraction relative to the standard model by ${exposure_A} in Arm A, with paired seed-bootstrap interval ${exposure_ci_A}, and ${exposure_B} in Arm B, with interval ${exposure_ci_B}. Other descriptors were not uniformly favorable: Arm A's hydrophobic patch proxy increased. These geometric descriptors do not establish expression, folding or solubility. [All model comparisons](results/model_comparison_effects.tsv) retain the matched-seed estimates and uncertainty.

## Natural-sequence coverage

Coverage counts the fraction of eligible sampled accessions represented by exact interface haplotypes. It is conditional on the retrieved cohort, exclusions and interface definition; it is not worldwide prevalence or demonstrated protective efficacy.

| Quantity | Design and evaluation cohort | Diversity coverage cohort |
|---|---:|---:|
| Retrieved accessions | ${original_retrieved} | ${retrieved} |
| Eligible accessions | ${original_included} | ${eligible} |
| Development accessions | ${original_development} | ${development} |
| Interface positions | ${original_positions} | ${expanded_positions} |
| Development structural states | ${original_states} | ${balanced_states} |
| Exact development haplotype coverage | ${original_coverage} | ${balanced_coverage} |

The diversity retrieval contains ${unique_e2} distinct mapped E2 sequences and ${excluded} excluded accessions. A frequency-based global panel of ${global_states} states covers ${global_coverage}, with lower coverage in genotype 1. Exploratory genotype balancing extends the panel to ${balanced_states} states, with overall coverage of ${balanced_coverage} and a minimum of ${balanced_minimum} within each included recorded genotype. Panel construction uses sequence counts, without candidate or assay scores. The two cohorts differ in accession membership and interface definition, so their percentages do not estimate the effect of panel size alone.

![Eligible accession coverage as states are added](followup/figures/coverage.png)

*Coverage within the design and diversity development cohorts. Candidate performance is evaluated on the design cohort's held-out states; the diversity panel characterizes sequence coverage rather than additional candidate performance. [Coverage definitions and counts](followup/coverage_summary.json), [final panel](followup/balanced_development_panel.tsv).*

![Coverage within each recorded genotype](followup/figures/genotype_coverage.png)

*Genotype coverage is conditional on structural eligibility. Genotype 8 is excluded from the coordinate panel; eligible unknown-genotype accessions remain explicitly unclassified. [Per-genotype counts](followup/balanced_genotype_coverage.tsv).*

The exclusions matter. All four metadata-confirmed genotype 8 records contain interface insertions that this fixed-backbone side-chain model cannot represent. No genotype 8 structural coverage is claimed. There are ${unknown_genotype} eligible accessions with unknown genotype. Interface positions 415-417 lack coordinates in the archived E2 chain: haplotypes use ${expanded_positions} positions, while states represent only ${modeled_positions} observed positions. Broader sequence coverage therefore remains incomplete structural coverage.

A reserve of ${future_reserve} accessions occupies novel connected components separated from the analyzed data. It is excluded from modeling, scoring and panel selection. The design cohort's held-out set supports its fixed candidate comparison; it is not independent validation data for the diversity analysis, which includes those accessions in its development cohort. The reserve is too small to establish broad candidate performance, and an independently collected future cohort is needed.

## Interpretation and next steps

The study establishes a reproducible sequence-design comparison and identifies why its physical scoring framework is insufficient. It also prepares a broader, explicitly limited natural-diversity panel. It does **not** demonstrate that any generated sequence binds HCV E2, neutralizes HCV, folds correctly or functions as a therapeutic.

The greatest methodological uncertainty is receptor conformation. Alternative rigid fits can worsen clashes, and a humanized bound-template model can look geometrically cleaner while still missing experimental controls. Side-chain preparation, absent glycans, fixed backbones and assay-specific biology add uncertainty. The least impressive result remains central: the frozen scorer performs worse than the simple majority-direction baseline, and escape-aware selection does not improve the reported held-out worst-score measure.

The next priorities are a conformational ensemble constrained by experimental structures, matched WT/mutant preparation, receptor fold-integrity diagnostics and independently curated controls assessed under a separately frozen protocol. Broader geographic sampling and a method that can represent interface insertions are needed before expanding genotype coverage claims. Experimental binding and soluble-protein measurements would determine whether any computational candidate supports the receptor-decoy hypothesis. Human CD81 interactions, specificity, immune effects and safety remain untested.

## Reproduction and data

The primary implementation is CPU-based. The locked environment targets Python 3.14 on Windows with Git Bash; other platforms require a separate software audit. No large GPU predictor, molecular dynamics stack or de novo backbone generator is required. [Supplementary tables and stage commands](followup/README.md), [software provenance](results/software_manifest.tsv) and [score definitions](docs/score_dictionary.md) describe the inputs and outputs.

For the offline analytical fixture and tests:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install numpy biopython scipy matplotlib psutil pyyaml pytest
.\.venv\Scripts\python.exe -m pipeline.cli configure --yes
.\.venv\Scripts\python.exe -m pipeline.cli smoke
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe scripts/verify_publication.py --check-only
```

The minimal fixture dependencies above do not constitute the locked scientific environment. The full scientific stages require the pinned software, downloaded source snapshots and model weights, all excluded from Git. [Scientific reproduction](followup/README.md#scientific-reproduction) gives the design, structural and coverage stage commands. Run `python scripts/build_report.py` to render the study report from its result tables.

The metered design and evaluation workflow reached ${peak_rss_gb} GB peak sampled process-tree RSS and ${peak_project_gb} GB peak measured project footprint. Complete peak-resource measurements are unavailable for the structural diagnostics, so these values do not describe the whole analysis. The resource guard protects a 1 GB non-project disk reserve and a 13 GB project ceiling. [Publication verification](publication/README.md) defines the measured checks and their scope.

Raw downloads, environments, weights, intermediate coordinate populations and state caches are ignored. Necessary controls, frozen candidates, final tables and expected regeneration hashes are retained. The [figure index](figures/README.md), [design and evaluation traceability](results/readme_traceability.tsv), [structural and coverage traceability](publication/readme_traceability.tsv) and [report provenance](publication/unified_report_provenance.json) connect claims to source data.

## Sources and license

Primary sources include [Kumar et al.](https://doi.org/10.1038/s41586-021-03913-5), [Yang et al.](https://doi.org/10.1096/fj.15-272880), [Higginbottom et al.](https://doi.org/10.1128/JVI.74.8.3642-3649.2000), [Drummer et al. 2002](https://doi.org/10.1128/JVI.76.21.11143-11147.2002), [Drummer et al. 2005](https://pubmed.ncbi.nlm.nih.gov/15670777/), [Bertaux and Dragic](https://doi.org/10.1128/JVI.80.10.4940-4948.2006) and [Flint et al.](https://doi.org/10.1128/JVI.00104-06). The evidence table retains assay-specific provenance.

Original code uses the [MIT license](LICENSE). Third-party software and weights are retrieved independently. ProteinMPNN's code license and bundled-weight terms are audited separately; a separate weight grant is not assumed. EvoEF2's MIT license file conflicts with its academic-use README language, so its source and binary are not redistributed. [Data and software license reviews](publication/data_licenses.tsv) cover coordinates, sequence records and article copyright. Source article prose and figures are not redistributed.
