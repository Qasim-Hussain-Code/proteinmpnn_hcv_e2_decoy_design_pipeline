# Methods and interpretation

## How to read this experiment

A coordinate file is an experimental observation with missing residues and uncertain atoms. Sequence design asks which amino acids the inverse-folding model considers compatible with that coordinate scaffold. Its sequence probability is a different quantity from a physical interaction score. Neither establishes experimental binding for a computational candidate.

The controls are therefore tested first. A directional label from a soluble-E2 assay is compared with a mutant-minus-WT EvoEF2 output. Pseudoparticle entry observations remain separate because membrane receptor trafficking and later entry events can change their interpretation. Qualitative labels are not turned into invented measurements. A failure of the interface scorer restricts all later claims even when one selection method produces lower model scores.

## Coordinate preparation

The mmCIF entity source organism identifies the receptor species. The expression host does not. Entity and label-chain mappings are resolved from the downloaded file, then a receptor/E2 pair is selected by the number of geometric contacts with lexical tie breaking. The biological assembly operation metadata are retained; the selected pair is an observed contacting pair in the asymmetric unit. Antibody chains are removed from the protein-only comparison. Their presence in the experimental crystallization complex is a limitation.

Human scaffold chain selection is lexical among verified human CD81 entities. The global affine-gap sequence alignment matches observed residues, and all matched C-alpha coordinates define the least-squares rigid transform. No residues are removed to improve RMSD. Human numbering is independently checked against UniProt P60033. Observed missing coordinates remain missing. The humanized protein is a hybrid model, not a newly determined experimental complex. Its backbone coordinates are not minimized or docked.

The main interface is defined on the unrepaired humanized model, independently of conservation. Repair may change sidechain contacts without retroactively changing the main interface definition. Sensitivity cutoffs remain annotations, not a way to choose an attractive outcome. All removed nonprotein residue identities, atom counts, label chains and proximity to selected E2 are recorded. Geometric glycan proximity is not a complete glycan accessibility calculation.

## Viral sequences and numbering

GenBank queries request a bounded genotype-enriched convenience sample. Record annotations, rather than the query term alone, determine genotype labels. Failed extraction and unreliable mapping retain explicit exclusion reasons. Annotated CDS translations are preferred to guessing a reading frame. Polyprotein neighborhood alignment maps E2 onto H77; observed experimental E2 residues are mapped independently. A sidechain mutation always checks the actual experimental-template residue first, which may differ from H77. H77 numbering and template identity are separate fields.

The default clustering joins identical and one-substitution interface haplotypes into connected components. An entire component remains in one split. This prevents exact and near-haplotype leakage but may create a large component, so the accession split may be far from the nominal ratio. Stratification uses modal genotype/subtype only when enough components exist; smaller strata enter a deterministic pooled fallback. This fallback and the achieved counts are retained rather than silently rebalancing the test set.

Discovery panels use genotype round robin, choosing the most frequent unrepresented haplotype within each genotype. The target coverage is an arbitrary engineering threshold subject to the state cap. Sensitivity applies the same deterministic algorithm to alternate targets. Test haplotypes are not built or scored before the candidate freeze passes its hash validation.

## Sidechains and scoring

The verified upstream EvoEF2 source supports no run-count CLI option despite its README examples. RepairStructure uses its source default of one sequential optimization pass; backbone-dependent BuildMutant uses ten sequential rotamer passes. These paths do not call the stochastic sequence-design routine. Library hashes and the local static compilation command are recorded. The upstream Windows executable required an unavailable runtime; no runtime stack was installed.

Each designed receptor is mutated and repaired in the reference hybrid. Each viral state is independently mutated and repaired from that same hybrid. Candidate receptor coordinates are then combined with the state's E2 coordinates. There is no state-specific receptor repacking. This limits conformation adaptation but gives an interpretable identical procedure for comparisons. ComputeBinding uses chain split E,R. Monomer total energy and monomer SASA are separate descriptors. Score output precision is limited by upstream printed totals.

Clash filtering allows a predeclared arbitrary excess over the WT atom-pair clash count. Composition metrics remain distributions rather than unsupported universal good-protein thresholds. Pareto analysis is annotated and does not secretly alter the lexicographic selection rule. The two selection methods use one shared unique candidate pool and equal shortlist size. Overlapping selections remain shared; they are not counted as independent candidate experiments.

## Seeds, freezing, and uncertainty

ProteinMPNN's nonzero seed initializes Python, NumPy and PyTorch in the verified runner. The effective seed is checked in the original FASTA header. Model, arm, temperature and independent generation seed are retained for every sample, including duplicates. The pilot is excluded from the candidate pool. A throughput rule fixes the full budget before structural candidate scoring.

The freeze binds configuration, software manifest, discovery-score hashes and selected sequences. Frozen candidates cannot be edited once held-out scores are available. A file-open audit hook rejects test-file reads during discovery stages before a valid freeze exists. Split construction itself must inspect all raw sequences to form nonleaking components; that necessary preprocessing is distinct from candidate selection and remains explicit.

The primary strategy effect compares mean candidate worst paired deltas. A paired state bootstrap resamples structural states with replacement, preserving the same draw for both selection sets. The percentile interval describes this sampled model panel; related haplotypes and mixed genotype labels limit independence. Genotype-specific summaries are descriptive. No small P value or model score is treated as biological proof.
