# Score dictionary

All untested sequences are computational candidates. The experimental mutation gate failed in this run.

## ProteinMPNN_score

Mean -ln P(sequence residue | backbone, decoding order, preceding sequence) over masked redesigned positions; upstream _scores and mask*chain_M*chain_M_pos. unitless, natural-log NLL.

Model sequence compatibility. Does not demonstrate binding, affinity, folding or solubility.

## interaction_score

ComputeBinding split=E,R: weighted complex energy minus separated-chain contributions with coordinates unchanged. EvoEF2 score units; physical unit not asserted.

Compare this protein-only fixed-conformation model output. Does not estimate KD or rigorous binding free energy; failed mutant gate.

## stability_score

ComputeStability on extracted receptor monomer; weighted total terms including residue reference terms. EvoEF2 score units.

Model total-energy descriptor. Not measured melting point, expression, folding probability or stability.

## sasa_angstrom2

Shrake-Rupley atom sphere sampling, probe 1.4 A, 100 sphere points; monomer residue sum. angstrom squared.

Solvent-exposure geometry. Not measured solubility.

## hydrophobic_sasa_fraction

Sum monomer residue SASA for AVILMFWY / total receptor SASA. unitless.

Exposure proxy. Not measured solubility or aggregation.

## hydrophobic_patch_proxy_angstrom2

Maximum sum hydrophobic-residue SASA within 8 A CA sphere around a receptor residue. angstrom squared.

Patch-size proxy. Not experimentally measured hydrophobic patch or aggregation.

## net_charge_ph7

ProteinAnalysis charge_at_pH(7), residue and terminal Henderson-Hasselbalch model. elementary-charge proxy.

Composition descriptor. Not measured electrostatics of the folded molecule.

## pI

ProteinAnalysis.isoelectric_point, bisection of charge model. pH.

Composition descriptor. Not measured pI.

## hydrophobic_fraction

Count AVILMFWY / sequence length. unitless.

Composition descriptor. Not measured solubility.

## sequence_identity_percent

100 * exact residue matches to human scaffold / scaffold length. percent.

Sequence novelty relative to this scaffold. Not candidate quality.

## mutation_count

Number of residues differing from WT observed human scaffold. count.

Sequence change count. Not improved function.

## longest_unchanged_stretch

Longest contiguous sequence-index run identical to WT. residue count.

Sequence novelty descriptor. Not fold similarity.

## shannon_entropy

-sum p(aa)*log2 p(aa); excludes gaps and ambiguous aa; accession-weighted discovery only. bits.

Observed sequence diversity. Not population prevalence or escape phenotype.

## clash_count

Number of interchain heavy-atom pairs at distance <2.0 A. atom-pair count.

Geometric incompatibility descriptor. Not energetic penalty or experimentally observed clash.

## contact_count

Number of interchain heavy-atom pairs at distance <=5.0 A. atom-pair count.

Geometric interface descriptor. Not independent contacts or biological affinity.

## ca_rmsd_angstrom

Square root mean squared distance after least-squares rigid CA superposition. angstrom.

Experimental conformation mismatch. Not model confidence or fold validation.

## median_score

Median EvoEF2 interaction score across equally weighted structural states. EvoEF2 score units.

Model central tendency. Not binding across genotypes.

## worst_score

Maximum EvoEF2 interaction score across equally weighted structural states. EvoEF2 score units.

Model worst-state behavior. Not worst biological phenotype.

## best_score

Minimum EvoEF2 interaction score across states. EvoEF2 score units.

Model best-state behavior. Not biological performance.

## paired_delta_wt

Candidate interaction score minus WT score on same viral state. EvoEF2 score units.

Paired model comparison. Not experimental mutation thermodynamics.

## worst_delta_wt

Maximum paired candidate-minus-WT delta across states. EvoEF2 score units.

Computational prioritization. Not experimentally established escape resistance.

## median_delta_wt

Median paired candidate-minus-WT delta across states. EvoEF2 score units.

Paired model summary. Not affinity.

## iqr

75th minus 25th percentile interaction scores across states. EvoEF2 score units.

Model variability. Not biological uncertainty.

## std

Population standard deviation of interaction scores across states. EvoEF2 score units.

Model variability. Not independent experiment error.

## fraction_better_wt

Fraction of equally weighted states with paired delta <0. unitless.

Model comparison. Not binding frequency.

## fraction_within_wt_tolerance

Fraction states with paired delta <=0.5 arbitrary score tolerance. unitless.

Model comparison. Not biological noninferiority.

## rank_spearman

Spearman correlation of discovery and held-out worst paired-delta ranks across frozen candidates. unitless.

Rank stability in this sample. Not quantitative experimental validation.

## worst_delta_effect_escape_minus_single

Mean candidate maximum paired WT delta for escape set minus mean maximum for single set. EvoEF2 score units.

Between-selection-strategy model effect. Not biological improvement.
