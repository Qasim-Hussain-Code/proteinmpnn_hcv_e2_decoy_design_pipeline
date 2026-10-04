"""Extra native-WT repair sensitivity on the immutable held-out candidate set."""
import json
import shutil
import numpy as np
from .common import ROOT,MeasuredStage,read_tsv,write_json,write_tsv
from .scoring import evo,combine,energy,predicted_direction
from .selection import require_freeze

def main():
    require_freeze()
    with MeasuredStage('wt_preparation_sensitivity','One additional native WT RepairStructure pass; frozen held-out scores unchanged'):
        work=ROOT/'data/work/wt_preparation_sensitivity'; work.mkdir(parents=True,exist_ok=True)
        base=ROOT/'data/processed/humanized_repaired.pdb'
        evo('RepairStructure',base,work)
        native=work/'humanized_repaired_Repair.pdb'
        target=ROOT/'data/processed/controls/wild_type_second_repair.pdb'; shutil.copyfile(native,target)
        full_reference=energy(target)
        scores=read_tsv(ROOT/'results/heldout_scores.tsv')
        original={r['haplotype_id']:float(r['interaction_score']) for r in scores if r['candidate_id']=='wild_type'}
        alt={}; rows=[]
        for state in read_tsv(ROOT/'results/e2_state_manifest.tsv'):
            if state['split']!='held_out' or state['haplotype_id']=='reference': continue
            combined=work/'native_combined.pdb'; combine(target,ROOT/state['path'],combined)
            score=energy(combined); alt[state['haplotype_id']]=score
            rows.append(dict(haplotype_id=state['haplotype_id'],genotype=state['genotype'],primary_wt_score=original[state['haplotype_id']],
                             extra_repair_wt_score=score,wt_score_change=score-original[state['haplotype_id']],
                             interpretation='Same E2 coordinates; extra native receptor repair; no candidate reselection'))
        write_tsv(ROOT/'results/wt_preparation_sensitivity.tsv',rows)
        frozen=read_tsv(ROOT/'results/candidate_freeze.tsv'); candidates=[]
        for candidate in frozen:
            values=[float(r['interaction_score'])-alt[r['haplotype_id']] for r in scores if r['candidate_id']==candidate['candidate_id']]
            candidates.append(dict(candidate_id=candidate['candidate_id'],selection_rule=candidate['selection_rule'],worst_delta_wt=max(values),median_delta_wt=float(np.median(values))))
        write_tsv(ROOT/'results/wt_preparation_candidate_metrics.tsv',candidates)
        strategy={s:np.mean([r['worst_delta_wt'] for r in candidates if s in r['selection_rule']]) for s in ['single_state','escape_aware']}
        benchmark=read_tsv(ROOT/'results/cd81_mutation_benchmark.tsv')
        directional={r['mutation']:r for r in benchmark if 'entry' not in r['assay_context'] and r['experimental_direction'] in {'enhanced','reduced'}}
        correct=sum(predicted_direction(float(r['mutant_score'])-full_reference)==r['experimental_direction'] for r in directional.values())
        write_json(ROOT/'results/wt_preparation_summary.json',dict(primary_protocol='Native WT uses repaired reference; mutated sequences use BuildMutant and additional postmutation RepairStructure',
            sensitivity='One additional native WT repair pass; exact viral-state coordinates and frozen candidate outputs retained',
            native_reference_score_extra_repair=full_reference,alternative_directional_concordance=correct,directional_controls=len(directional),
            alternative_worst_effect_escape_minus_single=float(strategy['escape_aware']-strategy['single_state']),
            maximum_absolute_wt_state_score_change=max(abs(r['wt_score_change']) for r in rows),
            primary_freeze_and_benchmark_unchanged=True,interpretation='Exploratory preparation sensitivity; cannot rescue failed primary validation or establish experimental binding'))

if __name__=='__main__': main()
