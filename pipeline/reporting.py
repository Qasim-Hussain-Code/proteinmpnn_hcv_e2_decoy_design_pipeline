"""Figures and prose generated only from completed, auditable run artifacts."""
from __future__ import annotations
from collections import defaultdict
import json
from pathlib import Path
import re
import numpy as np
from .common import ROOT, config, now, read_tsv, sha256, write_json, write_tsv

def candidate_reporting():
    """Readable joins and seed-level comparisons, without changing the freeze."""
    frozen=read_tsv(ROOT/'results/candidate_freeze.tsv')
    held={r['candidate_id']:r for r in read_tsv(ROOT/'results/heldout_metrics.tsv')}
    interface={int(r['position']) for r in read_tsv(ROOT/'config/fixed_positions.tsv') if float(r['structural_distance_to_e2'])<=config()['interface_cutoff']}
    flat=[]; generalization=[]
    generated_annotations=[]
    for candidate in read_tsv(ROOT/'results/generated_sequences.tsv.gz'):
        mutations=candidate['mutation_list'].split(',') if candidate['mutation_list'] else []
        ninterface=sum(int(m[1:-1]) in interface for m in mutations)
        generated_annotations.append(dict(candidate_id=candidate['candidate_id'],interface_mutation_count=ninterface,
                                          scaffold_mutation_count=len(mutations)-ninterface,definition='Unrepaired hybrid heavy-atom 5 A receptor interface'))
    write_tsv(ROOT/'results/generated_sequence_annotations.tsv.gz',generated_annotations)
    for candidate in frozen:
        discovery=json.loads(candidate['discovery_metrics']); bio=json.loads(candidate['biophysics']); h=held[candidate['candidate_id']]
        columns={k:v for k,v in candidate.items() if k not in {'discovery_metrics','biophysics'}}
        mutations=candidate['mutation_list'].split(',') if candidate['mutation_list'] else []
        ninterface=sum(int(m[1:-1]) in interface for m in mutations)
        columns.update(interface_mutation_count=ninterface,scaffold_mutation_count=len(mutations)-ninterface,
                       selection_rationale=candidate['selection_rule'],**bio)
        for prefix,values in [('discovery',discovery),('heldout',h)]:
            columns.update({prefix+'_'+k:v for k,v in values.items() if k not in {'candidate_id','selection_rule','genotype_stratified'}})
        columns['warning']=json.loads((ROOT/'results/benchmark_summary.json').read_text())['warning']
        flat.append(columns)
        generalization.append(dict(candidate_id=candidate['candidate_id'],selection_rule=candidate['selection_rule'],
            discovery_worst_delta_wt=discovery['worst_delta_wt'],heldout_worst_delta_wt=h['worst_delta_wt'],
            worst_delta_change=float(h['worst_delta_wt'])-float(discovery['worst_delta_wt']),
            discovery_median_delta_wt=discovery['median_delta_wt'],heldout_median_delta_wt=h['median_delta_wt'],
            median_delta_change=float(h['median_delta_wt'])-float(discovery['median_delta_wt']),
            interpretation='Different state panels; descriptive generalization change, not a biological effect'))
    write_tsv(ROOT/'results/candidate_summary.tsv',flat)
    write_tsv(ROOT/'results/strategy_generalization.tsv',generalization)
    comparisons=[]
    for strategy in ['single_state','escape_aware']:
        rows=[r for r in held.values() if strategy in r['selection_rule']]
        comparisons.append(dict(strategy=strategy,number_candidates=len(rows),
            mean_candidate_worst_delta_wt=np.mean([float(r['worst_delta_wt']) for r in rows]),
            mean_candidate_median_delta_wt=np.mean([float(r['median_delta_wt']) for r in rows]),
            largest_candidate_worst_delta_wt=max(float(r['worst_delta_wt']) for r in rows),
            mean_fraction_better_wt=np.mean([float(r['fraction_better_wt']) for r in rows]),
            mean_fraction_within_wt_tolerance=np.mean([float(r['fraction_within_wt_tolerance']) for r in rows]),
            interpretation='Equal frozen candidates and held-out states; overlapping selections retained'))
    write_tsv(ROOT/'results/heldout_strategy_comparison.tsv',comparisons)
    model=read_tsv(ROOT/'results/model_comparison.tsv'); seeds=read_tsv(ROOT/'results/seed_variance.tsv')
    effects=[]
    for arm in ['A','B']:
        for source,metrics in [(model,['mean_hydrophobic_sasa_fraction','mean_hydrophobic_patch_proxy','mean_stability_score','mean_charge','mean_pI']),
                               (seeds,['mean_hydrophobic_fraction','mean_identity','mean_nll'])]:
            for metric in metrics:
                standard={r['seed']:float(r[metric]) for r in source if r['design_arm']==arm and r['model']=='standard'}
                soluble={r['seed']:float(r[metric]) for r in source if r['design_arm']==arm and r['model']=='soluble'}
                paired=sorted(standard.keys()&soluble.keys()); delta=np.array([soluble[s]-standard[s] for s in paired])
                rng=np.random.default_rng(config()['bootstrap_seed'])
                boots=delta[rng.integers(0,len(delta),(config()['bootstrap_replicates'],len(delta)))].mean(axis=1)
                low,high=np.quantile(boots,[.025,.975])
                effects.append(dict(design_arm=arm,metric=metric,seeds=','.join(paired),number_paired_seeds=len(paired),
                    standard_seed_mean=np.mean([standard[s] for s in paired]),soluble_seed_mean=np.mean([soluble[s] for s in paired]),
                    soluble_minus_standard=float(delta.mean()),ci95_low=float(low),ci95_high=float(high),
                    bootstrap_unit='matched generation seed within arm; temperatures pooled within seed',
                    bootstrap_replicates=config()['bootstrap_replicates'],bootstrap_seed=config()['bootstrap_seed'],
                    caveat='Small seed sample; descriptors on one scaffold; not expression or solubility measurements'))
    write_tsv(ROOT/'results/model_comparison_effects.tsv',effects)

def score_dictionary():
    import Bio
    sw=read_tsv(ROOT/'results/software_manifest.tsv')
    version={r['tool']:r['version'] for r in sw}
    definitions=[
      ('ProteinMPNN_score','ProteinMPNN','Mean -ln P(sequence residue | backbone, decoding order, preceding sequence) over masked redesigned positions; upstream _scores and mask*chain_M*chain_M_pos','lower','unitless, natural-log NLL','Official v_48_020 standard or soluble training','Model sequence compatibility','Does not demonstrate binding, affinity, folding or solubility'),
      ('interaction_score','EvoEF2','ComputeBinding split=E,R: weighted complex energy minus separated-chain contributions with coordinates unchanged','lower','EvoEF2 score units; physical unit not asserted','EvoEF2 design parameterization and rotamer libraries','Compare this protein-only fixed-conformation model output','Does not estimate KD or rigorous binding free energy; failed mutant gate'),
      ('stability_score','EvoEF2','ComputeStability on extracted receptor monomer; weighted total terms including residue reference terms','lower','EvoEF2 score units','EvoEF2 design parameterization','Model total-energy descriptor','Not measured melting point, expression, folding probability or stability'),
      ('sasa_angstrom2','biopython','Shrake-Rupley atom sphere sampling, probe 1.4 A, 100 sphere points; monomer residue sum','descriptor','angstrom squared','Geometric solvent probe model','Solvent-exposure geometry','Not measured solubility'),
      ('hydrophobic_sasa_fraction','biopython','Sum monomer residue SASA for AVILMFWY / total receptor SASA','descriptor','unitless','Predeclared hydrophobic alphabet and geometry','Exposure proxy','Not measured solubility or aggregation'),
      ('hydrophobic_patch_proxy_angstrom2','biopython','Maximum sum hydrophobic-residue SASA within 8 A CA sphere around a receptor residue','descriptor','angstrom squared','Arbitrary local geometric neighborhood','Patch-size proxy','Not experimentally measured hydrophobic patch or aggregation'),
      ('net_charge_ph7','biopython','ProteinAnalysis charge_at_pH(7), residue and terminal Henderson-Hasselbalch model','descriptor','elementary-charge proxy','Biopython pKa parameters','Composition descriptor','Not measured electrostatics of the folded molecule'),
      ('pI','biopython','ProteinAnalysis.isoelectric_point, bisection of charge model','descriptor','pH','Biopython pKa parameters','Composition descriptor','Not measured pI'),
      ('hydrophobic_fraction','pipeline','Count AVILMFWY / sequence length','descriptor','unitless','Predeclared alphabet','Composition descriptor','Not measured solubility'),
      ('sequence_identity_percent','pipeline','100 * exact residue matches to human scaffold / scaffold length','descriptor','percent','Human 3X0E observed segment mapped to P60033','Sequence novelty relative to this scaffold','Not candidate quality'),
      ('mutation_count','pipeline','Number of residues differing from WT observed human scaffold','descriptor','count','Exact residue comparison','Sequence change count','Not improved function'),
      ('longest_unchanged_stretch','pipeline','Longest contiguous sequence-index run identical to WT','descriptor','residue count','Exact residue comparison','Sequence novelty descriptor','Not fold similarity'),
      ('shannon_entropy','pipeline','-sum p(aa)*log2 p(aa); excludes gaps and ambiguous aa; accession-weighted discovery only','higher = more diverse','bits','Bounded GenBank convenience sample','Observed sequence diversity','Not population prevalence or escape phenotype'),
      ('clash_count','pipeline','Number of interchain heavy-atom pairs at distance <2.0 A','lower','atom-pair count','Arbitrary fixed geometric threshold','Geometric incompatibility descriptor','Not energetic penalty or experimentally observed clash'),
      ('contact_count','pipeline','Number of interchain heavy-atom pairs at distance <=5.0 A','descriptor','atom-pair count','Primary cutoff fixed before candidate scores','Geometric interface descriptor','Not independent contacts or biological affinity'),
      ('ca_rmsd_angstrom','biopython','Square root mean squared distance after least-squares rigid CA superposition','lower','angstrom','Observed sequence-matched CA atoms','Experimental conformation mismatch','Not model confidence or fold validation'),
      ('median_score','pipeline','Median EvoEF2 interaction score across equally weighted structural states','lower','EvoEF2 score units','State aggregation; no composite across unlike metrics','Model central tendency','Not binding across genotypes'),
      ('worst_score','pipeline','Maximum EvoEF2 interaction score across equally weighted structural states','lower','EvoEF2 score units','State aggregation','Model worst-state behavior','Not worst biological phenotype'),
      ('best_score','pipeline','Minimum EvoEF2 interaction score across states','lower','EvoEF2 score units','State aggregation','Model best-state behavior','Not biological performance'),
      ('paired_delta_wt','pipeline','Candidate interaction score minus WT score on same viral state','lower','EvoEF2 score units','Identical state and preparation protocol','Paired model comparison','Not experimental mutation thermodynamics'),
      ('worst_delta_wt','pipeline','Maximum paired candidate-minus-WT delta across states','lower','EvoEF2 score units','Predeclared escape-aware lexicographic selection','Computational prioritization','Not experimentally established escape resistance'),
      ('median_delta_wt','pipeline','Median paired candidate-minus-WT delta across states','lower','EvoEF2 score units','State aggregation','Paired model summary','Not affinity'),
      ('iqr','pipeline','75th minus 25th percentile interaction scores across states','descriptor','EvoEF2 score units','State aggregation','Model variability','Not biological uncertainty'),
      ('std','pipeline','Population standard deviation of interaction scores across states','descriptor','EvoEF2 score units','State aggregation','Model variability','Not independent experiment error'),
      ('fraction_better_wt','pipeline','Fraction of equally weighted states with paired delta <0','higher','unitless','Comparison with native human scaffold','Model comparison','Not binding frequency'),
      ('fraction_within_wt_tolerance','pipeline','Fraction states with paired delta <=0.5 arbitrary score tolerance','higher','unitless','Engineering tolerance fixed before selection','Model comparison','Not biological noninferiority'),
      ('rank_spearman','scipy','Spearman correlation of discovery and held-out worst paired-delta ranks across frozen candidates','higher','unitless','Candidate-level rank comparison','Rank stability in this sample','Not quantitative experimental validation'),
      ('worst_delta_effect_escape_minus_single','pipeline','Mean candidate maximum paired WT delta for escape set minus mean maximum for single set','lower','EvoEF2 score units','Equal candidates per strategy; shared candidates retained','Between-selection-strategy model effect','Not biological improvement'),
    ]
    rows=[dict(score_name=n,tool=t,version=version.get(t,'original code 0.1.0'),mathematical_or_algorithmic_definition=d,
               direction=direction,units_or_unitless=u,training_or_parameterization_context=context,what_it_can_support=yes,what_it_cannot_support=no) for n,t,d,direction,u,context,yes,no in definitions]
    write_tsv(ROOT/'results/score_dictionary.tsv',rows)
    doc=['# Score dictionary','', 'All untested sequences are computational candidates. The experimental mutation gate failed in this run.','']
    for row in rows:
        doc.extend(['## '+row['score_name'],'',row['mathematical_or_algorithmic_definition']+'. '+row['units_or_unitless']+'.',
                    '',row['what_it_can_support']+'. '+row['what_it_cannot_support']+'.',''])
    (ROOT/'docs').mkdir(exist_ok=True); (ROOT/'docs/score_dictionary.md').write_text('\n'.join(doc),encoding='utf-8')

def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
    folder=ROOT/'figures'; folder.mkdir(exist_ok=True)
    manifest=[]
    def save(fig,number,title,finding):
        filename=f'figure_{number}.png'; fig.tight_layout(); fig.savefig(folder/filename); fig.savefig(folder/f'figure_{number}.pdf'); plt.close(fig)
        manifest.append(dict(figure=number,file='figures/'+filename,title=title,finding=finding,generation_stage='figures',sha256=sha256(folder/filename)))
    from .structures import read_structure,interface
    model=read_structure(ROOT/'data/processed/humanized_unrepaired.pdb')
    fig=plt.figure(figsize=(10,4)); ax=fig.add_subplot(121,projection='3d')
    for chain,color,label in [('E','#547aa5','Experimental E2 orientation'),('R','#bb6b48','Aligned human CD81')]:
        xyz=np.array([r['CA'].coord for r in model[0][chain] if 'CA' in r]); ax.plot(*xyz.T,color=color,lw=1,label=label)
    from Bio.PDB import MMCIFParser
    metadata=json.loads((ROOT/'data/processed/humanization.json').read_text())
    source=MMCIFParser(QUIET=True,auth_chains=False).get_structure('7MWX',str(ROOT/'data/raw/7MWX.cif'))
    xyz=np.array([r['CA'].coord for r in source[0][metadata['tamarin_chain']] if 'CA' in r]); ax.plot(*xyz.T,color='#839788',lw=1,label='Experimental bound tamarin CD81')
    ax.set(xlabel='x (A)',ylabel='y (A)',zlabel='z (A)',title='Modeled hybrid, experimental backbones'); ax.legend(fontsize=7)
    ax=fig.add_subplot(122); human=read_tsv(ROOT/'results/humanization_metrics.tsv')[0]; repair=read_tsv(ROOT/'results/repair_audit.tsv')[0]
    ax.bar(['Original tamarin contacts','Hybrid contacts','Repaired contacts'],[int(human['contacts_before']),int(human['contacts_after']),int(repair['contacts_after'])],color=['#547aa5','#bb6b48','#839788'])
    ax.set(ylabel='Heavy atom pairs at 5 A',title=f'CA RMSD {float(human["ca_rmsd_angstrom"]):.2f} A; initial clashes {human["initial_steric_clashes"]}')
    save(fig,1,'Humanization and interface audit','Substantial conformation mismatch and residual repaired clashes limit interpretation.')
    con=[r for r in read_tsv(ROOT/'results/e2_position_conservation.tsv') if r['is_interface']=='True']; fig,(ax,axg)=plt.subplots(2,1,figsize=(10,6))
    ax.bar(range(len(con)),[float(r['shannon_entropy']) for r in con],color='#547aa5'); ax.set_xticks(range(len(con)),[r['h77_position'] for r in con],rotation=60)
    ax.set(xlabel='H77 E2 position (discovery only)',ylabel='Shannon entropy (bits)',title='Geometric interface diversity in a convenience sample')
    genotypes=sorted(json.loads(con[0]['genotype_specific_frequencies']))
    matrix=[]; labels=[]
    for genotype in genotypes:
        counts=[json.loads(r['genotype_specific_frequencies'])[genotype] for r in con]
        values=[]
        for count in counts:
            total=sum(v for aa,v in count.items() if aa in 'ACDEFGHIKLMNPQRSTVWY')
            values.append(-sum((v/total)*np.log2(v/total) for aa,v in count.items() if aa in 'ACDEFGHIKLMNPQRSTVWY') if total else 0)
        matrix.append(values); labels.append(genotype+' (n='+str(sum(counts[0].values()))+')')
    im=axg.imshow(matrix,aspect='auto',cmap='viridis',vmin=0)
    axg.set_yticks(range(len(labels)),labels); axg.set_xticks(range(len(con)),[r['h77_position'] for r in con],rotation=60)
    axg.set(xlabel='H77 position',ylabel='Annotated genotype',title='Within-genotype entropy; unknown annotations retained'); fig.colorbar(im,ax=axg,label='Bits')
    save(fig,2,'Discovery E2 interface conservation','Some interface positions vary in the retrieved discovery sequences; frequencies are sampling-dependent.')
    bench=read_tsv(ROOT/'results/cd81_mutation_benchmark.tsv'); unique={}
    for r in bench:
        if 'entry' not in r['assay_context']:
            unique[(r['mutation'],r['experimental_direction'])]=r
    b=list(unique.values()); fig,ax=plt.subplots(figsize=(10,4)); colors={'enhanced':'#3c8c73','reduced':'#bf604b','approximately_neutral':'#838383'}
    for i,r in enumerate(b):
        ax.scatter(i,float(r['delta_score']),c=colors[r['experimental_direction']],s=55,marker='o' if r['concordant']=='true' else 'x')
    ax.axhline(0,color='black',lw=.6); ax.set_xticks(range(len(b)),[r['mutation'] for r in b],rotation=60); ax.set(ylabel='Mutant minus WT EvoEF2 score',title='Experimental direction is categorical; x = disagreement')
    for d,c in colors.items(): ax.scatter([],[],c=c,label=d.replace('_',' '))
    ax.legend(fontsize=7); save(fig,3,'Experimental mutation score benchmark','The interface scoring gate failed; candidate model outputs cannot support improved-binding claims.')
    metrics=read_tsv(ROOT/'results/discovery_metrics.tsv'); generated={r['candidate_id']:r for r in read_tsv(ROOT/'results/generated_sequences.tsv.gz')}
    fig,ax=plt.subplots(figsize=(6,4))
    for arm,color in [('A','#547aa5'),('B','#bb6b48')]:
        rows=[r for r in metrics if generated[r['candidate_id']]['design_arm']==arm]
        ax.scatter([float(r['reference_score']) for r in rows],[float(r['worst_delta_wt']) for r in rows],s=10,alpha=.6,c=color,label='Arm '+arm)
    ax.set(xlabel='Reference-state interaction score',ylabel='Worst discovery candidate minus WT score'); ax.legend()
    frozen=read_tsv(ROOT/'results/candidate_freeze.tsv')
    for strategy,marker,color in [('single_state','o','#16365c'),('escape_aware','s','#267858')]:
        selected={r['candidate_id'] for r in frozen if strategy in r['selection_rule']}
        rows=[r for r in metrics if r['candidate_id'] in selected]
        ax.scatter([float(r['reference_score']) for r in rows],[float(r['worst_delta_wt']) for r in rows],s=70,facecolors='none',edgecolors=color,marker=marker,label=strategy.replace('_',' '))
    ax.legend(fontsize=8)
    save(fig,4,'Single-state and multi-state discovery tradeoff','Single-state score and worst paired delta measure different model properties.')
    frozen=read_tsv(ROOT/'results/candidate_freeze.tsv'); scores=read_tsv(ROOT/'results/discovery_scores.tsv.gz')
    ids=['wild_type']+[r['candidate_id'] for r in frozen]; stateids=sorted({r['haplotype_id'] for r in scores}); lookup={(r['candidate_id'],r['haplotype_id']):float(r['interaction_score']) for r in scores if r['score_status']=='passed'}
    matrix=np.array([[lookup[c,s]-lookup['wild_type',s] for s in stateids] for c in ids]); fig,ax=plt.subplots(figsize=(11,4)); limit=max(abs(matrix.min()),abs(matrix.max()),1)
    im=ax.imshow(matrix,aspect='auto',cmap='coolwarm',vmin=-limit,vmax=limit); ax.set_yticks(range(len(ids)),ids); ax.set_xticks(range(len(stateids)),[s[-6:] for s in stateids],rotation=70); fig.colorbar(im,ax=ax,label='Paired EvoEF2 delta versus WT')
    save(fig,5,'Frozen candidates across discovery E2 states','WT is the zero baseline on every state; positive deltas favor WT in the model.')
    held=read_tsv(ROOT/'results/heldout_scores.tsv'); fm=read_tsv(ROOT/'results/heldout_metrics.tsv'); fig,ax=plt.subplots(figsize=(10,4)); rng=np.random.default_rng(42)
    for i,row in enumerate(fm):
        deltas=[float(r['paired_delta_wt']) for r in held if r['candidate_id']==row['candidate_id'] and r['score_status']=='passed']
        color='#3c8c73' if 'escape_aware' in row['selection_rule'] else '#547aa5'
        ax.scatter(i+rng.uniform(-.12,.12,len(deltas)),deltas,s=22,alpha=.65,c=color)
        ax.scatter(i,max(deltas),s=70,marker='_',c='black')
    ax.axhline(0,color='black',lw=.7); ax.set_xticks(range(len(fm)),[r['candidate_id']+'\n'+r['selection_rule'].replace(',','/').replace('_',' ') for r in fm],rotation=55,fontsize=7)
    ax.set(ylabel='Held-out paired EvoEF2 delta versus WT',title='Individual held-out states; black tick = candidate worst state')
    save(fig,6,'Frozen held-out evaluation','The headline effect is a model-score comparison under a failed experimental validation gate.')
    bio={r['candidate_id']:r for r in read_tsv(ROOT/'results/candidate_biophysics.tsv')}; fig,ax=plt.subplots(figsize=(6,4))
    for mode,color in [('standard','#547aa5'),('soluble','#bb6b48')]:
        rows=[r for c,r in generated.items() if c in bio and r['ProteinMPNN_model']==mode]
        ax.scatter([float(r['sequence_identity_percent']) for r in rows],[float(bio[r['candidate_id']]['hydrophobic_sasa_fraction']) for r in rows],s=10,alpha=.6,label=mode,c=color)
    ax.set(xlabel='Identity to human scaffold (%)',ylabel='Hydrophobic monomer SASA fraction'); ax.legend()
    save(fig,7,'Sequence and exposure tradeoff','Soluble-model differences are geometric proxies, not expression or solubility measurements.')
    latest={r['stage']:r for r in read_tsv(ROOT/'logs/resource_usage.tsv') if r['exit_status']=='0'}
    resources=list(latest.values()); fig,axs=plt.subplots(1,3,figsize=(13,6)); names=[r['stage'] for r in resources]
    for ax,column,label in zip(axs,['elapsed_seconds','peak_rss_bytes','disk_peak_bytes'],['Time (s)','Peak sampled process-tree RSS (GB)','Peak project bytes (GB)']):
        scale=1 if column=='elapsed_seconds' else 1e9; ax.barh(range(len(names)),[float(r[column])/scale for r in resources],color='#547aa5'); ax.set_yticks(range(len(names)),names,fontsize=7); ax.set_xlabel(label)
    save(fig,8,'Measured resource profile','Scientific stages and a fresh locked CPU installation are measured; original bootstrap telemetry remains unavailable.')
    write_tsv(ROOT/'results/figure_manifest.tsv',manifest)
    captions=['# Scripted figures','', 'Regenerate with `python -m pipeline.cli figures` after held-out evaluation. PNG files support README viewing; PDF files support export.','']
    for row in manifest:
        captions.extend(['## Figure '+str(row['figure'])+'. '+row['title'],'',row['finding'],'',f'![{row["title"]}](figure_{row["figure"]}.png)',''])
    (folder/'README.md').write_text('\n'.join(captions),encoding='utf-8')

def report():
    score_dictionary()
    benchmark=json.loads((ROOT/'results/benchmark_summary.json').read_text()); diversity=json.loads((ROOT/'results/diversity_summary.json').read_text()); held=json.loads((ROOT/'results/heldout_summary.json').read_text())
    human=read_tsv(ROOT/'results/humanization_metrics.tsv')[0]; pilot=json.loads((ROOT/'results/mpnn_pilot.json').read_text()); resources=read_tsv(ROOT/'logs/resource_usage.tsv')
    generated=read_tsv(ROOT/'results/generated_sequences.tsv.gz'); funnel=read_tsv(ROOT/'results/design_filter_funnel.tsv'); panel=read_tsv(ROOT/'data/processed/discovery_panel.tsv')
    candidate_reporting()
    peakram=max(int(r['peak_rss_bytes']) for r in resources); peakdisk=max(int(r['disk_peak_bytes']) for r in resources); runtime=sum(float(r['elapsed_seconds']) for r in resources if r['exit_status']=='0' and r['stage'] not in {'scoring_pilot'})
    design_time=float([r for r in resources if r['stage']=='design' and r['exit_status']=='0'][0]['elapsed_seconds'])
    summary=dict(benchmark=benchmark,diversity=diversity,heldout=held,humanization_rmsd=float(human['ca_rmsd_angstrom']),generated=len(generated),
                 unique_passed=next(int(r['count']) for r in funnel if r['transition']=='unique_passed'),peak_rss_bytes=peakram,peak_project_bytes=peakdisk,
                 stage_elapsed_sum_seconds=runtime,discovery_panel_coverage=float(panel[0]['achieved_coverage']),failure_count=len(read_tsv(ROOT/'logs/failures.tsv')))
    summary.update(generation_elapsed_seconds=design_time,generation_projection_error_percent=100*(design_time/pilot['projected_elapsed_seconds']-1))
    write_json(ROOT/'results/run_summary.json',summary)
    issues=[dict(issue='Experimental directional scoring gate failed',evidence=f'{benchmark["concordant_directional"]}/{benchmark["directional_binding_controls"]} unique directional controls; F186L effectively neutral',
                 why_it_matters='Known receptor effects are not reliably recovered',effect_on_claims='All candidate results are model outputs; no improved-binding interpretation'),
            dict(issue='Humanization conformation mismatch',evidence=f'CA RMSD {human["ca_rmsd_angstrom"]} A; initial clashes {human["initial_steric_clashes"]}',why_it_matters='Scoring may reflect scaffold placement artifacts',effect_on_claims='Fixed hybrid comparison only'),
            dict(issue='Compact discovery panel misses target',evidence=f'{summary["discovery_panel_coverage"]} achieved at 24 states versus target 0.9',why_it_matters='Rare observed haplotypes are underrepresented',effect_on_claims='No broad population coverage claim'),
            dict(issue='Cluster split observation imbalance',evidence=f'{diversity["held_out"]}/{diversity["included"]} held-out observations',why_it_matters='Near-haplotype clusters have unequal size',effect_on_claims='Nominal cluster split is not an 80/20 accession split'),
            dict(issue='Generation timing projection underestimated actual duration',evidence=f'Projected {pilot["projected_elapsed_seconds"]} s; observed {design_time} s',
                 why_it_matters='Repeated model loads and concurrent installation may contribute; effects were not isolated',effect_on_claims='Report actual duration; throughput cap was projected, not a guaranteed wall-clock bound'),
            dict(issue='Bootstrap installation resource telemetry gap',evidence='Initial bootstrap preceded instrumentation; canonical locked fresh installation later reproduced with process-tree and project telemetry',why_it_matters='Original bootstrap peak cannot be reconstructed',effect_on_claims='Claims use measured scientific stages and canonical reproduced installation, not an invented original peak')]
    write_tsv(ROOT/'results/scientific_discomfort.tsv',issues)
    trace=[]
    def number(claim,value,source,column):
        trace.append(dict(readme_claim=claim,value=value,source_file=source,source_column_or_metric=column,generation_stage='report'))
        return str(value)
    correct=number('Directional controls recovered',benchmark['concordant_directional'],'results/benchmark_summary.json','concordant_directional')
    controls=number('Distinct directional controls',benchmark['directional_binding_controls'],'results/benchmark_summary.json','directional_binding_controls')
    rows=number('Experimental assay observations',benchmark['assay_rows'],'results/benchmark_summary.json','assay_rows')
    natural=number('Retrieved GenBank records',diversity['retrieved'],'results/diversity_summary.json','retrieved')
    eligible=number('Eligible natural records',diversity['included'],'results/diversity_summary.json','included')
    excluded=number('Excluded natural records',diversity['excluded'],'results/diversity_summary.json','excluded')
    ndesign=number('Generated computational candidates',len(generated),'results/run_summary.json','generated')
    nunique=number('Unique passing candidates',summary['unique_passed'],'results/run_summary.json','unique_passed')
    rmsd=number('Humanization CA RMSD',f'{float(human["ca_rmsd_angstrom"]):.2f}','results/humanization_metrics.tsv','ca_rmsd_angstrom')
    effect=number('Held-out worst paired-delta effect',f'{held["worst_delta_effect_escape_minus_single"]:.2f}','results/heldout_summary.json','worst_delta_effect_escape_minus_single')
    ci=number('Held-out percentile interval',f'[{held["ci95"][0]:.2f}, {held["ci95"][1]:.2f}]','results/heldout_summary.json','ci95')
    baseline=number('Majority-direction baseline controls',benchmark['majority_direction_baseline_correct'],'results/benchmark_summary.json','majority_direction_baseline_correct')
    direction_ci=number('Directional concordance interval',f'[{benchmark["concordance_ci95"][0]:.3f}, {benchmark["concordance_ci95"][1]:.3f}]','results/benchmark_summary.json','concordance_ci95')
    interpretation='Escape-aware selection lowered the held-out worst paired model score.' if held['worst_delta_effect_escape_minus_single']<0 else 'Escape-aware selection did not lower the held-out worst paired model score.'
    model_effects=read_tsv(ROOT/'results/model_comparison_effects.tsv')
    model_text=[]
    for arm in ['A','B']:
        item=next(r for r in model_effects if r['design_arm']==arm and r['metric']=='mean_hydrophobic_sasa_fraction')
        delta=number('Arm '+arm+' soluble-minus-standard exposure difference',f'{float(item["soluble_minus_standard"]):.4f}','results/model_comparison_effects.tsv','soluble_minus_standard; arm='+arm+'; metric=mean_hydrophobic_sasa_fraction')
        interval=number('Arm '+arm+' exposure seed-bootstrap interval',f'[{float(item["ci95_low"]):.4f}, {float(item["ci95_high"]):.4f}]','results/model_comparison_effects.tsv','ci95_low/ci95_high; arm='+arm)
        model_text.append(f'Arm {arm}: {delta}, paired seed-bootstrap interval {interval}.')
    model_sentence=' '.join(model_text)
    hstates=number('Held-out E2 states',held['states'],'results/heldout_summary.json','states')
    coverage=number('Achieved discovery coverage percent',f'{summary["discovery_panel_coverage"]*100:.1f}','results/run_summary.json','discovery_panel_coverage')
    dseq=number('Discovery accessions',diversity['discovery'],'results/diversity_summary.json','discovery')
    hseq=number('Held-out accessions',diversity['held_out'],'results/diversity_summary.json','held_out')
    failures=number('Recorded failures',summary['failure_count'],'results/run_summary.json','failure_count')
    ram=number('Measured scientific-stage RSS GB',f'{peakram/1e9:.3f}','results/run_summary.json','peak_rss_bytes')
    disk=number('Measured peak project GB',f'{peakdisk/1e9:.3f}','results/run_summary.json','peak_project_bytes')
    seconds=number('Summed successful stage time seconds',f'{runtime:.1f}','results/run_summary.json','stage_elapsed_sum_seconds')
    generation_seconds=number('Observed generation duration',f'{design_time:.1f}','results/run_summary.json','generation_elapsed_seconds')
    projected_generation=number('Projected generation duration',f'{pilot["projected_elapsed_seconds"]:.1f}','results/mpnn_pilot.json','projected_elapsed_seconds')
    allvalues={'Pilot count':(pilot['count'],'results/mpnn_pilot.json','count'),'Pilot seconds':(round(pilot['elapsed_seconds'],1),'results/mpnn_pilot.json','elapsed_seconds'),
               'Discovery state limit':(24,'results/diversity_summary.json','discovery_states'),'Master seed':(config()['master_seed'],'results/smoke_test.json','seed'),
               'Primary interface cutoff':(config()['interface_cutoff'],'results/e2_interface_residues.tsv','cutoff_angstrom'),
               'Bootstrap replicates':(held['bootstrap_replicates'],'results/heldout_summary.json','bootstrap_replicates'),
               'Frozen count':(held['frozen'],'results/heldout_summary.json','frozen')}
    for claim,(value,source,column) in allvalues.items(): number(claim,value,source,column)
    for claim,value,column in [('Command thread setting',2,'threads'),('Command RAM GB',14,'ram_bytes / 1e9'),('Command disk ceiling GB',13,'requested cap'),('Command seed',20261004,'seed')]:
        number(claim,value,'results/configuration.tsv',column)
    number('Run-count default',1,'results/repair_audit.tsv','num_of_runs')
    for n in range(1,9): number('Figure identifier',n,'results/figure_manifest.tsv','figure')
    text=f'''# EvoEF2 fails the CD81 directional benchmark in an escape-aware sequence-design experiment

## Summary

The interface scoring procedure recovered {correct} of {controls} distinct directional soluble-E2 mutation controls. The curated evidence contains {rows} assay observations. This failed the predeclared directional validation gate, so every designed sequence remains a computational candidate and all structural scores are model outputs.

The run retrieved {natural} GenBank records and retained {eligible}. Official standard and soluble ProteinMPNN generated {ndesign} sequences on the experimental human CD81 scaffold; {nunique} unique sequences passed sequence constraints. Escape-aware selection changed the mean candidate worst paired score by {effect} EvoEF2 score units relative to single-state selection across {hstates} held-out E2 states. The paired state-bootstrap interval was {ci}. {interpretation} Lower scores are favored within this model. The score benchmark failure prevents interpreting this effect as better binding.

## Background

HCV E2 engages the large extracellular loop of CD81 during entry. A soluble receptor-derived construct is a receptor-mimic hypothesis, while membrane CD81 also participates in later entry events. Published soluble-E2 binding and pseudoparticle entry effects can disagree. Inverse folding samples amino-acid sequences conditioned on coordinates. This experiment performs constrained de novo sequence design on an experimental backbone; ProteinMPNN does not generate a new fold here.

Natural E2 variation creates a multi-state selection problem. The comparison uses one generated pool and two transparent selection rules. Single-state selection minimizes the reference score. Escape-aware selection minimizes the largest candidate-minus-WT score across discovery states, then its median. WT is scored under the same state geometry. Unlike a weighted sum, these rules do not mix NLL, charge and exposure into one arbitrary quantity.

## Data

The experimental complex [7MWX](https://www.rcsb.org/structure/7MWX) contains Saguinus oedipus (tamarin) CD81. Homo sapiens in its expression-system metadata does not identify the receptor species. The human scaffold is [3X0E](https://www.rcsb.org/structure/3X0E). The human E2:CD81 complex is a modeled hybrid constructed from these two experimental structures. The CA alignment RMSD was {rmsd} A. Original downloaded mmCIF files are kept untouched and excluded from Git; checksum and chain/entity audits are tracked.

Primary ground truth comes from Higginbottom, Drummer, Bertaux and Dragic, and Flint. The supplied attribution for the Different domains paper was corrected using its primary record. Assay directions remain qualitative; no plotted values were digitized. See [ground truth](config/cd81_mutation_ground_truth.tsv) for source-specific construct, assay, strain uncertainty, replicate information and excluded secondary summaries.

The bounded GenBank sample excluded {excluded} records for translation, ambiguity, interface-coverage or modeling failures. H77 AF009606.1 provides one-based polyprotein numbering. Exact E2 sequences retain accession multiplicity. The cluster split assigned {dseq} discovery and {hseq} held-out accessions. Connected interface haplotypes differing at no more than one site remain together. Unequal connected-component sizes prevent an accession-balanced split. The compact discovery panel represented {coverage}% of discovery observations, below the arbitrary target. Unknown genotypes remain explicitly unclassified.

## Pipeline

Each stage is independently callable through the CLI. Retrieval stores original bytes and checksums. Humanization uses sequence matching and rigid CA superposition, followed by lightweight sidechain repair. Natural sequences are aligned to H77 using Biopython affine-gap alignment; observed structure residues are mapped separately. This avoids installing a second aligner. Conservation uses discovery sequences only. The split builder quarantines test sequences before design. A file-open audit guard and explicit discovery-only function inputs enforce the boundary.

ProteinMPNN fixes the E2 sequence and designs only receptor chain R. Arm A preserves the geometric interface and experimental recognition constraints. Arm B permits the remaining reliable contacts to vary. Both preserve native disulfide cysteines and reject additional cysteines. Standard and soluble models receive matched seeds, temperatures and budgets. Effective seeds are verified from upstream output headers. The pilot generated {pilot['count']} sequences in {pilot['elapsed_seconds']:.1f} seconds; its sequences are excluded from the comparison pool. The full budget was fixed from pilot throughput before structural candidate scores. The projection was {projected_generation} seconds, while generation took {generation_seconds} seconds. Repeated loading and concurrent installation may explain part of this difference; their effects were not isolated.

Scoring builds each receptor on the reference complex once and recombines it with separately modeled E2 sidechains. State-specific receptor repacking is omitted. This makes the cross-state conformational assumption explicit. Candidate freeze hashes bind sequences, configuration, software and discovery inputs before held-out evaluation. Pareto exposure/score tradeoffs are annotated separately. The clash tolerance is arbitrary. The freeze retains each strategy's selections, including overlap.

```bash
bash scripts/00_configure.sh --threads 2 --ram 14 --disk 13 --seed 20261004 --yes
bash scripts/02_install.sh
bash run_all.sh --mode core --parallel
bash scripts/19_verify.sh
```

## Results

Experimental validation comes first. EvoEF2 failed the benchmark despite recovering the enhanced T163A direction. The exact binomial interval for directional concordance is {direction_ci}; related assay backgrounds limit its independence interpretation. An always-reduced direction baseline would recover {baseline} of {controls} controls. F186L was effectively neutral under the declared score tolerance, contrary to the soluble-E2 observations. Experimental assay units were not pooled and no continuous experimental correlation was claimed.

![Experimental mutation benchmark](figures/figure_3.png)

The modeled human conformation differs substantially from the bound tamarin conformation. Sidechain repair cannot resolve uncertainty about the backbone orientation.

![Humanization audit](figures/figure_1.png)

Discovery interface entropy varies by H77 position. Sampling and the compact panel limit any statement about circulating diversity.

![Interface conservation](figures/figure_2.png)

Reference-state preference and worst discovery paired scores differ. The two strategies share an identical generation budget and candidate pool.

![Discovery tradeoff](figures/figure_4.png)

The frozen heatmap includes WT as a paired zero comparator in every state. Scores above zero favor WT under the model.

![Frozen discovery states](figures/figure_5.png)

Held-out candidates were evaluated only after the candidate table and its manifest were written. The effect estimate is {effect}, with interval {ci}. {interpretation} Each point is a viral structural state. The small test sample and correlated haplotype components limit inference. The final table carries the failed-benchmark warning next to every computational candidate. [Strategy comparisons](results/heldout_strategy_comparison.tsv) retain both median and worst-case paired deltas, and [generalization changes](results/strategy_generalization.tsv) retain discovery-to-test differences for each frozen computational candidate.

![Held-out evaluation](figures/figure_6.png)

The measured soluble-minus-standard hydrophobic SASA fraction differences were {model_sentence} These compare matched generation seeds within each arm. [Model comparisons](results/model_comparison_effects.tsv) retain seed-level uncertainty for exposure, patch size, composition, charge and total monomer energy. Lower hydrophobic exposure cannot demonstrate expression, folding or solubility. The [candidate summary](results/candidate_summary.tsv) places sequences, mutation lists, geometric metrics and discovery/test statistics in flat columns for inspection.

![Sequence exposure tradeoff](figures/figure_7.png)

The recorded run contains {failures} failures, including a missing runtime for the upstream Windows EvoEF2 executable, unsupported upstream README options, and long-path failure. Local static compilation and verified-source defaults resolved those execution issues. Scientific stages and the reproduced locked installation used peak sampled process-tree RSS of {ram} GB and peak measured project footprint of {disk} GB. Summed successful stage durations were {seconds} seconds; stages that overlapped are not summed wall-clock time.

![Resource profile](figures/figure_8.png)

## Repository structure

`pipeline/` contains the implementation. `scripts/` provides Bash entry points. `config/` records thresholds, ground truth, fixed residues, upstream commits and split assignments. `results/` contains derived tables and freeze hashes. `figures/` contains scripted raster figures and exportable PDFs. `data/processed/` keeps the compact experimental hybrid, states and frozen models. `data/raw/`, `data/work/`, `.venv/` and `vendor/` are ignored. `tests/fixtures/` contains an offline synthetic example. `logs/` records scientific-stage resources and failures. [Score definitions](docs/score_dictionary.md) state what each output can support.

## Usage

Create a local environment with Python and install the locked dependencies. This host required its explicit Python executable because the default command resolved to an older runtime. For the offline smoke test, a CPU tensor runtime, external repositories and model weights are unnecessary.

```powershell
python -m venv .venv
.\\.venv\\Scripts\\python.exe -m pip install numpy biopython scipy matplotlib psutil pyyaml pytest
.\\.venv\\Scripts\\python.exe -m pipeline.cli configure --yes
.\\.venv\\Scripts\\python.exe -m pipeline.cli smoke
.\\.venv\\Scripts\\python.exe -m pytest -q
```

For the scientific run, use Git Bash and the installation command above. The resource guard measures current free disk and protects a non-project reserve. Installed environments, model weights, caches and generated files are counted. The measured footprint was {disk} GB. The original bootstrap preceded instrumentation, so its peak cannot be reconstructed. A canonical locked installation was reproduced in a disposable environment with telemetry and wheel hashes. Independent stages have resource records. Exact environment versions and model hashes are in the manifests.

Git, Bash and a C++ compiler are host prerequisites. Scientific entry points provide `--help`. `--mode smoke` runs bundled analytical fixtures, while `--from heldout` resumes from the immutable freeze. Fresh runs score serially; `--parallel` opts into the configured worker count after the scoring memory pilot. Core and full both honor the pilot-approved compact budget, so full cannot expand past the measured ceiling. Existing outputs are checked before reuse. Atomic table writes preserve completed outputs if a later write fails.

To begin an independent scientific analysis, run `python -m pipeline.new_run --destination data/work/new_analysis`. The helper refuses an existing destination, copies code and predeclared settings, and initializes Git with `.gitignore` first. Set up that directory's own environment and follow the scientific commands there. It retrieves contemporary query results and records its own freeze. The archived result tables in this repository belong to the present run. The exact locked environment targets Python 3.14 on Windows; another platform or compiler requires a separate software audit, and an executable hash mismatch is reported rather than silently accepted.

## Limitations

The model combines an unbound human receptor conformation with the experimentally bound tamarin orientation. The large alignment RMSD is the greatest structural uncertainty. Fixed backbones omit E2 conformational flexibility; only reliably mapped sidechain substitutions are modeled. Coordinates absent from the experiments are not invented. Protein-only scoring removes glycans and other heteroatoms, whose removed identities are recorded. Complete-virion accessibility is not represented.

EvoEF2 interaction outputs are not experimental affinity, KD or rigorous binding free energy. The experimental score gate failed. ProteinMPNN NLL measures conditional sequence compatibility. Stability and SASA outputs do not measure folding, expression or solubility. The literature includes heterogeneous soluble, cellular and entry assays, and some viral strain identifiers are not established in the directly reported methods. The convenience sequence sample and incomplete panel introduce sampling bias. Small, unbalanced held-out sets and connected haplotype dependence weaken the bootstrap interpretation. Genotype-stratified outputs retain unknown and mixed-genotype states.

No generated sequence has experimental binding measurements. Soluble CD81 differs from membrane CD81 in oligomerization, trafficking and later entry biology. Its interactions with other human proteins, immune effects, specificity, safety and pharmacology are unknown. This project uses no molecular dynamics, docking or de novo backbone generation. Original bootstrap telemetry is unavailable, and shellcheck availability is reported by verification.

The least certain conclusion is that the observed held-out strategy score difference would persist under a different receptor conformation. The humanization mismatch and failed experimental directional benchmark allow that difference to be driven by placement and sidechain-packing artifacts.

## Data availability

Public identifiers and accepted NCBI queries are recorded in [datasets](config/datasets.tsv), [download hashes](results/download_manifest.tsv), [sequence exclusions](results/hcv_sequence_manifest.tsv) and [structure provenance](results/structure_provenance.tsv). Cached accession snapshots reproduce this run; fresh default-order query results may change. Each external source can be retrieved through the corresponding CLI stage:

```bash
bash scripts/stage.sh structures
bash scripts/stage.sh humanize
bash scripts/stage.sh ground_truth
bash scripts/stage.sh sequences
bash scripts/stage.sh diversity
```

## Citation

The verified primary citations are [Kumar et al.](https://doi.org/10.1038/s41586-021-03913-5), [Yang et al.](https://doi.org/10.1096/fj.15-272880), [Higginbottom et al.](https://doi.org/10.1128/JVI.74.8.3642-3649.2000), [Drummer et al.](https://doi.org/10.1128/JVI.76.21.11143-11147.2002), [Bertaux and Dragic](https://doi.org/10.1128/JVI.80.10.4940-4948.2006), and [Flint et al.](https://doi.org/10.1128/JVI.00104-06). Software sources, exact commits and usage terms are recorded in [software provenance](results/software_manifest.tsv).

## License

Original code is MIT licensed. Third-party code and weights are retrieved independently and excluded from Git. ProteinMPNN supplies MIT-licensed code and bundled weights; no distinct model-weight license was found, so a separate grant is not assumed. EvoEF2's MIT license file conflicts with its academic-use README language. Its local binary and source are not redistributed. PDB and GenBank usage and article-specific copyright are recorded in [license audit](results/license_audit.tsv). Source article prose and figures are not redistributed.
'''
    # DOI/year, structure and figure identifiers are identifiers, not derived
    # quantitative claims; configuration values are separately traceable.
    write_tsv(ROOT/'results/readme_traceability.tsv',trace)
    (ROOT/'README.md').write_text(text,encoding='utf-8')
    (ROOT/'docs/report.md').write_text(text,encoding='utf-8')
