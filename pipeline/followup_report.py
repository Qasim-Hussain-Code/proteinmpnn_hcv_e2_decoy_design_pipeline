"""Controls, figures, manuscript supplement, and verification of the follow-up."""
from __future__ import annotations
from collections import Counter
import json
import platform
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import ROOT, now, read_tsv, write_tsv, write_json, sha256, footprint
from .followup import OUT, protocol, snapshot, repair, invoke, energy, combine, predicted_direction, mutations
from .structures import read_structure
from Bio.SeqUtils import seq1

def controls():
    declaration=OUT/'technical_controls.json'
    if not declaration.exists():write_json(declaration,dict(timestamp=now(),rationale='Observed E2 repacking sensitivity and genotype 1 undercoverage',wt_control='One additional WT repair, measured both with its repacked E2 and recombined original E2; diagnostic only',coverage_extension='Retain global panel and extend to 90% within every included recorded genotype; no score-dependent choice',posthoc=True))
    scores=read_tsv(OUT/'energy_terms.tsv');lookup={(r['model'],r['mutation'],r['mode']):r for r in scores};truth=read_tsv(OUT/'ground_truth.tsv')
    labels={(r['mutation'],r['normalized_direction']) for r in truth if r['included']=='true' and 'entry' not in r['assay_type'] and r['normalized_direction'] in {'enhanced','reduced'} and '+' not in r['mutation'] and r['followup_role']=='original_evidence'}
    out=[]
    for row in read_tsv(OUT/'model_audit.tsv'):
        ident=row['model'];base=ROOT/row['path'];new=base.with_name('wt_second.pdb');repair(base,new,ident+'_wt_second');fixed=new.with_name('wt_second_fixed.pdb');combine(new,base,fixed)
        baseline=float(lookup[ident,'WT','repacked']['Total'])
        for mode,path in [('repacked',new),('fixed_e2',fixed)]:
            second=energy(path,ident+'_wt_second_'+mode)['Total']
            n=sum((ident,m,mode) in lookup for m,l in labels)
            correct=sum(predicted_direction(float(lookup[ident,m,mode]['Total'])-second,.5)==l for m,l in labels if (ident,m,mode) in lookup)
            out.append(dict(model=ident,mode=mode,one_repair_wt=baseline,two_repair_wt=second,wt_shift=second-baseline,correct_with_second_wt=correct,total=n,interpretation='Additional preparation control, not a new validation gate'))
    write_tsv(OUT/'wt_repair_sensitivity.tsv',out)
    from .followup_diversity import balanced_panel
    balanced_panel()

def figures():
    folder=OUT/'figures';folder.mkdir(exist_ok=True);models=[r['model'] for r in read_tsv(OUT/'model_audit.tsv')]
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'figure.dpi':150})
    def save(fig,name):
        fig.savefig(folder/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(folder/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    scores=read_tsv(OUT/'benchmark_summary.tsv');rows=[r for r in scores if float(r['tolerance'])==.5]
    fig,ax=plt.subplots(figsize=(8,4));x=np.arange(len(models))
    for offset,mode,color in [(-.2,'repacked','#5975a4'),(.2,'fixed_e2','#dd8452')]:
        v={r['model']:r for r in rows if r['mode']==mode};ax.bar(x+offset,[int(v[m]['original_single_correct']) for m in models],.4,label=mode,color=color)
    ax.axhline(10,color='black',ls='--',label='Always-reduced baseline (10/11)');ax.set_ylim(0,11.5);ax.set_xticks(x,models,rotation=30,ha='right');ax.set_ylabel('Correct original directional controls / 11');ax.legend(fontsize=8);ax.set_title('Structural sensitivity does not establish a reliable scorer');save(fig,'benchmark')
    labels=sorted({r['mutation'] for r in read_tsv(OUT/'ground_truth.tsv') if r['included']=='true' and 'entry' not in r['assay_type']})
    data=read_tsv(OUT/'benchmark.tsv');lookup={(r['model'],r['mutation']):float(r['delta']) for r in data if r['mode']=='fixed_e2'}
    matrix=np.array([[lookup.get((m,mut),np.nan) for mut in labels] for m in models]);fig,ax=plt.subplots(figsize=(12,4));im=ax.imshow(matrix,cmap='RdBu_r',vmin=-15,vmax=15,aspect='auto');ax.set_xticks(range(len(labels)),labels,rotation=60,ha='right');ax.set_yticks(range(len(models)),models);fig.colorbar(im,ax=ax,label='Mutant minus WT score (display clipped at +/-15)');ax.set_title('Diagnostic scores under fixed E2 geometry; lower is model-favored');save(fig,'control_deltas')
    audit=read_tsv(OUT/'model_audit.tsv');v={r['model']:r for r in rows if r['mode']=='fixed_e2'};fig,axes=plt.subplots(1,2,figsize=(10,3.6));axes[0].bar(models,[int(r['clashes_after']) for r in audit]);axes[0].set_ylabel('Inter-chain heavy-atom pairs <2 A');axes[1].bar(models,[float(v[m]['wt_vdwrep']) for m in models]);axes[1].set_ylabel('WT inter-chain repulsive score term')
    for ax in axes:ax.tick_params(axis='x',rotation=45)
    fig.suptitle('Rigid fits retain substantial geometry and repulsion uncertainty');save(fig,'geometry')
    fig,ax=plt.subplots(figsize=(7,4))
    for filename,label in [('archival_coverage_curve.tsv','Archived 20-position development panel'),('coverage_curve.tsv','Expanded 37-position development panel')]:
        rows=read_tsv(OUT/filename);ax.plot([int(r['states']) for r in rows],[float(r['coverage']) for r in rows],label=label)
    ax.axhline(.9,color='black',ls='--');ax.axvline(24,color='grey',ls=':');ax.set_xlabel('Number of natural interface haplotype states');ax.set_ylabel('Fraction of sampled development accessions');ax.legend(fontsize=8);ax.set_ylim(0,1.02);save(fig,'coverage')
    original={r['genotype']:r for r in read_tsv(OUT/'genotype_coverage.tsv')};balanced=read_tsv(OUT/'balanced_genotype_coverage.tsv');g=[r['genotype'] for r in balanced];x=np.arange(len(g));fig,ax=plt.subplots(figsize=(7,3.8));ax.bar(x-.2,[float(original[r]['coverage']) for r in g],.4,label='Global 90% panel');ax.bar(x+.2,[float(r['coverage']) for r in balanced],.4,label='Extended panel');ax.axhline(.9,color='black',ls='--');ax.set_xticks(x,g);ax.set_ylabel('Exact haplotype accession coverage');ax.set_xlabel('Recorded genotype (8 excluded by interface indels)');ax.legend(fontsize=8);ax.set_ylim(0,1.08);save(fig,'genotype_coverage')

def report():
    snapshot();figures();coverage=json.loads((OUT/'coverage_summary.json').read_text());balanced=json.loads((OUT/'balanced_panel_summary.json').read_text());models=read_tsv(OUT/'model_audit.tsv');summary=read_tsv(OUT/'benchmark_summary.tsv');bench=read_tsv(OUT/'benchmark.tsv');wt=read_tsv(OUT/'wt_repair_sensitivity.tsv');counts=Counter(r['reason'] for r in read_tsv(OUT/'sequence_manifest.tsv') if r['status']=='excluded');fixed={r['model']:r for r in summary if float(r['tolerance'])==.5 and r['mode']=='fixed_e2'};rep={r['model']:r for r in summary if float(r['tolerance'])==.5 and r['mode']=='repacked'}
    table='\n'.join(f"| {r['model']} | {int(r['clashes_after'])} | {float(fixed[r['model']]['wt_vdwrep']):.2f} | {rep[r['model']]['original_single_correct']}/11 | {fixed[r['model']]['original_single_correct']}/11 | {fixed[r['model']]['external_reduced_correct']}/4 |" for r in models)
    source='https://pmc.ncbi.nlm.nih.gov/articles/PMC111874/'
    text=f'''# Follow-up: scoring remains unreliable; broader coverage is prepared

The requested scoring, structural, and viral-coverage follow-up is complete. The original experiment is preserved byte for byte by [the snapshot](original_snapshot.tsv). No original candidate is re-ranked, and no new designed receptor sequence is generated. These are exploratory diagnostic comparisons with already consulted literature labels.

## Findings

The original model still recovers **5/11** distinct directional binding controls. A human-sequence model on the first bound tamarin backbone recovers **8/11**, but the second experimental bound pose recovers **3/11**. The first model misses all **four** additional binding-loss challenge mutations from [Drummer 2005](https://pubmed.ncbi.nlm.nih.gov/15670777/). Always predicting reduced binding recovers **10/11** original directional labels. Changing the model or the score threshold after seeing these labels cannot establish validation. The failed original experiment remains failed.

The original WT contains a **1.09 A** T163-to-E2 contact after repair and a **93.15** inter-chain van der Waals repulsion term. In the follow-up fixed-E2 T163A comparison, its score delta is **-6.16**, with a repulsive-term delta of **-9.35**. This is consistent with clash removal contributing to the apparent favorable score; it does not prove that the experimentally enhanced interaction has the same mechanism. D196 is **7.70 A** from E2 in the archived model, beyond EvoEF2's 6 A interaction cutoff. Folding, dimerization, receptor dynamics, glycans, membrane context, and assay expression effects are not represented by this interface score.

Sidechain preparation is an additional confound. Repacking both proteins versus recombining every mutant receptor with the exact WT E2 coordinates changes several model/control classifications. WT receives one repair in the initial model, whereas a mutant receives BuildMutant and another repair. [An additional WT repair control](wt_repair_sensitivity.tsv) measures this asymmetry explicitly; it is not a new validation gate. Model-specific shifts and all score terms are retained. Scores are model units, not measured affinities or calibrated free energies.

## Evidence correction

The original curation incorrectly assigned two soluble GST-LEL rows to single F186L and single E188K. [The primary Methods and Results]({source}) describe the **F186L+E188K double mutant** for the soluble construct. The follow-up excludes those two misassigned rows and adds the double-mutant observation. Valid full-length cell-surface single-mutant observations are retained. The corrected original evidence has **28 included assay rows**; adding six source-external challenge observations gives **34**. Correcting this error does not change the original 11 distinct single-mutation directional labels or their 5/11 result. See [the correction audit](curation_audit.tsv) and [the corrected table](ground_truth.tsv). The abstract/Results discrepancy for D196E is retained as assay-context uncertainty. The original files remain an archival record, not the corrected evidence table.

The six new observations are explicitly abstract-level evidence. K124T and V146E retained binding, without establishing WT-equivalent affinity; they are not mislabeled neutral or counted in directional accuracy. The other four report binding loss in recombinant LEL. Full-length and soluble phenotypes are not pooled. Mutations affecting disulfides and dimerization challenge the score's scope.

## Structural comparison

Eight prespecified models use archived 3X0E, both 3X0E chains, both [1G8Q](https://www.rcsb.org/structure/1G8Q) chains, [5TCX](https://www.rcsb.org/structure/5TCX), and both experimental [7MWX](https://www.rcsb.org/structure/7MWX) binding pairs. All human entities and mapped native identities are verified. Missing residues are recorded and are never filled with invented coordinates. The core-fit sensitivity excludes the variable head by the ranges specified in [the locked protocol](protocol.json). Bound-template models use observed tamarin backbone coordinates with five human substitutions and are **not experimentally determined human complexes**. Their missing receptor termini are documented.

| Model | Clashes <2 A | WT repulsive term | Controls, repacked | Controls, fixed E2 | Additional loss controls, fixed E2 |
|---|---:|---:|---:|---:|---:|
{table}

All native WT models retain the two native disulfide geometries. Alternative human rigid fits can worsen clashes markedly; lower CA RMSD alone is not a reliable quality criterion. Zero clashes in a humanized bound-template model do not establish its accuracy or its affinity. No best model is selected using these benchmark outcomes.

![Structural scoring sensitivity](figures/benchmark.png)

## Viral coverage

The bounded retrieval expanded from 309 to **{coverage['retrieved']} unique accessions**, with **{coverage['included']} eligible**, **{coverage['excluded']} excluded**, and **{coverage['unique_e2']} distinct mapped E2 sequences**. The interface definition expands from 20 to **37 H77 positions**, combining the archived human geometry, both native bound interfaces, and the five already curated functional positions. Queries request at most 120 default-order records per stratum and explicitly inspect patent references; the ignored original patent query filter is not reused. Query strings, UIDs, dates, source bytes, and hashes are retained. This is a convenience sample, not a worldwide prevalence estimate.

The original 24-state panel covered **39.44%** of its development accessions; reaching 90% on that same archived cohort requires **{coverage['archival_states_for_90_percent']} states** under the same selection rule. The expanded 37-position development panel requires **{coverage['states_for_90_percent']} states**, covering **{coverage['achieved_coverage']:.2%}** of {coverage['development']} development accessions. Its genotype 1 coverage is only 72.13%. A documented secondary extension to **{balanced['states']} states** raises overall coverage to **{balanced['coverage']:.2%}** and the minimum coverage within each included recorded genotype to **{balanced['minimum_genotype_coverage']:.2%}**. All of these coordinate models are built and hashed. The panel extension uses sequence counts only, with no candidate or assay score selection.

Coverage is conditional on structural eligibility. **Genotype 8 remains outside this coordinate panel**: its four metadata-confirmed records contain interface insertions that a fixed-backbone sidechain-substitution model cannot represent. Those accessions and reasons are retained; no claim of genotype 8 coverage is made. There are **22 unknown-genotype eligible accessions**. Three expanded interface positions, **415-417**, lack coordinates in the archived E2 chain; every state records this missing structural coverage. Thus sequence haplotype coverage at 37 sites and coordinate coverage at the 34 observed sites are distinct. The prepared states are partial structural proxies, not complete models of all natural E2 differences.

The excluded-record counts are: {dict(sorted(counts.items()))}. Strict mapping and indel exclusions can bias the retained cohort, so broad protective efficacy cannot be inferred from these percentages.

A future reserve contains **{coverage['future_reserve']} newly retrieved accessions** in novel connected components, separated from all previously retrieved accessions and all old interface haplotypes within Hamming distance one. These reserve sequences are not modeled, scored, or used to select the panel. The old held-out data are already exposed and remain development data here. A six-accession reserve is too small for a convincing broad-coverage claim; an independent future cohort remains necessary.

![Accession coverage curve](figures/coverage.png)

![Coverage by recorded genotype](figures/genotype_coverage.png)

## Reproduction and artifacts

Use the existing locked CPU environment and pinned EvoEF2 executable. No new dependency stack or GPU predictor is installed. From the project root, run:

```powershell
.\\.venv\\Scripts\\python.exe -m pipeline.followup prepare
.\\.venv\\Scripts\\python.exe -m pipeline.followup models
.\\.venv\\Scripts\\python.exe -m pipeline.followup sequences
.\\.venv\\Scripts\\python.exe -m pipeline.followup coverage
.\\.venv\\Scripts\\python.exe -m pipeline.followup benchmark
.\\.venv\\Scripts\\python.exe -m pipeline.followup states
.\\.venv\\Scripts\\python.exe -m pipeline.followup controls
.\\.venv\\Scripts\\python.exe -m pipeline.followup report
.\\.venv\\Scripts\\python.exe -m pipeline.followup verify
.\\.venv\\Scripts\\python.exe -m pytest -q
```

Online retrieval needs network access. Original raw inputs and the pinned software are prerequisites, as in the original experiment; source caches are outside Git. Retrieval hashes enforce the downloaded snapshot. EvoEF2 receives short local filenames to avoid its upstream path buffer limitation. Separate process work directories prevent simultaneous state and receptor jobs from overwriting inputs. Repair/mutation caches are retained. The initial code hash, final code manifest, and control-addition timestamp distinguish prespecified choices from subsequent technical fixes and coverage extension. This repository adds a supplement; it does not rewrite the original manuscript.

The primary tables are [benchmark results](benchmark.tsv), [energy terms](energy_terms.tsv), [model audits](model_audit.tsv), [sequence QC](sequence_manifest.tsv), [coverage summary](coverage_summary.json), [balanced panel](balanced_development_panel.tsv), and [state hashes](state_manifest.tsv). Five figures are supplied as PNG and PDF. Verification checks archive immutability, model sequence identities, score logs, state identities, coordinate coverage limitations, leakage boundaries, and per-genotype target coverage. See [verification](verification.json) for the measured result. Resource measurements from the original experiment remain untouched; follow-up stage timings printed during execution and the three-state pilot are separate. Earlier follow-up stages did not continuously meter RSS, and no unmeasured peak is claimed.

## Next work

The scoring problem requires a method that represents receptor conformational change and fold integrity, with matched WT/mutant preparation and glycan/complex context, followed by validation on independently curated controls. Receptor-only fold/stability diagnostics and an ensemble constrained by experimental structures are the next computational steps. A change in this score or the score threshold is insufficient. Collecting an independent, geographically broader cohort and modeling interface indels are necessary before genotype 8 can enter a structural breadth claim. Future receptor ranking should use a separately frozen protocol and a fresh independent test set. Experimental binding and soluble-protein behavior would then determine whether any therapeutic-decoy hypothesis is supported.
'''
    (OUT/'README.md').write_text(text,encoding='utf-8')
    write_tsv(OUT/'code_manifest.tsv',[dict(path=str(p.relative_to(ROOT)),sha256=sha256(p)) for p in [ROOT/'pipeline/followup.py',ROOT/'pipeline/followup_diversity.py',ROOT/'pipeline/followup_report.py',ROOT/'tests/test_followup.py',OUT/'protocol.json']])
    outputs=[p for p in OUT.rglob('*') if p.is_file() and p.name not in {'artifact_manifest.tsv','verification.json'} and not p.name.startswith('verification_')]
    write_tsv(OUT/'artifact_manifest.tsv',[dict(path=str(p.relative_to(ROOT)),sha256=sha256(p),bytes=p.stat().st_size) for p in sorted(outputs)])
    print('report and five PDF/PNG figures written',flush=True)

def verify():
    from .followup import lock
    lock()
    checks=[]
    def check(name,passed,detail=''):
        checks.append(dict(check=name,passed=bool(passed),detail=detail))
    snapshot();check('original_archive_hashes',True,len(read_tsv(OUT/'original_snapshot.tsv')))
    for r in read_tsv(OUT/'artifact_manifest.tsv'):check('artifact:'+r['path'],sha256(ROOT/r['path'])==r['sha256'])
    source=json.loads((ROOT/'data/processed/humanization.json').read_text());human=''.join((ROOT/'data/raw/P60033.fasta').read_text().splitlines()[1:])
    for row in read_tsv(OUT/'model_audit.tsv'):
        p=ROOT/row['path'];model=read_structure(p);check('native_wt:'+row['model'],sha256(p)==row['wt_sha256'] and all(seq1(r.resname)==human[r.id[1]-1] for r in model[0]['R']))
    scores=read_tsv(OUT/'energy_terms.tsv');check('score_count',len(scores)==360,len(scores))
    for r in scores:
        ident=r['model'];mutation=r['mutation'];key=mutation.replace('+','_');folder=OUT/'models'/ident
        p=folder/('wt.pdb' if mutation=='WT' else key+('_fixed' if r['mode']=='fixed_e2' else '')+'.pdb');model=read_structure(p)
        tag=ident+'_wt' if mutation=='WT' else ident+'_'+key+('_fixed' if r['mode']=='fixed_e2' else '')
        terms=energy(p,tag)
        check(f'logged_score:{ident}:{mutation}:{r["mode"]}',all(abs(float(r.get(k,0))-v)<1e-9 for k,v in terms.items()))
        if mutation!='WT':
            check(f'mutation:{ident}:{mutation}:{r["mode"]}',all(seq1(model[0][c][(' ',pos,' ')].resname)==b for c,pos,a,b in mutations(mutation)))
        if r['mode']=='fixed_e2':
            actual=[line for line in p.read_text().splitlines() if line.startswith('ATOM') and line[21]=='E'];expected=[line for line in (folder/'wt.pdb').read_text().splitlines() if line.startswith('ATOM') and line[21]=='E'];check(f'fixed_e2:{ident}:{mutation}',actual==expected)
    panel=read_tsv(OUT/'balanced_development_panel.tsv');state_rows=read_tsv(OUT/'state_manifest.tsv');check('all_balanced_states_built',{r['haplotype_id'] for r in panel}=={r['haplotype_id'] for r in state_rows})
    for r in state_rows:
        p=ROOT/r['path'];s=read_structure(p);check('state_hash:'+r['haplotype_id'],sha256(p)==r['sha256'])
        check('state_mutations:'+r['haplotype_id'],all(seq1(s[0]['E'][(' ',int(text[1:-1]),' ')].resname)==text[-1] for text in json.loads(r['mutations'])))
        check('state_missing_coordinates:'+r['haplotype_id'],json.loads(r['unmodeled_h77_positions'])==[415,416,417])
    dev=json.loads((OUT/'development_sequences.json').read_text());reserve=json.loads((OUT/'quarantine/future_reserve.json').read_text());old=read_tsv(ROOT/'results/hcv_sequence_manifest.tsv');check('no_reserved_archived_accessions',not ({r['accession'] for r in old}&{r['accession'] for r in reserve}))
    check('no_reserved_development_component',not ({r['component'] for r in dev}&{r['component'] for r in reserve}))
    check('no_reserved_development_sequence',not ({r['sequence_hash'] for r in dev}&{r['sequence_hash'] for r in reserve}))
    check('no_reserved_near_development_haplotype',all(sum(a!=b for a,b in zip(x['haplotype'],y['haplotype']))>1 for x in reserve for y in dev))
    coverage=json.loads((OUT/'coverage_summary.json').read_text());indices=[coverage['expanded_interface_positions'].index(p) for p in coverage['original_interface_positions']];oldseq=json.loads((ROOT/'data/processed/discovery_sequences.json').read_text())+json.loads((ROOT/'data/processed/quarantine/heldout_sequences.json').read_text());check('no_reserved_near_exposed_original_haplotype',all(sum(x['haplotype'][i]!=a for i,a in zip(indices,y['haplotype']))>1 for x in reserve for y in oldseq))
    check('genotype_coverage_target',all(float(r['coverage'])>=.9 for r in read_tsv(OUT/'balanced_genotype_coverage.tsv')))
    check('unmodeled_genotype8_explicit','8' in coverage['genotypes_missing'])
    # A fresh executable check anchors cached score parsing to current coordinate bytes.
    for ident,filename,expected_tag in [('archive','wt.pdb','archive_wt'),('archive','T163A_fixed.pdb','archive_T163A_fixed'),('bound_ae','wt.pdb','bound_ae_wt'),('g8q_b','F186L.pdb','g8q_b_F186L')]:
        p=OUT/'models'/ident/filename;expected=energy(p,expected_tag)['Total'];actual,_=invoke('ComputeBinding',p,'verification_'+expected_tag,['--split=E,R'])
        import re
        values=re.findall(r'^\s*Total\s*=\s*([-+\d.eE]+)',actual,re.M)
        check('fresh_score:'+expected_tag,len(values)==1 and abs(float(values[0])-expected)<.01)
    passed=all(c['passed'] for c in checks)
    write_json(OUT/'verification.json',dict(timestamp=now(),passed=passed,checks=len(checks),failed=[c for c in checks if not c['passed']],notes='Data checks; no biological validation claim',project_bytes=footprint(),free_disk_bytes=__import__('shutil').disk_usage(ROOT).free))
    if not passed:raise RuntimeError('Follow-up verification failed')
    print('passed verification checks',len(checks),flush=True)
