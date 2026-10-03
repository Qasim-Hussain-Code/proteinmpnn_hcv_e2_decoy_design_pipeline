"""Curated factual extraction with source-local evidence anchors, not plot digits."""
import re
from .common import ROOT, fetch, now, read_tsv, write_tsv

def ground_truth():
    sources=[('higginbottom2000','Higginbottom et al. J Virol 2000','10.1128/JVI.74.8.3642-3649.2000','PMC111874'),
             ('drummer2002','Drummer, Wilson and Poumbourios. J Virol 2002','10.1128/JVI.76.21.11143-11147.2002','PMC136624'),
             ('bertaux2006','Bertaux and Dragic. J Virol 2006','10.1128/JVI.80.10.4940-4948.2006','PMC1472091'),
             ('flint2006','Flint et al. J Virol 2006','10.1128/JVI.00104-06','PMC1642177')]
    evidence={}
    for sid,citation,doi,pmc in sources:
        path=fetch(f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',ROOT/f'data/raw/{pmc}.html',pmc)
        text=re.sub('<[^>]+>',' ',path.read_text(encoding='utf-8'))
        if 'CD81' not in text or doi.lower() not in text.lower():
            raise RuntimeError(f'Primary article retrieval not verified: {pmc}')
        evidence[sid]=(citation,doi,pmc,text)
    rows=[]
    def add(sid,mutation,direction,construct,assay,strain,measurement,replicates,notes='',include=True):
        citation,doi,pmc,text=evidence[sid]
        pos=int(mutation[1:-1])
        mapping=read_tsv(ROOT/'results/human_residue_mapping.tsv')
        m=next((r for r in mapping if int(r['full_length_position'])==pos),None)
        if m and m['wild_type']!=mutation[0]:
            raise RuntimeError('Ground truth/reference numbering mismatch')
        rows.append(dict(source_id=sid,citation=citation,doi_or_pmid=doi,cd81_construct=construct,assay_type=assay,
                         hcv_e2_or_virus_strain=strain,mutation=mutation,wild_type_residue=mutation[0],position_full_length_cd81=pos,
                         position_structure=m['structure_original_position'] if m else '',measurement=measurement,
                         measurement_units='qualitative unless explicitly stated',reported_direction=direction,normalized_direction=direction,
                         replicate_information=replicates,notes=notes, included=str(include and m is not None).lower(),
                         exclusion_reason='' if include and m else 'Not independently classifiable or coordinate unavailable',
                         primary_source=f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',evidence_anchor=notes))
    for mut,direction,measurement in [('T163A','enhanced','greater E2 signal'),('F186L','reduced','no E2 signal'),
                                     ('E188K','approximately_neutral','similar E2 signal'),('D196E','reduced','intermediate E2 signal')]:
        add('higginbottom2000',mut,direction,'GST human CD81 LEL','soluble E2 binding EIA','E2_661; strain reported in source methods',measurement,
            'Fig 4/5: see source; no numeric replicate count inferred','Figures 4-5; thiocyanate index is not KD')
        add('higginbottom2000',mut,direction,'full-length human CD81 in rat KM3','cell-surface soluble E2 binding','E2_661',measurement,
            'Figures 6-7 representative of two independent experiments','Figure 6; expression and trafficking differ')
    for mut in ['L162P','I182F','N184Y','F186S']:
        add('drummer2002',mut,'reduced','MBP-LEL 113-201 dimer','soluble E2 binding EIA','E2_661myc',
            'decreased signal' if mut=='N184Y' else 'lost signal','Figure 2B; no replicate count inferred','Figure 2B; no plot digitization')
        add('drummer2002',mut,'reduced','full-length CD81 with C8 tag in CHO-K1','cell-surface soluble E2 binding','E2_661myc','lost signal',
            'two independent transfections','Figure 2D; L162P surface detection confounded' if mut=='L162P' else 'Figure 2D',include=mut!='L162P')
    for mut in ['K171A','I181A','I182A','F186A']:
        for assay in ['cell-surface soluble E2 binding','HCVpp entry']:
            add('bertaux2006',mut,'reduced','full-length human CD81 in HepG2',assay,'source E1E2/soluble E2 construct; see methods',
                'significantly decreased','Figure 3; triplicate independent experiments','Figure 3A-B; biochemical and entry outcomes kept separate')
    for mut in ['N184A','E188A','D196A']:
        for assay in ['cell-surface soluble E2 binding','HCVpp entry']:
            add('bertaux2006',mut,'approximately_neutral','full-length human CD81 in HepG2',assay,'source E1E2/soluble E2 construct; see methods',
                'no major receptor phenotype','Figure 3; see source','Results and discussion; neutral is qualitative')
    for mut in ['I182F','N184Y','F186S']:
        add('flint2006',mut,'reduced','full-length human CD81 in HepG2','HCVpp entry','H77/H','moderate reduction; residual entry',
            'Figure 6B; see source','Entry reduction does not imply abolition; no numeric values digitized')
    for mut in ['F186L','E188K','D196E']:
        add('bertaux2006',mut,'not_classifiable','prior studies summarized in discussion','HCVpp entry literature summary',
            'not assigned from review passage','wild-type entry reported for prior work','secondary citation within primary article',
            'Excluded: verify cited original entry study before inclusion',include=False)
    write_tsv(ROOT/'config/cd81_mutation_ground_truth.tsv',rows)
    write_tsv(ROOT/'results/literature_audit.tsv',[dict(source_id=sid,citation=cit,doi=doi,pmcid=pmc,date_checked=now(),
              correction='Prompt misattributed this title/DOI to Flint; actual authors Bertaux and Dragic' if sid=='bertaux2006' else '',
              extraction='Manual factual curation checked against primary full text; qualitative; no plot digitization') for sid,cit,doi,pmc in sources])
