"""Reproducible experimental glycan proximity and literature interface annotations."""
import json
import numpy as np
from Bio.PDB import MMCIFParser
from scipy.spatial.distance import cdist
from .common import ROOT,fetch,read_tsv,write_tsv
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from .structures import column

def annotations():
    metadata=json.loads((ROOT/'data/processed/humanization.json').read_text())
    source=MMCIFParser(QUIET=True,auth_chains=False).get_structure('7MWX',str(ROOT/'data/raw/7MWX.cif'))
    glycan_names={'NAG','BMA','MAN','FUC','GAL','SIA','NDG'}
    glycans=[a.coord for chain in source[0] for r in chain if r.resname in glycan_names for a in r.get_atoms() if a.element not in {'H','D'}]
    if not glycans:
        raise RuntimeError('No observed experimental glycans; proximity table not applicable')
    rows=[]
    for residue in source[0][metadata['e2_chain']]:
        if residue.id[0]!=' ':
            continue
        xyz=[a.coord for a in residue.get_atoms() if a.element not in {'H','D'}]
        distance=float(cdist(xyz,glycans).min())
        rows.append(dict(e2_structure_position=residue.id[1],minimum_experimental_glycan_distance_angstrom=distance,
                         glycan_proximity_flag=distance<=5,
                         interpretation='Proximity to any observed experimental glycan; does not establish steric shielding on a complete virion'))
    write_tsv(ROOT/'results/e2_glycan_proximity.tsv',rows)
    article=fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC1563869/',ROOT/'data/raw/PMC1563869.html','PMC1563869')
    if '10.1128/JVI.00271-06'.lower() not in article.read_text(encoding='utf-8').lower():
        raise RuntimeError('E2 literature interface source DOI not verified')
    interface={int(r['h77_position']) for r in read_tsv(ROOT/'results/e2_interface_residues.tsv') if float(r['cutoff_angstrom'])==5}
    rows=[dict(h77_position=p,wild_type=aa,experimental_annotation='Published critical CD81-binding position; not reclassified by geometric cutoff',
               inside_primary_humanized_interface=p in interface,source='Owsianka et al. J Virol 2006',doi='10.1128/JVI.00271-06',
               primary_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC1563869/') for p,aa in [(420,'W'),(527,'Y'),(529,'W'),(530,'G'),(535,'D')]]
    write_tsv(ROOT/'results/e2_literature_interface_annotations.tsv',rows)
    human=MMCIFParser(QUIET=True,auth_chains=False).get_structure('3X0E',str(ROOT/'data/raw/3X0E.cif'))
    cif=MMCIF2Dict(str(ROOT/'data/raw/3X0E.cif')); bonds=[]
    keys=['conn_type_id','ptnr1_label_asym_id','ptnr2_label_asym_id','ptnr1_auth_seq_id','ptnr2_auth_seq_id']
    for kind,a,b,p,q in zip(*(column(cif,'_struct_conn.'+key) for key in keys)):
        if kind=='disulf' and a==b==metadata['human_chain']:
            p,q=int(p),int(q)
            distance=float(np.linalg.norm(human[0][a][p]['SG'].coord-human[0][b][q]['SG'].coord))
            bonds.append(dict(source_pdb='3X0E',label_chain=a,position_1=p,position_2=q,experimental_sg_distance_angstrom=distance,
                              evidence='Experimental mmCIF _struct_conn disulf record; native cysteine pair preserved'))
    if not bonds:
        raise RuntimeError('Selected human scaffold has no verified disulfide records')
    write_tsv(ROOT/'results/disulfide_provenance.tsv',bonds)
