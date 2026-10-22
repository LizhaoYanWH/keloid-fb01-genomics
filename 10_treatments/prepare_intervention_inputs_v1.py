"""Input-only admission: public units/columns/axes/coverage, no response statistics."""
from pathlib import Path
import os,json,hashlib,gzip,csv,io,re,collections,shutil,resource
P=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
D=P/'mechanisms_20261021/literature_intervention_v1';S=D/'sources';O=Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
def need(x,m):
 if not x:raise RuntimeError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def soft(p):
 out=[];r=None
 for l in gzip.open(p,'rt'):
  l=l.rstrip('\n')
  if l.startswith('^SAMPLE ='):r={'GSM':l.split(' = ',1)[1]};out.append(r)
  elif r is not None and l.startswith('!Sample_') and ' = ' in l:k,v=l.split(' = ',1);r.setdefault(k,[]).append(v)
 return out
def main():
 os.sched_setaffinity(0,[min(os.sched_getaffinity(0))]);os.environ['OPENBLAS_NUM_THREADS']='1'
 import pandas as pd,numpy as np
 need(not(D/'ADMISSION_COMPLETE.json').exists(),'Input admission already exists')
 need(min(shutil.disk_usage('/root').free,shutil.disk_usage('/root/autodl-tmp').free)>=5*2**30,'Two5GiB')
 HG=O/'shared_inputs/programs_01_03/sources/hgnc_complete_set.txt';need(sha(HG)=='12363916be36f8f6f91879b0faade4eeaac973722de7136ec61329b77f2c2f91','HGNC changed')
 h=pd.read_csv(HG,sep='\t',dtype=str,keep_default_na=False);h=h.loc[h.status.eq('Approved')&h.ensembl_gene_id.ne('')].copy();unique=h.ensembl_gene_id.value_counts();h=h.loc[h.ensembl_gene_id.map(unique).eq(1)]
 h[['ensembl_gene_id','symbol','entrez_id','locus_type']].to_csv(D/'all_unique_approved_HGNC_gene_axis.tsv',sep='\t',index=False)
 ext=P/'extensions_20261021/external_RNA_v1/FREEZE.json';e=json.loads(ext.read_text());need(sha(ext)=='8795fed83ba83fe7141e57e8aa16b6f28390ebaac7a3a412f17ed097e1c38b7a','Original shared25 freeze changed')
 shared=e['fallback_shared_FB01_nonGS_module'];need(len(shared)==25,'Original25 required')
 gene_count=h.loc[h.entrez_id.ne(''),'entrez_id'].value_counts();mapping={r.entrez_id:r.ensembl_gene_id for r in h.itertuples() if r.entrez_id and gene_count[r.entrez_id]==1}
 programs=[{'program':'fixed_shared_nonGS_FB01_25','source_id':g,'ensembl_gene_id':g,'weight':1.,'mapping_status':'FROZEN_ENSEMBL'} for g in shared]
 for g,w in e['original_GS_weights'].items():programs.append({'program':'original_MAGMA1000_Z_RNA','source_id':g,'ensembl_gene_id':mapping.get(g,''),'weight':w,'mapping_status':'UNIQUE_HGNC_ENTREZ_TO_ENSEMBL' if g in mapping else 'UNRESOLVED_NA_NO_BACKFILL'})
 pd.DataFrame(programs).to_csv(D/'fixed2_program_members_and_weights.tsv',sep='\t',index=False)
 can=e['original9_candidates'];pd.DataFrame(can).to_csv(D/'fixed9_candidates.tsv',sep='\t',index=False)
 metadata=[];studies={};axes={};inputs=[HG,ext]
 for g in ['GSE232079','GSE232435']:
  sp=S/(g+'_family.soft.gz');records=soft(sp);cp=S/('GSE232079_counts.txt.gz' if g=='GSE232079' else 'GSE232435_counts_ready.txt.gz');x=pd.read_csv(cp,sep='\t',index_col=0)
  need(x.index.is_unique and np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all(),'Count matrix identity/finite/nonnegative')
  need(np.array_equal(x.to_numpy(),np.rint(x.to_numpy())),'Count matrix not original integer-valued; do not round')
  axes[g]=set(x.index.astype(str));need(len(records)==x.shape[1],'GSM/column cardinality')
  for r in records:
   ch=dict(z.split(': ',1) for z in r['!Sample_characteristics_ch1']);title=r['!Sample_title'][0]
   if g=='GSE232079':
    col=title;strain=ch['cell line'];treatment=ch['treatment'];rep=title.rsplit(' ',1)[1];passage=ch['passage'];kind='primary_keloid' if strain in ['PK1','PK2','PK3','PK4'] else 'adjacent_primary_normal' if strain=='PNF4' else 'commercial_keloid' if strain=='CRL-1762' else 'commercial_normal';patient='PK4' if strain=='PNF4' else strain
   else:
    matches=[z for z in r['!Sample_description'] if z in x.columns];need(len(matches)==1,'TAK1 original count column not uniquely described');col=matches[0];strain=re.search(r'_(N[568])_FS$',title).group(1);treatment={'Control':'Control','TGF-Beta':'TGFbeta','HS-276+TGF-':'HS276_TGFbeta'}[ch['treatment']];rep='1';passage='not_reported';kind='primary_foreskin';patient=strain
   need(col in x.columns,'Count column/GSM exact title-description mismatch')
   metadata.append({'dataset':g,'GSM':r['GSM'],'count_column':col,'source_strain':strain,'biological_patient_or_source':patient,'material':kind,'treatment':treatment,'culture_repeat':rep,'passage':passage,'primary_donor_inference':str(kind in ['primary_keloid','primary_foreskin']),'pairing_status':'PAIRED_WITHIN_SOURCE','dose_time':'corin1uM24h_vs_DMSO' if g=='GSE232079' else 'TGFbeta10ngmL_then_HS27610uM_post24h_GEO;paper_dose_timing_context'})
  studies[g]={'shape':list(x.shape),'count_mode':'provider_raw_integer_valued_counts','min':float(x.to_numpy().min()),'max':float(x.to_numpy().max()),'column_mapping':'GSM exact title' if g=='GSE232079' else 'GSM Sample_description original internal column','full_assayed_gene_axis':len(axes[g]),'noncoding_axis_retained':True};inputs.extend([sp,cp])
 mp=S/'GSE7890_series_matrix.txt.gz';pre=[];body=[];inside=False
 for l in gzip.open(mp,'rt'):
  if l.startswith('!series_matrix_table_begin'):inside=True;continue
  if l.startswith('!series_matrix_table_end'):inside=False;continue
  if inside:body.append(l)
  else:pre.append(l)
 fields=[next(csv.reader([l],delimiter='\t')) for l in pre if l.startswith('!Sample_')];get=lambda k:[a[1:] for a in fields if a[0]==k]
 ids=get('!Sample_geo_accession')[0];titles=get('!Sample_title')[0];character=get('!Sample_characteristics_ch1');matrix=pd.read_csv(io.StringIO(''.join(body)),sep='\t',index_col=0)
 need(list(matrix.columns)==ids and len(ids)==19 and matrix.index.is_unique and np.isfinite(matrix.to_numpy()).all(),'HC matrix19 identity/finite')
 need((matrix.to_numpy()>0).all(),'Positive provider intensity required for predetermined log2')
 for i,gsm in enumerate(ids):
  ch=dict(a[i].split(': ',1) for a in character);disease='keloid' if ch['Disorder']=='keloid' else 'normal_scar';tr='HC' if ch['Condition']=='+hydrocortisone' else 'NoHC';strain=ch['strain']
  metadata.append({'dataset':'GSE7890','GSM':gsm,'count_column':gsm,'source_strain':strain,'biological_patient_or_source':strain,'material':disease,'treatment':tr,'culture_repeat':'three_original_cultures_pooled_before_one_array','passage':'low_passage','primary_donor_inference':'True','pairing_status':'UNPAIRED_NO_HC_RETAINED_NOT_PAIRED_TEST' if strain=='58' else 'PAIRED_WITHIN_STRAIN','dose_time':'HC1.5uM_to_confluency10days'})
 annot=S/'GPL570.annot.gz';start=False;lines=[]
 for l in gzip.open(annot,'rt'):
  if l.startswith('!platform_table_begin'):start=True;continue
  if l.startswith('!platform_table_end'):start=False;continue
  if start:lines.append(l)
 a=pd.read_csv(io.StringIO(''.join(lines)),sep='\t',dtype=str,keep_default_na=False)
 need('Gene ID' in a.columns and 'ID' in a.columns,'Official GPL570 Gene ID/ID columns')
 amap=[]
 for r in a.to_dict('records'):
  ent=r['Gene ID'].strip();ens=mapping.get(ent,'') if re.fullmatch(r'[0-9]+',ent) else ''
  amap.append({'probe_id':r['ID'],'original_Gene_ID':ent,'ensembl_gene_id':ens,'status':'UNIQUE_APPROVED_MAPPING' if ens else 'MULTIPLE_OR_UNRESOLVED_GENE_ID_NA'})
 pd.DataFrame(amap).to_csv(D/'GPL570_fixed_all_probe_to_HGNC_Ensembl.tsv.gz',sep='\t',index=False)
 axes['GSE7890']={r['ensembl_gene_id'] for r in amap if r['ensembl_gene_id'] and r['probe_id'] in matrix.index}
 studies['GSE7890']={'shape':list(matrix.shape),'provider_declares_RMA':all('RMA' in z for z in get('!Sample_data_processing')[0]),'actual_export_positive_intensity_min':float(matrix.to_numpy().min()),'actual_export_positive_intensity_max':float(matrix.to_numpy().max()),'analysis_input_mode':'log2 positive provider normalized expression export, no second normalization/no CEL/RMA refit','full_unique_HGNC_mapped_gene_axis':len(axes['GSE7890']),'probe_rule':'drop ambiguous or unapproved original Gene ID; average ALL uniquely mapped probes in log2 space, not highest effect probe; retain all original probes in source','pairing':'5KL and4normal_scar matched strains; untreated58 retained only descriptive; numbering not pairing'}
 inputs.extend([mp,annot,S/'GSE7890_series_metadata.txt'])
 pd.DataFrame(metadata).to_csv(D/'all56_GSM_sample_metadata.tsv',sep='\t',index=False)
 coverage=[]
 for g,axis in axes.items():
  for pr in programs:coverage.append({'dataset':g,'kind':'program','name':pr['program'],'source_id':pr['source_id'],'ensembl_gene_id':pr['ensembl_gene_id'],'weight':pr['weight'],'assayed':bool(pr['ensembl_gene_id'] in axis)})
  for ca in can:coverage.append({'dataset':g,'kind':'candidate','name':ca['approved_symbol'],'source_id':ca['gene_id'],'ensembl_gene_id':ca['gene_id'],'weight':None,'assayed':ca['gene_id'] in axis})
 pd.DataFrame(coverage).to_csv(D/'all_program_and9candidate_assay_coverage.tsv',sep='\t',index=False)
 semantic={'PMC11995714_corin_BioC.xml':'RETRIEVAL_NOT_ARTICLE: HTTP200 HTML No result; not full text. GEO raw data/protocol independently read; PMC indexed full-body Methods/results/FigS4 reviewed.','PMC11995714_corin_fullText.xml':'EuropePMC HTTP500 preserved, no fullXML available by attempted API.','GSE261116':'Five treatment samples each replicate1, donor pairing not exposed in GEO; not immediate donor-inference analysis; no count matrices downloaded.'}
 rows=pd.DataFrame(metadata);need(len(rows)==56,'28+9+19')
 need(set(rows.loc[(rows.dataset=='GSE232079')&rows.primary_donor_inference.eq('True'),'source_strain'])=={'PK1','PK2','PK3','PK4'},'Exactly4clinicalprimarypatients')
 stats={g:{'programs':{pr:sum(x['dataset']==g and x['kind']=='program' and x['name']==pr and x['assayed'] for x in coverage) for pr in ['fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA']},'assayed_candidates':[x['name'] for x in coverage if x['dataset']==g and x['kind']=='candidate' and x['assayed']]} for g in axes}
 out={'status':'INPUT_ADMISSION_COMPLETE_NO_RESPONSE_STATISTICS','source_sha256':sha(Path(__file__)),'studies':studies,'assay_coverage':stats,'exact_sample_metadata_sha256':sha(D/'all56_GSM_sample_metadata.tsv'),'biological_units':'Corin4patientpairs after withinpatient/condition culture pooling; TAK1 threeprimaryforeskin source-lines; HC5KL+4ordinaryscar strainpairs, each array already triplicateculturepool. No cells/culture repeats counted as patients.','semantic_retrieval_status':semantic,'source_data_bytes':sum(q.stat().st_size for q in S.iterdir() if q.is_file()),'inputs':[rec(p) for p in inputs],'outputs':[rec(p) for p in D.iterdir() if p.is_file() and p.suffix in ['.tsv','.gz']],'no_differential_expression_or_drug_response_computed':True}
 (D/'ADMISSION_COMPLETE.json').write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n');print(json.dumps({'status':out['status'],'studies':studies,'coverage':stats,'metadata_sha256':out['exact_sample_metadata_sha256'],'ADMISSION_sha256':sha(D/'ADMISSION_COMPLETE.json')}))
if __name__=='__main__':main()
