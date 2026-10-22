#!/usr/bin/env python3
"""Capture-specific QC diagnostics and sparse scDblFinder exchange on AutoDL."""
import argparse,datetime,gzip,json,os,platform
from pathlib import Path
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io,sparse

def mad(x): return float(np.median(np.abs(x-np.median(x))))
def low_bound(x):
    values=np.log1p(x);return float(max(1,np.expm1(np.median(values)-4*1.4826*mad(values))))
def high_bound(x):
    values=np.log1p(x);return float(np.expm1(np.median(values)+4*1.4826*mad(values)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--capture',action='append');p.add_argument('--study');a=p.parse_args()
    if platform.system()!='Linux' or os.getenv('FIVE_STUDIES_COMPUTE')!='autodl':raise SystemExit('AutoDL only')
    paths=sorted((a.project_root/'inputs/processed').glob('GSE*/*/raw.h5ad'))
    if a.study:paths=[x for x in paths if x.parent.parent.name==a.study]
    if a.capture:paths=[x for x in paths if x.parent.name in a.capture]
    for path in paths:
        destination=a.project_root/'results/scRNA_qc'/path.parent.name
        if (destination/'prepare_receipt.json').exists():continue
        destination.mkdir(parents=True,exist_ok=True);obj=ad.read_h5ad(path);X=obj.X.tocsr();X.eliminate_zeros()
        total=np.asarray(X.sum(axis=1)).ravel();genes=np.diff(X.indptr)
        mt=obj.var.gene_symbol.str.upper().str.startswith('MT-').to_numpy()
        pct=np.asarray(X[:,mt].sum(axis=1)).ravel()/np.maximum(total,1)*100
        low_g,low_c=low_bound(genes),low_bound(total)
        mt_mad=mad(pct);mt_limit=min(100,float(np.median(pct)+4*1.4826*mt_mad)) if mt_mad>0 else float(np.quantile(pct,.99))
        qc=obj.obs.copy();qc['total_counts']=total;qc['n_genes_by_counts']=genes;qc['pct_counts_mt']=pct
        qc['low_genes']=genes<low_g;qc['low_counts']=total<low_c;qc['high_mt']=(pct>mt_limit) if mt.any() else False
        qc['high_counts_flag']=total>high_bound(total);qc['high_genes_flag']=genes>high_bound(genes)
        qc['preliminary_pass']=~(qc.low_genes|qc.low_counts|qc.high_mt)
        qc['exclusion_reason']=['|'.join(k for k in ['low_genes','low_counts','high_mt'] if bool(qc.iloc[i][k])) for i in range(len(qc))]
        # Stress/apoptosis are recorded as descriptive marker detection, never default deletion rules.
        for label,panel in {'stress':['FOS','JUN','HSPA1A','HSPA1B','DDIT3'],'apoptosis':['BAX','BBC3','PMAIP1','CASP3']}.items():
            present=obj.var.gene_symbol.isin(panel).to_numpy();qc[label+'_marker_UMI']=np.asarray(X[:,present].sum(axis=1)).ravel();qc[label+'_assayed_markers']=int(present.sum())
        qc.to_csv(destination/'cell_qc_audit.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
        keep=qc.preliminary_pass.to_numpy();selected=X[keep,:].T.tocoo()
        with gzip.open(destination/'counts_for_doublets.mtx.gz','wb') as stream:io.mmwrite(stream,selected,field='integer')
        qc.loc[keep].to_csv(destination/'cells_for_doublets.tsv',sep='\t',index_label='cell_id')
        obj.var.to_csv(destination/'genes_for_doublets.tsv',sep='\t',index_label='source_row_key')
        fig,axes=plt.subplots(1,3,figsize=(10,3))
        for ax,values,label,threshold in zip(axes,[np.log1p(genes),np.log1p(total),pct],['log1p genes','log1p counts','Mitochondrial %'],[np.log1p(low_g),np.log1p(low_c),mt_limit]):
            ax.hist(values,bins=60,color='#316B9B');ax.axvline(threshold,color='#C24A3A');ax.set_xlabel(label)
        fig.suptitle(path.parent.name);fig.tight_layout();fig.savefig(destination/'qc_distributions.pdf');plt.close(fig)
        receipt={'capture':path.parent.name,'input':str(path),'input_cells':len(qc),'preliminary_pass':int(keep.sum()),
                 'thresholds':{'min_genes':low_g,'min_counts':low_c,'max_mt_pct':mt_limit,'mitochondrial_genes_assayed':int(mt.sum())},
                 'rule':'capture-specific lower log1p genes/counts 4 scaled MAD; upper mt 4 scaled MAD (99% fallback if MAD=0); obvious zero counts excluded; upper genes/counts flagged, scDblFinder handles doublets',
                 'ambient_status':'published EmptyDrops/filtered matrix; no unfiltered empty-droplet matrix available, ambient inference limited and not silently deleted',
                 'stress_apoptosis':'flag only','program_effects_seen_at_threshold_choice':False,
                 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
        (destination/'prepare_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
