#!/usr/bin/env python
"""One bounded descriptive bridge; immutable old inputs, new system-disk output."""
import os
for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[k]="1"
import argparse, csv, gzip, hashlib, itertools, json, math, resource, shutil, sys, time
from pathlib import Path
import h5py, numpy as np, pandas as pd
from scipy.sparse import csr_matrix
P=Path("/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state")
O=Path("/root/autodl-tmp/Five_Studies_Execution_20261020")
HERE=P/"extensions_20261021/gene_state_bridge_v1"
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()
def jwrite(p,x):
    q=Path(str(p)+".tmp")
    q.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    q.replace(p)
def clean(v):
    if isinstance(v,dict): return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)): return [clean(x) for x in v]
    if isinstance(v,np.generic): v=v.item()
    if isinstance(v,float) and not math.isfinite(v): return None
    return v
def harray(node):
    if isinstance(node,h5py.Dataset):
        z=node[:]
        if z.dtype.kind in ("O","S","U"): return np.asarray([x.decode() if isinstance(x,bytes) else str(x) for x in z])
        return z
    if "codes" in node and "categories" in node:
        cats=harray(node["categories"]); codes=node["codes"][:]
        return np.asarray([cats[c] if c>=0 else "" for c in codes])
    raise ValueError("Unsupported h5ad array")
def mem_available():
    vals=[]
    kv=dict(line.split(":",1) for line in Path("/proc/meminfo").read_text().splitlines())
    vals.append(int(kv["MemAvailable"].split()[0])*1024)
    for base,lim,cur in [("/sys/fs/cgroup","memory.max","memory.current"),
                         ("/sys/fs/cgroup/memory","memory.limit_in_bytes","memory.usage_in_bytes")]:
        try:
            cap=(Path(base)/lim).read_text().strip()
            if cap!="max": vals.append(int(cap)-int((Path(base)/cur).read_text()))
        except FileNotFoundError: pass
    return min(vals)
def guard(startup=False):
    # 4GiB bounded RSS job + 8GiB floor + root-reserved 12GiB peer.
    required=(24 if startup else 20)*1024**3
    if mem_available()<required: raise RuntimeError("RESOURCE_HOLD: available RAM below job/floor/peer reserve")
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024>4*1024**3:
        raise RuntimeError("RESOURCE_HOLD: observed RSS exceeded4GiB")
    for p,reserve in [(Path("/root"),64*1024**2),(Path("/root/autodl-tmp"),0)]:
        if shutil.disk_usage(p).free<5*1024**3+reserve: raise RuntimeError("DISK_HOLD: actual5GiB disk floor")
def verify(p,s):
    p=Path(p)
    got=sha(p)
    if got!=s: raise RuntimeError("INPUT_SHA_MISMATCH "+str(p))
    return {"path":str(p),"bytes":p.stat().st_size,"sha256":got}
def freeze_inputs():
    frozen=HERE/"frozen_inputs_v1.json"
    f=json.loads(frozen.read_text())
    assert f["scope"]=="FINITE_MAIN_PP4_GENE_TO_FB01_DESCRIPTIVE_BRIDGE"
    rec=[]
    for x in f["files"]: rec.append(verify(x["path"],x["sha256"]))
    old=json.loads((P/"reviews/scDRS_production_preparation_v1/RNA_freeze_v1/COMPLETE.json").read_text())
    return f,old,rec
def candidate_tables(f):
    t=pd.read_csv(f["coloc_profile_path"],sep="\t")
    assert len(t)==2391
    t["primary_PP4"]=pd.to_numeric(t["primary_PP4"],errors="coerce")
    t["low_p12_PP4"]=pd.to_numeric(t["low_p12_PP4"],errors="coerce")
    a=t[t.primary_PP4.ge(.8)].copy()
    assert len(a)==11 and a.gene_id.nunique()==9
    a["low_prior_robust"]=a.low_p12_PP4.ge(.8)
    h=pd.read_csv(f["HGNC_path"],sep="\t",dtype=str,keep_default_na=False)
    h=h[h.status.eq("Approved")]
    records=[]
    for gene in sorted(a.gene_id.unique()):
        q=h[h.ensembl_gene_id.eq(gene)]
        records.append({"gene_id":gene,"HGNC_rows":len(q),
            "approved_symbol":q.iloc[0]["symbol"] if len(q)==1 else "",
            "locus_group":q.iloc[0]["locus_group"] if len(q)==1 else "",
            "locus_type":q.iloc[0]["locus_type"] if len(q)==1 else "",
            "HGNC_mapping_status":"ONE_APPROVED" if len(q)==1 else "AMBIGUOUS_OR_MISSING_NA",
            "primary_supported_units":int((a.gene_id==gene).sum()),
            "low_prior_robust_units":int(((a.gene_id==gene)&a.low_prior_robust).sum())})
    return t,a,pd.DataFrame(records)
def risk_direction(a,f,old_records):
    keys=set(map(tuple,a[["context","region_id","gene_id"]].itertuples(index=False,name=None)))
    leads={}; total=0
    # Stream entire original completed SNP table; no new ABF fit, no SNP cutoff.
    with gzip.open(f["primary_snp_path"],"rt") as z:
        rd=csv.DictReader(z,delimiter="\t")
        for row in rd:
            total+=1; key=(row["context"],row["region_id"],row["gene_id"])
            if total%1000000==0: guard()
            if key not in keys: continue
            pp=float(row["SNP_PP_H4_conditional_on_H4"])
            cur=leads.get(key)
            if cur is None or pp>float(cur["SNP_PP_H4_conditional_on_H4"]) or (pp==float(cur["SNP_PP_H4_conditional_on_H4"]) and row["variant_id_GRCh38"]<cur["variant_id_GRCh38"]):
                leads[key]=row
    manifest=pd.read_csv(f["harmonized_manifest_path"],sep="\t",dtype=str,keep_default_na=False)
    selected={}
    for key in keys:
        q=manifest[manifest.context.eq(key[0])&manifest.region_id.eq(key[1])&manifest.gene_id.eq(key[2])]
        assert len(q)==1
        path=q.iloc[0]["path"]; digest=q.iloc[0]["input_SHA256"]
        selected.setdefault(path,{"sha256":digest,"keys":set()})["keys"].add(key)
    matched={}
    for path,z in selected.items():
        old_records.append(verify(path,z["sha256"]))
        with gzip.open(path,"rt") as fp:
            for row in csv.DictReader(fp,delimiter="\t"):
                for key in z["keys"]:
                    if row["gene_id"]==key[2] and row["variant_id_GRCh38"]==leads[key]["variant_id_GRCh38"]:
                        assert key not in matched
                        matched[key]=row
    out=[]
    for key in sorted(keys):
        base=a[a.context.eq(key[0])&a.region_id.eq(key[1])&a.gene_id.eq(key[2])].iloc[0].to_dict()
        r=matched.get(key); lead=leads.get(key)
        base.update({"representative_SNP_selection":"maximum original primary SNP posterior conditional on H4; lexicographic tie",
                     "SNP_direction_is_not_MR":"single representative SNP orientation only; not effect of intervention",
                     "risk_allele":"","expression_direction_for_risk_allele":"NA","direction_status":"NA_NO_HARMONIZED_LEAD"})
        if r:
            base.update({"variant_id_GRCh38":r["variant_id_GRCh38"],
                         "SNP_PP_H4_conditional_on_H4":float(lead["SNP_PP_H4_conditional_on_H4"]),
                         "GWAS_beta_ALT_if_direction_resolved":r["GWAS_beta_ALT_if_direction_resolved"],
                         "QTL_beta_ALT":r["QTL_beta_ALT"],"direction_status":r["direction_status"],
                         "GWAS_source_row":r["GWAS_source_row"],"nominal_source_row":r["nominal_source_row"]})
            try: bg=float(r["GWAS_beta_ALT_if_direction_resolved"]); bq=float(r["QTL_beta_ALT"])
            except ValueError: bg=bq=float("nan")
            if math.isfinite(bg) and math.isfinite(bq) and bg!=0 and bq!=0:
                base["risk_allele"]=r["QTL_alt"] if bg>0 else r["QTL_ref"]
                base["expression_direction_for_risk_allele"]="increased" if bg*bq>0 else "decreased"
        out.append(base)
    return pd.DataFrame(out),total
def aggregate(study,old,genes):
    x=old["cohorts"][study]; inp=x["input"]
    meta_path=P/"reviews/scDRS_production_preparation_v1/RNA_freeze_v1"/study/"all_final_cell_metadata.tsv.gz"
    meta=pd.read_csv(meta_path,sep="\t",dtype=str,keep_default_na=False).set_index("cell_id")
    assert meta.index.is_unique
    genes=list(genes); sums={}; qc=[]; axisrows=[]
    documented=set(x["independent_patients"])
    with h5py.File(inp["counts_path"],"r") as h:
        obs=harray(h["obs"][h["obs"].attrs["_index"]])
        var=harray(h["var"][h["var"].attrs["_index"]])
        assert len(obs)==len(meta)==x["n_cells"] and pd.Index(obs).is_unique
        assert set(obs)==set(meta.index)
        meta=meta.loc[obs]
        unique,counts=np.unique(var,return_counts=True)
        assert max(counts)==1, "Duplicated original Ensembl axis: do not silently aggregate"
        pos={g:i for i,g in enumerate(var)}
        ix=[pos[g] for g in genes if g in pos]
        present=[g for g in genes if g in pos]
        syms=harray(h["var/gene_symbol"]) if "var/gene_symbol" in h else np.repeat("",len(var))
        for g in genes:
            axisrows.append({"study":study,"gene_id":g,"raw_axis_status":"ASSAYED" if g in pos else "NOT_ASSAYED_NA",
                "raw_axis_index":pos.get(g,""),"source_gene_symbol":syms[pos[g]] if g in pos else ""})
        countsnode=h["layers/counts"]
        assert countsnode.attrs["encoding-type"]=="csr_matrix"
        shape=tuple(countsnode.attrs["shape"])
        assert shape==(len(obs),len(var))
        indptr=countsnode["indptr"]
        rowsize=512
        for start in range(0,len(obs),rowsize):
            if start%(rowsize*32)==0: guard()
            end=min(start+rowsize,len(obs))
            ptr=indptr[start:end+1]; lo,hi=int(ptr[0]),int(ptr[-1])
            data=countsnode["data"][lo:hi]; ids=countsnode["indices"][lo:hi]
            assert np.isfinite(data).all() and (data>=0).all() and np.equal(data,np.rint(data)).all()
            mat=csr_matrix((data,ids,ptr-lo),shape=(end-start,len(var)))
            lib=np.asarray(mat.sum(axis=1)).ravel().astype(np.int64)
            prior=pd.to_numeric(meta.iloc[start:end]["full_input_UMI"]).to_numpy(dtype=np.int64)
            assert np.array_equal(lib,prior),"full library differs from frozen raw counts"
            small=mat[:,ix].toarray().astype(np.int64)
            for k in range(end-start):
                m=meta.iloc[start+k]; patient=m["independent_patient_id"]
                fb=m["lineage"]=="fibroblast"; state=m["FB_state"]
                if fb: assert state in x["FB_state_family"]
                # Retain all lineages, states, captures; documented patients are separate from unknown contexts.
                is_documented=patient in documented
                unit=patient if is_documented else "UNCODED_CONTEXT::"+m["capture_id"]
                patient_status="DOCUMENTED_PATIENT" if is_documented else "UNCODED_CONTEXT_NOT_PATIENT"
                groups=[("lineage",m["lineage"]),("all_cells","ALL")]
                if fb: groups += [("FB_state",state),("FB_contrast","FB01" if state=="FB01" else "otherFB"),("all_FB","ALL_FB")]
                for level,label in groups:
                    for resolution,unitkey in [("patient",unit),("capture",m["capture_id"])]:
                        key=(resolution,unitkey,patient,patient_status,m["capture_id"] if resolution=="capture" else "",level,label)
                        z=sums.setdefault(key,{"n":0,"lib":0,"raw":np.zeros(len(present),dtype=np.int64),
                            "detect":np.zeros(len(present),dtype=np.int64),"cell_cp10k":np.zeros(len(present),dtype=np.float64)})
                        z["n"]+=1; z["lib"]+=int(lib[k]); z["raw"]+=small[k];z["detect"]+=(small[k]>0)
                        if lib[k]>0: z["cell_cp10k"]+=small[k]/float(lib[k])*1e4
            del mat,data,ids,small
    result=[]
    for key,z in sorted(sums.items()):
        res,unit,patient,status,capture,level,label=key
        for g in genes:
            presentgene=g in present; j=present.index(g) if presentgene else -1
            result.append({"study":study,"resolution":res,"unit_id":unit,"independent_patient_id":patient,
                "unit_status":status,"capture_id":capture,"group_level":level,"group_label":label,
                "gene_id":g,"assay_status":"ASSAYED" if presentgene else "NOT_ASSAYED_NA","n_cells":z["n"],
                "full_raw_library_UMI":z["lib"],"gene_raw_UMI":int(z["raw"][j]) if presentgene else np.nan,
                "detected_cells":int(z["detect"][j]) if presentgene else np.nan,
                "detection_fraction":float(z["detect"][j]/z["n"]) if presentgene else np.nan,
                "pseudobulk_CP10K":float(z["raw"][j]/z["lib"]*1e4) if presentgene and z["lib"] else np.nan,
                "mean_cell_CP10K":float(z["cell_cp10k"][j]/z["n"]) if presentgene else np.nan})
    table=pd.DataFrame(result)
    patienttable=table[table.resolution.eq("patient")&table.group_level.eq("FB_contrast")&table.unit_status.eq("DOCUMENTED_PATIENT")]
    contrasts=[]
    for patient in x["independent_patients"]:
        for gene in genes:
            a=patienttable[patienttable.independent_patient_id.eq(patient)&patienttable.gene_id.eq(gene)&patienttable.group_label.eq("FB01")]
            b=patienttable[patienttable.independent_patient_id.eq(patient)&patienttable.gene_id.eq(gene)&patienttable.group_label.eq("otherFB")]
            assert len(a)==1 and len(b)==1
            a=a.iloc[0];b=b.iloc[0]; assayed=a.assay_status=="ASSAYED"
            contrasts.append({"study":study,"independent_patient_id":patient,"gene_id":gene,
              "assay_status":a.assay_status,"FB01_n_cells":int(a.n_cells),"otherFB_n_cells":int(b.n_cells),
              "FB01_pseudobulk_CP10K":a.pseudobulk_CP10K,"otherFB_pseudobulk_CP10K":b.pseudobulk_CP10K,
              "pseudobulk_CP10K_difference_FB01_minus_otherFB":a.pseudobulk_CP10K-b.pseudobulk_CP10K if assayed else np.nan,
              "FB01_mean_cell_CP10K":a.mean_cell_CP10K,"otherFB_mean_cell_CP10K":b.mean_cell_CP10K,
              "mean_cell_CP10K_difference_FB01_minus_otherFB":a.mean_cell_CP10K-b.mean_cell_CP10K if assayed else np.nan,
              "FB01_detected_fraction":a.detection_fraction,"otherFB_detected_fraction":b.detection_fraction,
              "detection_fraction_difference_FB01_minus_otherFB":a.detection_fraction-b.detection_fraction if assayed else np.nan,
              "contrast_status":"DESCRIPTIVE_NO_PATIENT_P" if assayed else "NOT_ASSAYED_NA"})
    return table,pd.DataFrame(contrasts),pd.DataFrame(axisrows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--run-root-approved",action="store_true");args=ap.parse_args()
    if not args.run_root_approved: raise SystemExit("Prepared source only; root-approved explicit flag required.")
    guard(startup=True); frozenpath=HERE/"frozen_inputs_v1.json"; approval=json.loads((HERE/"ROOT_APPROVAL_TO_RUN.json").read_text())
    assert approval["source_sha256"]==sha(__file__) and approval["freeze_sha256"]==sha(frozenpath)
    assert approval["approved"] is True
    out=HERE/"results"
    if out.exists(): raise RuntimeError("Output exists; no overwrite or repeat run")
    claim=HERE/"RUN_ONCE_CLAIM.json"
    with claim.open("x") as z: json.dump({"PID":os.getpid(),"start_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())},z)
    out.mkdir(); start=time.time()
    f,old,records=freeze_inputs()
    allunits,units,candidates=candidate_tables(f)
    units.to_csv(out/"all_primary_PP4_ge08_units.tsv",sep="\t",index=False)
    candidates.to_csv(out/"all_primary_PP4_ge08_unique_genes.tsv",sep="\t",index=False)
    directions,nprimary=risk_direction(units,f,records)
    directions.to_csv(out/"all_candidate_unit_representative_SNP_direction.tsv",sep="\t",index=False,na_rep="NA")
    allagg=[];allcontrast=[];allaxis=[]
    for study in ("GSE220300","GSE335482"):
        guard();a,b,c=aggregate(study,old,candidates.gene_id.tolist())
        allagg.append(a);allcontrast.append(b);allaxis.append(c)
    agg=pd.concat(allagg,ignore_index=True); con=pd.concat(allcontrast,ignore_index=True); axis=pd.concat(allaxis,ignore_index=True)
    agg.to_csv(out/"all_candidate_expression_complete_groups.tsv.gz",sep="\t",index=False,na_rep="NA",compression="gzip")
    con.to_csv(out/"all_candidates_all8_patient_FB01_contrasts.tsv",sep="\t",index=False,na_rep="NA")
    axis.to_csv(out/"all_candidate_original_count_axis_coverage.tsv",sep="\t",index=False)
    summary=[]
    for (study,gene),g in con.groupby(["study","gene_id"],sort=True):
        vals=g.pseudobulk_CP10K_difference_FB01_minus_otherFB
        summary.append({"study":study,"gene_id":gene,"documented_patients":len(g),
          "assayed_patients":int(vals.notna().sum()),"positive_patients":int((vals>0).sum()),
          "zero_patients":int((vals==0).sum()),"negative_patients":int((vals<0).sum()),
          "median_CP10K_difference":float(vals.median()) if vals.notna().any() else np.nan,
          "min_CP10K_difference":float(vals.min()) if vals.notna().any() else np.nan,
          "max_CP10K_difference":float(vals.max()) if vals.notna().any() else np.nan})
    pd.DataFrame(summary).merge(candidates,on="gene_id",how="left").to_csv(out/"all_candidate_descriptive_direction_summary.tsv",sep="\t",index=False,na_rep="NA")
    jwrite(out/"input_reverification.json",records)
    totalsize=sum(p.stat().st_size for p in out.iterdir())
    assert totalsize<64*1024**2
    files={p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in out.iterdir() if p.is_file()}
    complete={"status":"COMPLETE_FINITE_DESCRIPTIVE_GENE_STATE_BRIDGE_NO_NEW_GENE_SCORE_OR_CAUSAL_FIT",
      "PID":os.getpid(),"source_sha256":sha(__file__),"freeze_sha256":sha(frozenpath),
      "primary_units":len(units),"unique_primary_genes":len(candidates),
      "low_prior_robust_units":int(units.low_prior_robust.sum()),
      "all_primary_SNP_rows_read":nprimary,"RNA_cohorts":2,"documented_patients":8,
      "all_patient_candidate_contrasts":len(con),"complete_group_expression_rows":len(agg),
      "all_original_axes_preserved":True,"noncoding_not_Entrez_filtered":True,
      "patient_p_values_fitted":False,"old_scDRS_or_coloc_recomputed":False,
      "direction_interpretation":"harmonized representative SNP only; not gene causal effect or intervention direction",
      "elapsed_seconds":time.time()-start,"peak_RSS_bytes":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
      "files":files}
    jwrite(out/"COMPLETE.json",clean(complete));print(json.dumps(clean(complete),ensure_ascii=False))
if __name__=="__main__": main()

