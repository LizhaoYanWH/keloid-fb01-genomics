"""Render five integrated figures from complete, previously computed source tables.
No model fitting, significance testing, signature selection or scientific reruns.
All scientific rendering is confined to the existing AutoDL host.
"""
from pathlib import Path
import os, json, hashlib, shutil, platform, textwrap

ROOT = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
OUT = ROOT/'bmc_genomics_integration_20261021/figures_v4'
INPUTS = {
 'groups':'results/scDRS_fixed47_family_patient_descriptive_v1/all47_official_MC_one_BH_family_with_existing_annotation.tsv',
 'patients':'results/scDRS_fixed47_family_patient_descriptive_v1/all15_FB_states_known_patient_consistency_descriptive.tsv',
 'external':'extensions_20261021/external_RNA_v1/score_v1/all6_patient_FB01_context_score_descriptions.tsv',
 'ablation':'extensions_20261021/regional_exclusion_v1/results_v1/all45_B30_C15_sensitivity_tests_BH45.tsv',
 'time':'normal_repair_timecourse_20261021/results/all24_profile_program_scores.tsv',
 'time_summary':'normal_repair_timecourse_20261021/results/all8_equal_donor_descriptive_summaries.tsv',
 'genes':'extensions_20261021/gene_state_bridge_v1/results/all_primary_PP4_ge08_unique_genes.tsv',
 'coloc':'extensions_20261021/gene_state_bridge_v1/results/all_primary_PP4_ge08_units.tsv',
 'expression':'mechanisms_20261021/regulatory_axis_v1/results/all9_candidate_patient_expression_deltas.tsv',
 'TF':'mechanisms_20261021/regulatory_axis_v1/ECM_shadow_v1/results/fixed15_all3cohorts7contrasts_new_ECM_exclusion_view.tsv',
 'drug_summary':'mechanisms_20261021/literature_intervention_v1/results_v1/allSummary_complete.tsv.gz',
 'drug_delta':'mechanisms_20261021/literature_intervention_v1/results_v1/allDelta_complete.tsv.gz',
 'drug_camera':'mechanisms_20261021/literature_intervention_v1/results_v1/fixed4programs_all7contrasts_CAMERA_BH28.tsv.gz',
 'drug_fry':'mechanisms_20261021/literature_intervention_v1/results_v1/fixed4programs_all7contrasts_FRY_BH28.tsv.gz',
 'state_controls':'mechanisms_20261021/regulatory_axis_v1/state_marker_comparator_v1/results/all16_program_all7_contrast_descriptive_intervals.tsv.gz',
}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rec(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))

def main():
    assert platform.system()=='Linux' and platform.node().startswith('autodl-container-gh1ez37ml9')
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
        os.environ[key]='1'
    os.sched_setaffinity(0,[min(os.sched_getaffinity(0))])
    assert min(shutil.disk_usage('/root').free,shutil.disk_usage('/root/autodl-tmp').free)>5*2**30+64*2**20
    assert not (OUT/'FIGURES_RENDERED_COMPLETE.json').exists(), 'Existing figure outputs; no repeated render'
    OUT.mkdir(parents=True,exist_ok=False)
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm, ListedColormap
    from matplotlib.patches import Rectangle
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.labelsize':7,
        'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,'pdf.fonttype':42,'ps.fonttype':42,
        'image.cmap':'RdBu_r','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,
        'lines.linewidth':1.1,'savefig.facecolor':'white','figure.facecolor':'white'})
    data={k:pd.read_csv(ROOT/v,sep='\t') for k,v in INPUTS.items()}
    sources={k:rec(ROOT/v) for k,v in INPUTS.items()}
    assert len(data['groups'])==47 and len(data['patients'])==15 and len(data['ablation'])==45
    assert len(data['time'])==24 and len(data['genes'])==9 and len(data['state_controls'])==112
    blue='#0072B2';orange='#D55E00';green='#009E73';gray='#666666'
    cohort_colors={'GSE220300':blue,'GSE335482':orange,'GSE163973':green,'GSE241132_D0':gray}
    members=['TAB2','NEDD4','ASIP','LSP1','TNNT3','SETBP1','CARS1-AS1','SETBP1-DT','PRR33']
    names={r.gene_id:r.approved_symbol for r in data['genes'].itertuples()}
    caption_records=[]
    def panel(ax,letter,title):
        ax.set_title(title,loc='left',pad=8)
        ax.text(-.13,1.055,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
    def save(fig,stem,caption,keys):
        assert len(caption.split())<=300
        for ext in ['pdf','svg','png']:
            fig.savefig(OUT/(stem+'.'+ext),dpi=600 if ext=='png' else 300)
        fig.savefig(OUT/(stem+'.preview.png'),dpi=170)
        plt.close(fig)
        (OUT/(stem+'.caption.md')).write_text(caption+'\n')
        (OUT/(stem+'.alt.txt')).write_text(caption+'\n')
        manifest=dict(stem=stem,width_mm=170,source_script_sha256=sha(Path(__file__)),sources={k:sources[k] for k in keys},
            files=[rec(OUT/(stem+'.'+ext)) for ext in ['pdf','svg','png']],caption_words=len(caption.split()),actual_visual_QA='PENDING')
        (OUT/(stem+'.manifest.json')).write_text(json.dumps(manifest,indent=2)+'\n')
        caption_records.append(manifest)
    # Figure1: the complete primary family, without replacing it by selected FB01 points.
    fig=plt.figure(figsize=(170/25.4,180/25.4));gs=fig.add_gridspec(2,3,height_ratios=[.38,1],wspace=.85,hspace=.5)
    ax=fig.add_subplot(gs[0,:]);ax.axis('off')
    steps=[('Published EUR GWAS','Fixed top1,000 MAGMA genes'),('Public single-cell counts','47 complete lineage/state tests'),('Interpretation checks','Patients, ablations, external RNA')]
    for i,(a,b) in enumerate(steps):
        x=.01+i*.34
        ax.add_patch(Rectangle((x,.35),.3,.5,transform=ax.transAxes,facecolor=['#E8F1F7','#FFF0E8','#E7F4EE'][i],edgecolor=gray,lw=.6))
        ax.text(x+.15,.68,a,ha='center',va='center',fontsize=8,transform=ax.transAxes)
        ax.text(x+.15,.48,'\n'.join(textwrap.wrap(b,26)),ha='center',va='center',fontsize=7,transform=ax.transAxes)
        if i<2:ax.annotate('',xy=(x+.337,.6),xytext=(x+.303,.6),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=.8))
    ax.text(0,.04,'One genetic input; separate RNA sources. No individual genotypes or clinical risk prediction.',fontsize=7,transform=ax.transAxes)
    for i,c in enumerate(['GSE241132_D0','GSE220300','GSE335482']):
        ax=fig.add_subplot(gs[1,i]);d=data['groups'].loc[data['groups'].cohort.eq(c)].copy()
        labels=[str(r.group).replace('_',' ')+(' [state]' if r.group_family=='FB_state' else '') for r in d.itertuples()]
        yy=np.arange(len(d));xx=-np.log10(d.FDR_BH_all47_fixed_MC.to_numpy(float))
        for j,r in enumerate(d.itertuples()):
            mark='s' if r.group_family=='FB_state' else 'o';ax.scatter(xx[j],j,edgecolors=cohort_colors[c],marker=mark,s=24,facecolors=cohort_colors[c] if r.FDR_BH_all47_fixed_MC<.05 else 'none',lw=.8)
        ax.axvline(-np.log10(.05),ls='--',c=gray,lw=.7)
        ax.set(yticks=yy,yticklabels=[textwrap.fill(x,18,break_long_words=False,break_on_hyphens=False) for x in labels],xlabel=r'$-\log_{10}$(global 47-test q)',xlim=(-.08,2.55))
        ax.invert_yaxis();ax.grid(axis='x',color='#DDDDDD',lw=.4)
        panel(ax,'abc'[i],c.replace('GSE241132_D0','Normal D0'))
    fig.subplots_adjust(left=.115,right=.985,bottom=.075,top=.97)
    save(fig,'Fig1_Complete_genetic_cellular_context','Figure 1. Genetic information prioritizes cellular expression contexts. The diagram separates the fixed, RNA-independent genetic gene set from whole-cohort scoring and interpretation. All 47 original group tests are displayed: normal D0 lineages, GSE220300 lineages and nine fibroblast states, and GSE335482 lineages and six states. Squares indicate fibroblast states; circles indicate lineages. Filled symbols mark global 47-test Benjamini–Hochberg q<0.05; the dashed line is the same threshold. These are Monte Carlo associations against matched control gene sets, not patient-level susceptibility tests. The complete original family and contrary states are retained. FB01 was emphasized after examination of this exploratory family.',['groups'])
    # Reuse root-viewed Fig2 from v3; every plotted source value remains unchanged.
    for old in sorted((OUT.parent/'figures_v3').glob('Fig2_*')):
        shutil.copyfile(old,OUT/old.name)
    caption_records.append(json.loads(next(OUT.glob('Fig2_*.manifest.json')).read_text()))
    # Figure3: complete ablations and the only new scientific analysis, fixed normal temporal programmes.
    fig=plt.figure(figsize=(170/25.4,178/25.4));gs=fig.add_gridspec(2,2,height_ratios=[1.15,1],hspace=.65,wspace=.55)
    variants=['remove_ECM','remove_largest_RNA_blind_chromosome','remove_all22_published_window_union']
    assert set(data['ablation'].variant)==set(variants)
    imgs=[]
    for k,c in enumerate(['GSE220300','GSE335482']):
        ax=fig.add_subplot(gs[0,k]);d=data['ablation'].loc[data['ablation'].cohort.eq(c)]
        tab=d.pivot(index='state',columns='variant',values='combined_B30_C15_BH45_q').reindex(columns=variants).sort_index();values=-np.log10(tab.to_numpy(float))
        imgs.append(ax.imshow(values,cmap='Blues',vmin=0,vmax=2.6,aspect='auto'))
        ax.set(yticks=range(len(tab)),yticklabels=tab.index,xticks=range(3),xticklabels=['No ECM20','No chr1','No 22 windows'])
        ax.tick_params(axis='x',rotation=25)
        for iy in range(len(tab)):
            for ix in range(3):
                q=tab.iloc[iy,ix]
                if np.isfinite(q) and q<.05:ax.text(ix,iy,'•',ha='center',va='center',fontsize=9,c='white' if q<.01 else 'black')
        panel(ax,'ab'[k],c+' ablations')
    cb=fig.colorbar(imgs[0],cax=fig.add_axes([.925,.585,.012,.33]));cb.set_label(r'$-\log_{10}$(BH45 q); dot q<0.05')
    for k,pr in enumerate(['fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA']):
        ax=fig.add_subplot(gs[1,k]);d=data['time'].loc[data['time'].program.eq(pr)]
        for j,(donor,sub) in enumerate(d.groupby('patient_id',sort=True)):
            sub=sub.sort_values('day');ax.plot(range(4),sub.score,marker=['o','s','^'][j],c=[blue,orange,green][j],label=donor.split('::')[-1],ms=4)
        ax.set(xticks=range(4),xticklabels=['D0','D1','D7','D30'],xlabel='Normal repair visit (categorical)',ylabel='Within-assay mean rank fraction')
        ax.grid(axis='y',color='#DDDDDD',lw=.4);ax.legend(frameon=False,loc='best',fontsize=6)
        panel(ax,'cd'[k],'Fixed FB01-related 25' if k==0 else 'Fixed genetic 1,000 RNA proxy')
        ax.text(0,-.33,'3 donors; no P/FDR; all time points retained',transform=ax.transAxes,fontsize=6.5)
    fig.subplots_adjust(left=.11,right=.91,top=.95,bottom=.13)
    save(fig,'Fig3_Ablations_and_normal_repair','Figure 3. Genetic-set sensitivity and normal repair context. (a,b) Every existing fibroblast state is retained in three fixed ablations: removal of 20 direct ECM members, 92 chromosome-1 members, or 61 members within the 22 published windows. Color and dots use the existing joint 45-test q family; these tests are not patient comparisons. (c,d) The fixed shared 25 programme and original-Z-weighted 1,000-gene RNA proxy across D0, D1, D7 and D30 in three normal donors. All 12 existing author-labelled fibroblast pseudobulks are reused; no clustering or scDRS refit was performed. Rank(count)/G includes zeros and ties in all 26,158 measured HGNC genes. Coverage is 25/25 and 983/1,000. Visits are displayed at equal categorical spacing. These are RNA expression proxies, not changing genotypes, individual genetic risk or evidence of pathological persistence.',['ablation','time','time_summary'])
    # Figure4: all nine genetic candidates, prior sensitivity, and explicitly exploratory TF context.
    fig=plt.figure(figsize=(170/25.4,180/25.4));gs=fig.add_gridspec(2,2,height_ratios=[1,1],hspace=.6,wspace=.72)
    ax=fig.add_subplot(gs[0,0]);d=data['coloc'].copy();d['symbol']=d.gene_id.map(names);d=d.sort_values(['symbol','context']);yy=np.arange(len(d))
    ax.hlines(yy,d.low_p12_PP4,d.high_p12_PP4,color='#BBBBBB',lw=1.5)
    ax.scatter(d.primary_PP4,yy,c=blue,s=18,label='Primary prior');ax.scatter(d.low_p12_PP4,yy,c=orange,s=14,marker='s',label='Low shared prior')
    ax.axvline(.8,c=gray,lw=.7,ls='--');ax.set(xlim=(0,1.02),yticks=yy,yticklabels=[r.symbol+' | '+{'QTD000216':'FB','QTD000311':'NSkin','QTD000316':'SSkin'}[r.context] for r in d.itertuples()],xlabel='Colocalization PP4');ax.invert_yaxis();ax.legend(frameon=False,fontsize=6,loc='lower left')
    panel(ax,'a','All 11 supported units: prior sensitivity')
    ax=fig.add_subplot(gs[0,1]);d=data['expression'].loc[data['expression'].condition.eq('keloid')];tab=d.pivot(index='gene',columns='patient',values='log2_CP10K_plus0p1_delta').reindex(members)
    assert tab.shape==(9,11)
    im=ax.imshow(tab.to_numpy(float),cmap=plt.get_cmap('RdBu_r').with_extremes(bad='#BDBDBD'),norm=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=3),aspect='auto')
    ax.set(yticks=range(9),yticklabels=members,xticks=range(11),xticklabels=[z.split('::')[-1] for z in tab.columns]);ax.tick_params(axis='x',rotation=60)
    ax.text(0,-.31,'Columns: external 3 | 220 four | 335 four',transform=ax.transAxes,fontsize=6.1)
    cb=fig.colorbar(im,cax=fig.add_axes([.925,.59,.012,.33]));cb.set_label('log2(CP10K+0.1): state minus other FB')
    panel(ax,'b','All 9 candidates across 11 lesion donors')
    tf=data['TF'];tf=tf.loc[tf.kind.eq('state')&tf.TF.eq('FOSL2')&tf.condition.eq('keloid')]
    for k,method in enumerate(['ulm','mlm']):
        ax=fig.add_subplot(gs[1,k]);sub=tf.loc[tf.method.eq(method)].sort_values('sample');assert len(sub)==11
        colors=[cohort_colors[z] for z in sub.context];ax.scatter(sub.enrichment_t,range(11),c=colors,s=20,marker='o' if k==0 else 's')
        ax.axvline(0,c=gray,lw=.7);ax.set(yticks=range(11),yticklabels=[z.split('::')[-1] for z in sub['sample']],xlabel='Inferred footprint enrichment t');ax.invert_yaxis()
        panel(ax,'cd'[k],'FOSL2 after direct-target exclusions: '+method.upper())
        ax.text(0,-.25,'Exploratory context; not a genetic or causal edge',transform=ax.transAxes,fontsize=6.4)
    fig.subplots_adjust(left=.15,right=.89,top=.95,bottom=.12)
    save(fig,'Fig4_Candidate_and_regulatory_context','Figure 4. Genetic candidates and exploratory regulatory context. (a) All 11 primary-supported gene-context units, covering nine genes. Points show primary and low-shared-prior PP4, with low-to-high-prior ranges; FB, NSkin and SSkin denote cultured fibroblasts, non-sun-exposed skin and sun-exposed skin. Shared GTEx donors do not create independent replications. (b) All nine candidates across 11 lesion donors; colors show within-donor log2(CP10K+0.1) differences for FB01 or mapped FB01-like versus other fibroblasts. The color scale saturates at±3 without changing source values. (c,d) FOSL2 footprints under ULM and MLM after exclusion of direct genetic, mapping and ECM targets. These are inferred network scores, not measured TF binding or patient-level tests. FOSL2 was selected during exploration and is not the only coherent TF. No TAB2/NEDD4-to-FOSL2 causal edge is established.',['genes','coloc','expression','TF'])
    # Figure5: every prespecified intervention contrast and the full nonselectivity comparator family.
    contrasts=['TGFbeta_minus_Control','HS276_TGFbeta_minus_TGFbeta','HS276_TGFbeta_minus_Control','HC_minus_NoHC_keloid5','HC_minus_NoHC_normalScar4','HC_response_keloid_minus_normalScar','Corin_minus_DMSO_primary4']
    short=['TGFβ − Ctrl','HS276+TGFβ − TGFβ','HS276+TGFβ − Ctrl','HC − Ctrl (KL5)','HC − Ctrl (scar4)','HC response: KL − scar','Corin − DMSO (KL4)']
    fig=plt.figure(figsize=(170/25.4,190/25.4));gs=fig.add_gridspec(2,2,height_ratios=[.85,1.15],wspace=.7,hspace=.65)
    for k,pr in enumerate(['fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA']):
        ax=fig.add_subplot(gs[0,k]);d=data['drug_summary'].loc[data['drug_summary'].program.eq(pr)].set_index('contrast').reindex(contrasts)
        ax.errorbar(d['mean'],range(7),xerr=np.stack([d['mean']-d.CI95_lower,d.CI95_upper-d['mean']]),fmt='o',ms=3.5,c=blue if k==0 else orange,lw=.8,capsize=2)
        ax.axvline(0,c=gray,lw=.7);ax.set(yticks=range(7),yticklabels=short,xlabel='Within-assay programme rank change');ax.invert_yaxis();ax.grid(axis='x',color='#DDDDDD',lw=.4)
        panel(ax,'ab'[k],'Fixed 25' if k==0 else 'Fixed 1,000 RNA proxy')
    ax=fig.add_subplot(gs[1,:]);d=data['state_controls'];tab=d.pivot(index='program',columns='contrast',values='mean_delta').reindex(columns=contrasts)
    order=['ORIGINAL_SHARED_FB01_25']+[z for z in tab.index if z!='ORIGINAL_SHARED_FB01_25'];tab=tab.reindex(order);assert tab.shape==(16,7)
    im=ax.imshow(tab.to_numpy(float),cmap=plt.get_cmap('RdBu_r').with_extremes(bad='#BDBDBD'),norm=TwoSlopeNorm(vmin=-.08,vcenter=0,vmax=.08),aspect='auto')
    labels=[z.replace('ORIGINAL_SHARED_FB01_25','Original shared 25').replace('GSE220300','220').replace('GSE335482','335').replace('_',' ') for z in tab.index]
    ax.set(yticks=range(16),yticklabels=labels,xticks=range(7),xticklabels=['TGFβ−Ctrl','HS276−TGFβ','HS276−Ctrl','HC KL','HC scar','HC interaction','Corin KL']);ax.tick_params(axis='x',rotation=25)
    ax.add_patch(Rectangle((-.5,-.5),7,1,fill=False,edgecolor='black',lw=1));cb=fig.colorbar(im,cax=fig.add_axes([.90,.11,.012,.38]));cb.set_ticks([-.08,-.04,0,.04,.08]);cb.set_label('Mean rank change',fontsize=6.3,labelpad=2)
    panel(ax,'c','All 15 native-state comparators; original 25 outlined')
    fig.subplots_adjust(left=.25,right=.87,top=.95,bottom=.1)
    save(fig,'Fig5_Real_interventions_and_nonselectivity','Figure 5. Real interventions modify programmes without demonstrating FB01-selective therapy. (a,b) All seven prespecified comparisons for the shared 25 and genetic1,000 RNA programmes. Points and intervals are existing source-level means and small-sample descriptive 95% t intervals; the HC interaction compares the five keloid with four ordinary-scar source responses. Biological units are three primary foreskin sources for HS276, five and four scar strains for HC, and four primary keloid donors for Corin. They are not clinical responses. (c) Every native-state programme and the original 25 across all seven comparisons, using the existing complete descriptive table. Color scales saturate at±0.08; original values are retained. Other programmes also decline, limiting selectivity. Matching quality and programme overlap vary; the two new FB01 panels overlap 24/25 with the original and are not independent replications.',['drug_summary','drug_delta','drug_camera','drug_fry','state_controls'])
    receipt=dict(status='FIVE_INTEGRATED_FIGURES_RENDERED_VISUAL_QA_PENDING',source_sha256=sha(Path(__file__)),sources=sources,figures=caption_records,
        no_statistical_models_refit=True,no_existing_results_modified=True,physical_width_mm=170,PDF_embedded_fonts=True,PNG_dpi=600)
    (OUT/'FIGURES_RENDERED_COMPLETE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'figures':5,'directory':str(OUT)}))

if __name__=='__main__':main()
