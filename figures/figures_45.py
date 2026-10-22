"""Redraw completed candidate, regulatory and intervention results; no tests/refits."""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

GENES = ['TAB2','NEDD4','ASIP','LSP1','TNNT3','SETBP1','CARS1-AS1','SETBP1-DT','PRR33']
COHORTS = ['GSE220300','GSE335482','GSE163973']
CONTRASTS = ['TGFbeta_minus_Control','HS276_TGFbeta_minus_TGFbeta','HS276_TGFbeta_minus_Control','HC_minus_NoHC_keloid5','HC_minus_NoHC_normalScar4','HC_response_keloid_minus_normalScar','Corin_minus_DMSO_primary4']
SHORT = ['TGFβ vs Ctrl','HS276 add-on','HS276 vs Ctrl','HC: keloid','HC: ordinary scar','HC interaction','Corin: keloid']
PROGRAMS = ['fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA','fixed_Reactome_ECM_comparator','fixed_Reactome_CellCycle_comparator']
P_NAMES = ['Shared 25','MAGMA1000 RNA','ECM comparator','Cell-cycle comparator']
BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#009E73'; GREY='#777777'

def axat(fig, x,y,w,h):
    return fig.add_axes([x/170,y/190,w/170,h/190])

def clean_heat(ax,tab,cmap='RdBu_r',vmin=-3,vmax=3,fontsize=5.7):
    cm=plt.get_cmap(cmap).with_extremes(bad='#D9D9D9')
    a=tab.to_numpy(float)
    im=ax.imshow(a,aspect='auto',interpolation='nearest',cmap=cm,vmin=vmin,vmax=vmax)
    ax.set_yticks(np.arange(len(tab)));ax.set_yticklabels(tab.index,fontsize=fontsize)
    ax.set_xticks(np.arange(len(tab.columns)));ax.set_xticklabels(tab.columns,fontsize=fontsize)
    ax.tick_params(length=2,pad=1.8)
    for sp in ax.spines.values():sp.set_visible(False)
    return im

def cbar(fig,im,rect,label,ticks=None):
    ca=fig.add_axes([rect[0]/170,rect[1]/190,rect[2]/170,rect[3]/190])
    cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=ticks)
    cb.ax.tick_params(labelsize=5.3,length=2,pad=1)
    cb.set_label(label,fontsize=5.6,labelpad=1)
    return cb

def render_fig4(ctx):
    d=ctx.data
    fig=plt.figure(figsize=(170/25.4,190/25.4))
    # Complete primary-supported units; three existing priors, not confidence intervals.
    ax=axat(fig,24,137,54,43)
    names=dict(zip(d['genes'].gene_id,d['genes'].approved_symbol))
    q=d['coloc'].copy();q['symbol']=q.gene_id.map(names);q=q.sort_values(['symbol','context'])
    yy=np.arange(len(q));ax.hlines(yy,q.low_p12_PP4,q.high_p12_PP4,color='#B9B9B9',lw=1.4)
    ax.scatter(q.primary_PP4,yy,c=BLUE,s=13,label='Primary')
    ax.scatter(q.low_p12_PP4,yy,c=ORANGE,s=11,marker='s',label='Lower p12')
    ax.axvline(.8,c=GREY,ls='--',lw=.6)
    lab=[r.symbol+' / '+{'QTD000216':'FB','QTD000311':'NSkin','QTD000316':'SSkin'}[r.context] for r in q.itertuples()]
    ax.set(yticks=yy,yticklabels=lab,xlim=(0,1.02),xlabel='Colocalization PP4');ax.invert_yaxis()
    ax.tick_params(axis='y',labelsize=5.7)
    ax.legend(frameon=False,fontsize=5.5,ncol=2,loc='lower left',bbox_to_anchor=(-.26,-.54),handletextpad=.3,columnspacing=.7)
    ctx.panel(ax,'a','Candidate support depends on prior')
    # Every candidate and every lesion donor, including all reversals and NA.
    ax=axat(fig,98,137,65,43)
    sub=d['expression'].loc[d['expression'].condition.eq('keloid')]
    pats=sorted(sub.patient.unique());tab=sub.pivot(index='gene',columns='patient',values='log2_CP10K_plus0p1_delta').reindex(index=GENES,columns=pats)
    assert tab.shape==(9,11)
    im=clean_heat(ax,tab);ax.set_xticklabels([p.split('::')[-1] for p in pats],rotation=60,ha='right')
    for edge in [2.5,6.5]:ax.axvline(edge,color='white',lw=1.4)
    cbar(fig,im,(106,122,50,1.5),'State − other FB: log2(CP10K+0.1)',[-3,0,3])
    ctx.panel(ax,'b','Nine candidates in eleven lesion donors')
    tf=d['TF'].loc[d['TF'].kind.eq('state')&d['TF'].condition.eq('keloid')]
    # The complete fixed15, not only the coherent FOSL2 row.
    tfs=sorted(tf.TF.unique());assert len(tfs)==15
    for k,method in enumerate(['ulm','mlm']):
        ax=axat(fig,21+85*k,75,60,34)
        sub=tf.loc[tf.method.eq(method)]
        tab=sub.pivot(index='TF',columns='sample',values='enrichment_t').reindex(index=tfs,columns=pats)
        im=clean_heat(ax,tab,vmin=-5,vmax=5,fontsize=5.5)
        ax.set_xticklabels([p.split('::')[-1] for p in pats],rotation=60,ha='right')
        iy=tfs.index('FOSL2');ax.add_patch(Rectangle((-.5,iy-.5),11,1,fill=False,lw=.8,ec='black'))
        for edge in [2.5,6.5]:ax.axvline(edge,color='white',lw=1.1)
        ctx.panel(ax,'cd'[k],method.upper()+': all fixed regulatory footprints')
        cbar(fig,im,(36+85*k,61,33,1.5),'Footprint enrichment t (clipped ±5)',[-5,0,5])
    # Entire 72-row original request under each of its distinct coordinates.
    gf=d['GF'];ops=['delete_natural','delete_position_preserved']
    columns=[(s,c,o) for s in ['GSE220300','GSE335482'] for c in ['FB01','otherFB'] for o in ops]
    col_labels=[('220' if s=='GSE220300' else '335')+' '+('F1' if c=='FB01' else 'other')+'\n'+('del' if o==ops[0] else 'pos') for s,c,o in columns]
    for k,(field,title) in enumerate([('mean3seed_raw_S','Virtual deletion: raw projection'),('mean3seed_rank4_mean_centered_S','Virtual deletion: rank-four centered')]):
        ax=axat(fig,21+85*k,22,60,30)
        tab=gf.pivot(index='gene',columns=['study','context','operator'],values=field).reindex(index=GENES,columns=pd.MultiIndex.from_tuples(columns))
        a=tab.to_numpy(float);bound=float(np.nanmax(np.abs(a)))
        assert np.isfinite(bound) and bound>0
        im=clean_heat(ax,tab,vmin=-bound,vmax=bound)
        ax.set_xticklabels(col_labels,fontsize=5.1,rotation=40,ha='right')
        for iy,ix in zip(*np.where(~np.isfinite(a))):ax.text(ix,iy,'NA',ha='center',va='center',fontsize=4.7,color='#444444')
        for x in [1.5,3.5,5.5]:ax.axvline(x,color='white',lw=1.1)
        ctx.panel(ax,'ef'[k],title)
        exponent=int(np.floor(np.log10(bound)))
        cb=cbar(fig,im,(36+85*k,10.5,33,1.4),f'Embedding S (×10^{exponent}); 3-seed mean',[-bound,0,bound])
        cb.formatter=FuncFormatter(lambda x,pos,e=exponent:f'{x*10**(-e):.1f}');cb.update_ticks()
        cb.ax.xaxis.get_offset_text().set_visible(False)
    cap=('Figure 4. Genetic candidates, regulatory context and reference-dependent virtual perturbations. '
         '(a) All eleven primary-supported gene–context units; points show primary and lower shared priors and lines the original low-to-high-prior range, not a confidence interval. '
         '(b) All nine candidates in three external, four GSE220300 and four GSE335482 lesion donors; values saturate at ±3. '
         '(c,d) All fifteen fixed TFs under ULM and MLM after genetic/mapping/ECM-target exclusions, with FOSL2 outlined; inferred enrichment t is clipped at ±5 and is not TF binding or a patient-level P value. '
         '(e,f) All nine candidates, two cohorts, FB01/otherFB contexts and natural/position-preserved deletions from the complete 72-row request. Raw and rank-four mean-centered embedding projections along the lesional otherFB-minus-FB01 axis have separate symmetric scales. Gray/NA denotes unavailable full-four-patient support or unmeasurable expression, not biological zeros. Three seeds are model repeats, not biological replicates. No RNA restoration, selective dependency or directed candidate-to-TF treatment mechanism is established. del, natural deletion; pos, position-preserved deletion; F1, FB01.')
    ctx.save(fig,'Figure_4',cap,['genes','coloc','expression','TF','GF'])

def render_fig5(ctx):
    d=ctx.data;fig=plt.figure(figsize=(170/25.4,190/25.4))
    for k,pr in enumerate(PROGRAMS[:2]):
        ax=axat(fig,29+85*k,134,51,44)
        means=d['drug_summary'].loc[d['drug_summary'].program.eq(pr)].set_index('contrast').reindex(CONTRASTS)
        raw=d['drug_delta'].loc[d['drug_delta'].program.eq(pr)]
        for iy,co in enumerate(CONTRASTS):
            a=raw.loc[raw.contrast.eq(co)].sort_values('source_strain')
            offsets=np.linspace(-.13,.13,len(a)) if len(a)>1 else np.zeros(len(a))
            ax.scatter(a.delta,iy+offsets,s=8,facecolor='white',ec='#555555',lw=.45,zorder=3)
        ax.errorbar(means['mean'],np.arange(7),xerr=np.stack([means['mean']-means.CI95_lower,means.CI95_upper-means['mean']]),fmt='o',color=[BLUE,ORANGE][k],ms=3,lw=.8,capsize=1.8,zorder=4)
        ax.axvline(0,c=GREY,lw=.6);ax.set(yticks=np.arange(7),yticklabels=SHORT,xlabel='Within-assay rank change');ax.invert_yaxis();ax.grid(axis='x',color='#E5E5E5',lw=.35)
        ax.tick_params(axis='y',labelsize=5.5)
        ctx.panel(ax,'ab'[k],['Fixed shared 25: source-level responses','Fixed MAGMA1000: RNA response'][k])
        ax.text(0,-.29,'○ sources    ● mean + 95% interval',fontsize=5.1,transform=ax.transAxes)
    # Two original BH14 families remain primary; joint BH28 is sensitivity only.
    cmap=plt.get_cmap('RdBu_r').with_extremes(bad='#D9D9D9')
    for k,(key,qkey,name) in enumerate([('drug_camera','BH28_joint','CAMERA: competitive gene-set test'),('drug_fry','BH28_directional','FRY: self-contained directional test')]):
        qsource=d[key].copy()
        qsource['planned_q']=np.where(qsource.program.isin(PROGRAMS[:2]),qsource.BH14_primary2,qsource.BH14_ECM_growth_comparator)
        tab=qsource.pivot(index='program',columns='contrast',values='planned_q').reindex(index=PROGRAMS,columns=CONTRASTS)
        joint=qsource.pivot(index='program',columns='contrast',values=qkey).reindex(index=PROGRAMS,columns=CONTRASTS)
        direction=d[key].pivot(index='program',columns='contrast',values='Direction').reindex(index=PROGRAMS,columns=CONTRASTS)
        assert tab.shape==(4,7) and tab.notna().all().all()
        signed=-np.log10(tab.to_numpy(float))*np.where(direction.eq('Up'),1,-1)
        ax=axat(fig,29+85*k,91,51,20)
        im=ax.imshow(signed,cmap=cmap,vmin=-4,vmax=4,aspect='auto',interpolation='nearest')
        ax.set(yticks=range(4),yticklabels=P_NAMES,xticks=range(7),xticklabels=['TGFβ','HS add','HS vs C','HC KL','HC scar','HC int','Corin']);ax.tick_params(axis='y',labelsize=5.1,pad=1.5);ax.tick_params(axis='x',labelsize=4.9,rotation=45,pad=1.5)
        for iy in range(4):
            for ix in range(7):
                if tab.iloc[iy,ix]<.05:ax.text(ix,iy,'•',ha='center',va='center',fontsize=8,c='white' if abs(signed[iy,ix])>2 else 'black')
                elif joint.iloc[iy,ix]<.05:ax.scatter(ix,iy,marker='D',s=13,facecolors='none',edgecolors='black',linewidths=.65)
        ax.axhline(1.5,color='white',lw=1.1)
        ctx.panel(ax,'cd'[k],name)
        cbar(fig,im,(40+85*k,75,32,1.5),'Direction × −log10(planned BH14 q)',[-4,0,4])
    # The original16-state and 9-candidate families add biological contrast, not another fit.
    ax=axat(fig,29,24,68,41)
    tab=d['state_controls'].pivot(index='program',columns='contrast',values='mean_delta').reindex(columns=CONTRASTS)
    order=['ORIGINAL_SHARED_FB01_25']+[p for p in tab.index if p!='ORIGINAL_SHARED_FB01_25'];tab=tab.reindex(order)
    assert tab.shape==(16,7)
    im=ax.imshow(tab.to_numpy(float),cmap=cmap,norm=TwoSlopeNorm(vmin=-.08,vcenter=0,vmax=.08),aspect='auto',interpolation='nearest')
    labels=[x.replace('ORIGINAL_SHARED_FB01_25','Shared 25').replace('GSE220300','220').replace('GSE335482','335').replace('::',' ') for x in order]
    ax.set(yticks=range(16),yticklabels=labels,xticks=range(7),xticklabels=['TGFβ','HS add','HS vs C','HC KL','HC scar','HC int','Corin']);ax.tick_params(axis='y',labelsize=5.1,pad=1.5);ax.tick_params(axis='x',labelsize=5.1,rotation=45,pad=1.5)
    for iy,x in enumerate(order):
        if x=='ORIGINAL_SHARED_FB01_25' or x.endswith('::FB01'):ax.add_patch(Rectangle((-.5,iy-.5),7,1,fill=False,ec='black',lw=.55))
    ctx.panel(ax,'e','All native-state programs: broad response')
    cbar(fig,im,(44,10,43,1.5),'Mean program rank change',[-.08,0,.08])
    ax=axat(fig,117,24,46,41)
    candidates=d['drug_candidates']
    tab=candidates.pivot(index='approved_symbol',columns='contrast',values='logFC').reindex(index=GENES,columns=CONTRASTS)
    qtab=candidates.pivot(index='approved_symbol',columns='contrast',values='BH63_candidate_family').reindex(index=GENES,columns=CONTRASTS)
    im=clean_heat(ax,tab,vmin=-2,vmax=2,fontsize=5.4)
    ax.set_xticklabels(['TGFβ','HS add','HS vs C','HC KL','HC scar','HC int','Corin'],rotation=45,ha='right',fontsize=5.1)
    for iy in range(9):
        for ix in range(7):
            if np.isfinite(qtab.iloc[iy,ix]) and qtab.iloc[iy,ix]<.05:ax.text(ix,iy,'•',ha='center',va='center',fontsize=6,c='white' if abs(tab.iloc[iy,ix])>1 else 'black')
    ctx.panel(ax,'f','Candidate-gene RNA responses')
    cbar(fig,im,(124,10,32,1.5),'log2 fold change; dot: BH63 q<0.05',[-2,0,2])
    cap=('Figure 5. Controlled treatments modify programs without establishing selective therapy. '
         '(a,b) All seven fixed contrasts for the shared 25 and MAGMA1000 RNA proxy. Open points are the original independent source deltas; filled points and intervals are existing means and descriptive 95% t intervals. '
         'HS276 uses three primary foreskin sources; hydrocortisone (HC) uses five keloid and four ordinary-scar strains; Corin uses four primary keloid sources. The HC interaction is a difference of source-group responses, without invented individual interaction observations. '
         '(c,d) All four fixed programs under CAMERA and FRY. Color uses each original BH14 family (primary two programs or two comparators), saturating at ±4; dots indicate planned q<0.05, open diamonds joint-BH28-only support. Primary decisions are unchanged; tests have distinct null hypotheses. '
         '(e) All sixteen previously fixed programs, including fifteen native-state comparators; outlined rows identify the shared25 and both FB01 programs. Colors saturate at ±0.08. '
         '(f) All nine candidates under all seven contrasts; existing log2 fold changes saturate at ±2 and dots use the complete BH63 candidate family. Gray denotes unavailable values. Full numeric results remain in source tables. These are culture RNA responses, not clinical benefit, direct target engagement or FB01-selective restoration.')
    ctx.save(fig,'Figure_5',cap,['drug_summary','drug_delta','drug_camera','drug_fry','state_controls','drug_candidates'])
