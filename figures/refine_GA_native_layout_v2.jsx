#target illustrator
(function(){
var ROOT='/Users/lizhaoyan/Documents/文稿 - Lizhao的MacBook Pro/ChatGPT/生信/Five_Studies_Execution_20261020/projects/N02_keloid_genetic_cell_state/figure_design_20261022_v1';
var old=app.userInteractionLevel;
$.evalFile(new File(ROOT+'/source/illustrator_helpers_v1.jsxinc'));
try{
app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
var d=null;for(var i=0;i<app.documents.length;i++){try{if(app.documents[i].fullName.fsName===ROOT+'/exports/Graphical_Abstract_Adobe_sunburst_v1.ai')d=app.documents[i];}catch(ignore){}}
if(!d)d=app.open(new File(ROOT+'/exports/Graphical_Abstract_Adobe_sunburst_v1.ai'));
var changes=0;
for(i=0;i<d.textFrames.length;i++){
var t=d.textFrames[i],s=t.contents;
if(s==='Retained after ECM and locus exclusions'||s==='A reproducible cellular context\rfor scar-associated genetic signals'){t.left+=640*.019;changes++;}
if(s==='Candidate signals: TGFB1 / IGF2'||s==='Regulatory prediction: FOSL2'||s==='Non-ECM Up output reduced in silico'){t.top-=360*.020;changes++;}
}
if(changes!==5)throw new Error('Expected 5 targeted text spacing changes; got '+changes);
n02Audit(d,ROOT+'/qa/Graphical_Abstract_Adobe_native_audit_v2.json','graphical abstract spacing refinement',ROOT+'/assets/Graphical_Abstract_02_population_textless_sunburst_v4.png');
n02Export(d,ROOT+'/exports/Graphical_Abstract_Adobe_sunburst_v2',3840);
d.activate();app.redraw();n02Write(ROOT+'/qa/GA_LAYOUT_V2_COMPLETE.txt','Five native text spacing changes only.\nNo artwork regeneration or replacement.');
}catch(e){n02Write(ROOT+'/qa/GA_LAYOUT_V2_ERROR.txt',String(e)+'\nline='+e.line);}finally{app.userInteractionLevel=old;}
})();
