#target illustrator
// Standalone Fig1A has no panel letter. All six labels are added only at final composition.
(function(){
var D='/Users/lizhaoyan/Documents/文稿 - Lizhao的MacBook Pro/ChatGPT/生信/Five_Studies_Execution_20261020/projects/N02_keloid_genetic_cell_state/figure_design_20261022_v1';
$.evalFile(new File(D+'/source/illustrator_helpers_v1.jsxinc'));
var oldUI=app.userInteractionLevel;app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
function newLabel(doc,layer,s,left,top){var H=doc.artboards[0].artboardRect[1],t=n02Text(layer,s,left,0,9,true,'000000',H,false),a=t.textRange.characterAttributes;t.name='Panel '+s+' | code-assembled Arial Bold 9pt';a.horizontalScale=100;a.verticalScale=100;a.tracking=0;a.baselineShift=0;t.left=left;t.top=top;return t;}
try{
// Re-open preserved accepted source, never reuse the failed, unsaved v5 in-memory attempt.
var full=app.open(new File(D+'/exports/Figure_1_revised_A_sunburst_v3.ai'));full.activate();
n02Audit(full,D+'/qa/Figure1_native_audit_before_panel_fix_v6.json','Before code-assembled uniform panel labels','Preserved v3 source');
var targets=[],hiddenA=[],i;
for(i=0;i<full.textFrames.length;i++){var t=full.textFrames[i],a=t.textRange.characterAttributes,s=t.contents;
 if(s=='A'&&t.layer.name=='04 New A native Arial labels'&&a.size==11)targets.push({s:s,left:t.visibleBounds[0],top:t.visibleBounds[1],t:t});
 else if(/^[B-F]$/.test(s)&&a.textFont.name=='DejaVuSans-Bold'&&a.size==9)targets.push({s:s,left:t.visibleBounds[0],top:t.visibleBounds[1],t:t});
 else if(s=='A'&&t.layer.name=='01 Original B–F embedded PDF'&&a.size==9)hiddenA.push(t);
}
if(targets.length!=6||hiddenA.length!=1)throw new Error('Unexpected visible/hidden panel identifiers '+targets.length+'/'+hiddenA.length);
var originalLayer=full.layers.getByName('01 Original B–F embedded PDF');originalLayer.locked=false;
var layer=n02Layer(full,'05 Assembly labels A–F | Arial Bold 9pt | scale 100%');
for(i=0;i<targets.length;i++){var z=targets[i];newLabel(full,layer,z.s,z.left,z.top);z.t.remove();}
hiddenA[0].remove();originalLayer.locked=true;
n02Audit(full,D+'/qa/Figure1_Adobe_native_audit_v6.json','A-F added by composition code: Arial-BoldMT 9pt, scale100%, black','Original science unchanged except replacement of panel identifiers');
n02Export(full,D+'/exports/Figure_1_revised_A_sunburst_v6',3840);
var single=app.open(new File(D+'/exports/Fig1A_Adobe_sunburst_v1.ai'));single.activate();
n02Audit(single,D+'/qa/Fig1A_native_audit_before_panel_fix_v3.json','Before removal of standalone A label','Preserved standalone v1');
var sa=[];for(i=0;i<single.textFrames.length;i++)if(single.textFrames[i].contents=='A')sa.push(single.textFrames[i]);
if(sa.length!=1)throw new Error('Standalone A match ambiguous');sa[0].remove();
n02Audit(single,D+'/qa/Fig1A_Adobe_native_audit_v3.json','Standalone Fig1A without panel letter; identifiers belong to final composition','Same Sunburst artwork, native text and arrows');
n02Export(single,D+'/exports/Fig1A_Adobe_sunburst_v3',3840);
n02Write(D+'/qa/PANEL_LETTERS_V6_COMPLETE.txt','Adobe native composition complete. Standalone Fig1A has no A. Final Figure1 A-F each created by the same code/style: Arial-BoldMT 9pt, scale100%, black. All non-panel labels/data retained. Hidden old A removed. Original files retained.');full.activate();app.redraw();
}catch(e){n02Write(D+'/qa/PANEL_LETTERS_V6_ERROR.txt',String(e)+' line '+e.line);}finally{app.userInteractionLevel=oldUI;}
})();
