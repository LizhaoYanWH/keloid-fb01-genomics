"""Portable plotting wrapper. Reads completed data; never refits scientific models."""
from pathlib import Path
import argparse,json,hashlib,re

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--registry',default='config/figure_inputs.json');ap.add_argument('--out',default='outputs/figures');ap.add_argument('--figures',type=int,nargs='+',choices=range(1,6),default=[1,2,3,4,5]);a=ap.parse_args()
 root=Path(__file__).resolve().parents[1];r=json.loads((root/a.registry).read_text());out=(root/a.out).resolve()
 # Validate all requested registry entries before creating any output.
 missing=[]
 for key,x in r['sources'].items():
  p=root/x['path']
  if not p.exists():missing.append(str(p));continue
  assert p.stat().st_size==x['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256'],key
 if missing:raise FileNotFoundError('Restore unchanged completed figure inputs before plotting:\n'+'\n'.join(missing))
 out.mkdir(parents=True,exist_ok=False)
 import pandas as pd
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6.5,'axes.titlesize':7.2,'axes.labelsize':6.5,'xtick.labelsize':5.7,'ytick.labelsize':5.7,'legend.fontsize':5.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'xtick.major.width':.5,'ytick.major.width':.5,'xtick.major.size':2.5,'ytick.major.size':2.5,'lines.linewidth':.9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','axes.unicode_minus':True})
 class Context:
  def __init__(self):self.data={k:pd.read_csv(root/x['path'],sep='\t',float_precision='round_trip') for k,x in r['sources'].items()}
  def panel(self,ax,letter,title):
   ax.set_title(title,loc='left',pad=5,fontsize=7,fontweight='medium');offset=(3.5/(ax.figure.get_figwidth()*25.4))/ax.get_position().width;ax.text(-offset,1.04,letter.upper(),transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom',ha='left')
  def save(self,fig,stem,caption,sourcekeys):
   assert set(sourcekeys)<=set(r['sources']);fig.canvas.draw()
   for ext in ['pdf','svg','png']:fig.savefig(out/(stem+'.'+ext),dpi=600 if ext=='png' else 300,bbox_inches=None,facecolor='white')
   (out/(stem+'.caption.md')).write_text(caption+'\n');plt.close(fig)
 from figures_12 import render_fig1,render_fig2
 from figure_3 import render_fig3
 from figures_45 import render_fig4,render_fig5
 ctx=Context();f={1:render_fig1,2:render_fig2,3:render_fig3,4:render_fig4,5:render_fig5}
 for i in a.figures:f[i](ctx)
if __name__=='__main__':main()
