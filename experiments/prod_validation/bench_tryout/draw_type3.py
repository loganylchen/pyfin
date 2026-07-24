import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

w=0.62; gap=0.34; x=0; xs=[]
for i in range(7):
    xs.append(x); x+=w+gap
disp_x0=x; disp_w=3.2; e8_x=disp_x0+disp_w

def exon(ax,y,xx,ww,color): ax.add_patch(Rectangle((xx,y-0.16),ww,0.32,fc=color,ec="black",lw=1.2,zorder=3))
def conn(ax,y,x0,x1): ax.plot([x0,x1],[y,y],color="0.45",lw=1.0,zorder=1)
def caret(ax,y,x0,x1):
    xm=(x0+x1)/2; ax.plot([x0,xm,x1],[y,y+0.17,y],color="0.45",lw=1.0,zorder=1)

fig,ax=plt.subplots(figsize=(13,5.6))
rows=[("READ  (ground truth)",2.0,"#2c7fb8"),
      ("candidate BEST  (has intron)",1.0,"#31a354"),
      ("candidate #2  (no intron)",0.0,"#d95f0e")]
for label,y,col in rows:
    for i in range(7):
        exon(ax,y,xs[i],w,col)
        if i<6: conn(ax,y,xs[i]+w,xs[i+1])
    conn(ax,y,xs[6]+w,disp_x0)
    exon(ax,y,e8_x,w,col)
    ax.text(-0.45,y,label,ha="right",va="center",fontsize=11,fontweight="bold")

yR,yB,yN=2.0,1.0,0.0
caret(ax,yR,disp_x0,e8_x); ax.text((disp_x0+e8_x)/2,yR+0.30,"spliced out  ~21 kb\n(these bases are NOT in the RNA)",ha="center",va="bottom",fontsize=9,color="#2c7fb8")
caret(ax,yB,disp_x0,e8_x); ax.text((disp_x0+e8_x)/2,yB+0.26,"same splice  = matches READ",ha="center",va="bottom",fontsize=9,color="#31a354")
ax.add_patch(Rectangle((disp_x0,yN-0.16),disp_w,0.32,fc="#fdd0a2",ec="#d95f0e",lw=1.2,hatch="////",zorder=3))
ax.text((disp_x0+e8_x)/2,yN-0.28,"calls it EXON  (21 kb the read never has)",ha="center",va="top",fontsize=9,color="#d95f0e")

for xg in (disp_x0,e8_x):
    ax.plot([xg,xg],[-0.9,2.75],color="0.85",lw=0.8,ls="--",zorder=0)
ax.text(disp_x0,2.66,"7,985,285",ha="center",fontsize=7.5,color="0.4")
ax.text(e8_x,2.66,"8,006,211",ha="center",fontsize=7.5,color="0.4")
ax.text(xs[0],2.66,"7 shared exons (all three identical here)",ha="left",fontsize=8,color="0.4")

msg=("M1  (mappy AS: read re-aligned to each candidate TRANSCRIPT seq):  TIE\n"
     "     -> both candidates' exons match the read equally; AS identical -> tie set\n"
     "M2  (krill signal eventalign over the junction window):  NLL 2.0470 == 2.0470,  n_events 83 == 83\n"
     "     -> byte-identical. BLIND to whether the 21 kb is intron or exon.\n"
     "READ's own genomic CIGAR:  has the N-gap 7985285->8006211  ==>  it belongs to BEST.\n"
     "     -> the one thing that KNOWS the answer is the read alignment we already have.")
ax.text(-0.45,-1.35,msg,ha="left",va="top",fontsize=9.5,family="monospace",
        bbox=dict(boxstyle="round",fc="#f7f7f7",ec="0.6"))
ax.set_xlim(-3.4,e8_x+1.3); ax.set_ylim(-2.9,3.05); ax.axis("off")
ax.set_title("Type-3 blind spot: read c37c1e06 — structures differ by 21 kb, yet M1 & M2 score them identically",
             fontsize=12,fontweight="bold")
plt.tight_layout(); plt.savefig("type3_diagram.png",dpi=135,bbox_inches="tight")
print("saved")
