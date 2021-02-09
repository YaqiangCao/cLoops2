import pandas as pd
from matplotlib.patches import Arc
import matplotlib.patches as patches
from cLoops2.settings import *

mat = pd.read_csv("stat.txt",index_col=0,sep="\t")
s = mat["Trac"]
m = s.min()
fig, ax = pylab.subplots(figsize=(3,2))
ratio = s/s.sum()*100

#enhancer-other
a = 1 
b = 0.5
c = 0.5
k = "Enhancer-None"
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.25, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)

#enhancer-enhancer
k = "Enhancer-Enhancer"
c = 1.5 
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.2, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)

#enhancer-enhancer
k = "Enhancer-Promoter"
c = 2.5 
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.2, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)

#promoter-promoter
k = "Promoter-Promoter"
c = 3.5 
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.1, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)

#promoter-promoter
k = "Promoter-None"
c = 4.5 
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.2, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)


#none-none
k = "None-None"
a = 5
b = 1.2
c = 2.5 
ax.add_patch(
    Arc(
        (c, 0.2),
        a,
        b,
        theta1=0,
        theta2=180,
        edgecolor=colors[1],
        lw = s[k]/m/5
    ))
ax.text(c-0.2, b / 2, "%.2f%s"%(ratio[k],"%"), fontsize=6)

#plot elements
p = patches.Rectangle((-0.5,0.02),6.0,0.06, fill=True, color="k",alpha=0.1)
ax.add_patch(p)
p = patches.Rectangle( (-0.1,0.0),0.2,0.1,fill=True,color="k",alpha=0.8 )
ax.add_patch(p)
p = patches.Rectangle( (0.9,0.0),0.2,0.1,fill=True,color=colors[0],alpha=0.8 )
ax.add_patch(p)
p = patches.Rectangle( (1.9,0.0),0.2,0.1,fill=True,color=colors[0],alpha=0.8 )
ax.add_patch(p)
p = patches.Rectangle( (2.9,0.0),0.2,0.1,fill=True,color=colors[2],alpha=0.8 )
ax.add_patch(p)
p = patches.Rectangle( (3.9,0.0),0.2,0.1,fill=True,color=colors[2],alpha=0.8 )
ax.add_patch(p)
p = patches.Rectangle( (4.9,0.0),0.2,0.1,fill=True,color="k",alpha=0.8 )
ax.add_patch(p)
#annotate elements
ax.text(-0.3,-0.1,"Other",fontsize=8)
ax.text(1.1,-0.1,"Enhancer",fontsize=8)
ax.text(3.1,-0.1,"Promoter",fontsize=8)
ax.text(4.7,-0.1,"Other",fontsize=8)
#basic settings
ax.set_xlim([-1,6])
pylab.axis("off")
pylab.savefig("Trac_ep_stat.pdf")

