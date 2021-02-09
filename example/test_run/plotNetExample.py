import networkx as nx
from networkx.algorithms import bipartite
from networkx.drawing.layout import bipartite_layout
from cLoops2.settings import *
import matplotlib as mpl
from matplotlib.lines import Line2D


def getNet(f, chrom=None):
    G = nx.Graph()
    for line in open(f):
        line = line.split("\n")[0].split("\t")
        c = line[0].split("|")[0].split(":")[0]
        if chrom is not None and c != chrom:
            continue
        G.add_edge(line[0], line[2], type=line[1])
    for node in G.nodes():
        n = node.split("|")[-1]
        G.nodes[node]["type"] = n
    return G


def plotMostConnected(f, c="chr21"):
    n = f.split("/")[-1].split(".sif")[0]
    G = getNet(f, chrom=c)
    sgs = list(nx.connected_components(G))
    mg = None
    nodes = 0
    for sg in sgs:
        if len(sg) > nodes:
            mg = sg
            nodes = len(sg)
    mg = G.subgraph(mg)
    start = None
    end = None
    for tmp in list(mg.nodes()):
        sg2 = tmp.split("|")[0].split(":")[1].split("-")
        if start is None and end is None:
            start = int(sg2[0])
            end = int( sg2[1] )
        else:
            if int(sg2[0]) < start:
                start = int(sg2[0])
            if int(sg2[1]) > end:
                end = int( sg2[1] )
    ps = []
    es = []
    for node in mg.nodes:
        if mg.nodes[node]["type"] == "Promoter":
            ps.append(node)
        else:
            es.append(node)
    labels= { node:node.split("|")[1] for node in mg.nodes()}
    pos = nx.spring_layout(mg)
    #pos = nx.circular_layout(mg)
    nx.draw_networkx_nodes(mg,
                           pos,
                           nodelist=ps,
                           node_color=colors[0],
                           node_size=20)
    nx.draw_networkx_nodes(mg,
                           pos,
                           nodelist=es,
                           node_color=colors[1],
                           node_size=20)
    nx.draw_networkx_edges(mg,
                           pos,
                           edgelist=list(mg.edges()),
                           width=1,
                           edge_color="gray")   
    #nx.draw_networkx_labels(mg,pos,labels,font_size=6)
    pylab.axis("off")
    pylab.title( "%s:%s-%s"%(c,start,end))
    ax = mpl.pyplot.gca()
    legele = [Line2D([0],[0],marker='o',color=colors[0],label="promoter",markersize=5),
              Line2D([0],[0],marker='o',color=colors[1],label="enhancer",markersize=5),]

    ax.legend(handles=legele,loc="best")
    pylab.savefig(n + "_largest_components.pdf")


plotMostConnected("gm_loops_ep_net.sif","chr21")
