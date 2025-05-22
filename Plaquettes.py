import networkx.algorithms.isomorphism as iso
from scipy.spatial import KDTree
import matplotlib.pylab as plt
import networkx as nx
import numpy as np
import pickle

g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(20,20).txt','rb'))
pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(20,20).txt','rb'))
G = pickle.load(open('../Cen_Tiling_(k,R)=(20,20).txt','rb'))
Pos = pickle.load(open('../Cen_Pos_(k,R)=(20,20).txt','rb'))
Pos = {u:v*4 for u,v in Pos.items()}
PL_e = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(20,20).txt','rb'))

PL_BP = {}
PL = []
tree = KDTree(list(pos.values()))
for i in G.nodes():
    k = G.degree(i)
    dd, ii = tree.query(Pos[i], distance_upper_bound=3, k=k*2+4)
    ii = np.delete(np.array(ii),np.where(np.isinf(dd) | np.isnan(dd))[0]).tolist()
    gc = nx.induced_subgraph(g, ii)
    CC = list(nx.chordless_cycles(gc,k*2))
    if CC !=[]:
        lens = [len(j) for j in CC]
        gc = nx.induced_subgraph(gc,CC[lens.index(max(lens))])
        if(len(gc.nodes()) !=4):
            PL.append(list(gc.edges()))
            PL_BP[i] = list(gc.edges())

PL_e.extend(PL)
print(PL_e)
print(PL_BP)
pickle.dump(PL_BP, open("../ModifiedTiling/Cen_PlaqBPCharge_(k,R)=(20,20).txt", "wb"))
pickle.dump(PL_e, open("../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(20,20).txt", "wb"))

gr=nx.Graph()
for i in PL_e:
    gr.add_edges_from(i)

print(nx.is_isomorphic(gr,g))

nx.draw(gr, pos, node_size=5)
plt.show()

fdsdgfs

bn,tn = nx.bipartite.sets(G)
s = 0
print(len(bn),len(tn))
for i in bn:
    s += G.degree(i)
print(len(bn)+(s/2))
nx.draw(G,pos=Pos,node_size = 0) #,with_labels=True,font_size=8)
nx.draw(g,pos=pos,node_size = 5,edge_color='b')#,with_labels=True,font_size=8)
plt.show()


