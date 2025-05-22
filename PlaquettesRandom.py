import networkx.algorithms.isomorphism as iso
import matplotlib.pylab as plt
from math import atan2, degrees
import networkx as nx
import numpy as np
import operator
import pickle
import random

g = pickle.load(open('../RandomG/G_L=80_X1.txt','rb'))
pos = pickle.load(open('../RandomG/Pos_L=80_X1.txt','rb'))
#EB = pickle.load(open('../RandomG/PbcE_L=80_.txt','rb'))
#g.remove_edges_from(list(EB.keys()))

print(nx.is_bipartite(g))
print(nx.is_planar(g))
pos = nx.planar_layout(g)
nx.draw(g,pos,node_size=5) #,with_labels=True)

sortE = {}
for i in g.nodes():
    sortE[i] = []
    E = {}
    for j in list(g.neighbors(i)):
        angle = degrees(atan2(pos[j][1]-pos[i][1],pos[j][0]-pos[i][0]))
        E[(angle + 360) % 360] = (i,j)
    for j in sorted(E):
        sortE[i].append(E[j])
    sortE[i].append(sortE[i][0])

PL = []
E = list(g.edges()).copy()
while E != []:
    a,b = E[0]
    j = a
    pl = []
    while b != a:
        pl.append((j, b))
        j,b = sortE[b][sortE[b].index((b,j))+1]
    pl.append((j,b))
    pl = sorted(pl, key=operator.itemgetter(0))
    if pl not in PL:
        PL.append(pl)
        print(len(pl),pl)
    E.pop(0)
print(len(PL))
pickle.dump(PL,open('../RandomG/Plaqs_L=80_X.txt','wb'))
plt.show()

dfgfdgfd



