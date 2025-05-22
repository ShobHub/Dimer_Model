import matplotlib.pylab as plt
import networkx as nx
import numpy as np
import pickle

g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))
PL_BP = pickle.load(open('../ModifiedTiling/Cen_PlaqBPCharge_(k,R)=(10,10).txt','rb'))

Samples = pickle.load(open('../Mij_correF/Data_ModP/Samples/T=0.02_V=%.5f.txt' % -0.008, 'rb'))

DistV = {}
for i in list(g.edges())[::30]:
    cen1 = (pos[i[0]] + pos[i[1]]) / 2
    for j in list(g.edges())[::30]:
        cen2 = (pos[j[0]] + pos[j[1]]) / 2
        d = int(np.linalg.norm(cen1-cen2))
        if d!=0 and d not in DistV.keys():
            DistV[d] = [(i,j)]
        elif d in DistV.keys() and (j,i) not in DistV[d]:
            DistV[d].append((i,j))

def is_tuple_in_list(t, list_of_tuples):
    for tup in list_of_tuples:
        if frozenset(tup) == frozenset(t):
            return True
    return False

print(Samples[1]==Samples[-1],Samples[-1]==Samples[0],Samples[0]==Samples[25])

G = {}
for r in DistV.keys():
    s = []
    for i in Samples:
        c = 0
        for j in DistV[r]:
            if is_tuple_in_list(j[0], i) and is_tuple_in_list(j[1], i):
                c += 1
        s.append(c/len(DistV[r]))
    G[r] = np.mean(s)   #((-1)**r)*

print(G)
filter_G = {k: v for k, v in G.items() if v != 0}
sort_G = dict(sorted(filter_G.items()))
print(sort_G)

plt.scatter(sort_G.keys(),sort_G.values(),s=5)
plt.plot(sort_G.keys(),sort_G.values())
plt.xscale('log')
plt.yscale('log')
plt.show()

# nx.draw(g,pos,node_size=0,with_labels=True)
# plt.show()