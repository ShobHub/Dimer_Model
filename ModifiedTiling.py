import networkx as nx     
import matplotlib.pylab as plt
from scipy.spatial import KDTree
import numpy as np
import pickle

G = pickle.load(open('../Cen_Tiling_(k,R)=(20,20).txt','rb'))
pos = pickle.load(open('../Cen_Pos_(k,R)=(20,20).txt','rb'))
pos = {u:v*4 for u,v in pos.items()}
PL = pickle.load(open('../Cen_Plaqs_(k,R)=(20,20).txt','rb'))
print(len(PL))

def comE_PL(pl1,pl2):
    for u,v in pl1:
        if((u,v) in pl2 or (v,u) in pl2):
            return True
    return False

Cen = []
for i in PL:
    l = []
    xc = (pos[i[0][0]][0]+pos[i[1][0]][0]+pos[i[2][0]][0]+pos[i[3][0]][0])/ 4
    yc = (pos[i[0][0]][1] + pos[i[1][0]][1] + pos[i[2][0]][1] + pos[i[3][0]][1]) / 4
    l.append((xc, yc))
    a = np.linalg.norm(pos[i[0][0]]-pos[i[2][0]])
    b = np.linalg.norm(pos[i[1][0]]-pos[i[3][0]])
    if (a > b):
        l.append((pos[i[0][0]]-pos[i[2][0]])/a)
    else:
        l.append((pos[i[1][0]]-pos[i[3][0]])/b)
    if(abs(a-b)>3):
        l.append('Thin')
    else:
        l.append('Thick')
    Cen.append(l)
    # plt.plot(xc,yc,'ro')
    # plt.text(xc,yc,l[2])
    # plt.quiver(xc,yc,l[1][0],l[1][1],scale = 110)

Pos = {}
n = 0
E = []
cens = []
rects = {}
for i in Cen[:]:
    cens.append(i[0])
    u,v = i[0]
    if(i[2] == 'Thick'):
        l = 2.5
        w = 1.25
    else:
        l = 2
        w = 1
    Pos[4 * n] = i[0] + (l/2) * i[1] + (w/2) * np.array((i[1][1],-i[1][0]))
    Pos[4 * n + 1] = i[0] - (l/2) * i[1] + (w/2) * np.array((i[1][1],-i[1][0]))
    Pos[4 * n + 2] = i[0] - (l/2) * i[1] - (w/2) * np.array((i[1][1],-i[1][0]))
    Pos[4 * n + 3] = i[0] + (l/2) * i[1] - (w/2) * np.array((i[1][1],-i[1][0]))
    E.extend([(4*n,4*n+1),(4*n+1,4*n+2),(4*n+2,4*n+3),(4*n+3,4*n)])
    rects[i[0]] = [4*n,4*n+1,4*n+2,4*n+3]
    n += 1

Cens = cens.copy()
for i in Cens[:]:
    cens.remove(i)
    tree = KDTree(cens)
    _,iI = tree.query(i, k=4)
    ii = []
    for j in iI:
        if(comE_PL(PL[Cens.index(cens[j])],PL[Cens.index(i)])):
            ii.append(j)
    tree1 = KDTree([Pos[rects[i][0]],Pos[rects[i][1]],Pos[rects[i][2]],Pos[rects[i][3]]])
    for j in ii:
        rt = rects[cens[j]]
        dd,II = tree1.query([Pos[rt[0]],Pos[rt[1]],Pos[rt[2]],Pos[rt[3]]], k=1)
        ind = list(dd).index(min(dd))
        E.append((rects[i][II[ind]],rt[ind]))
    cens.append(i)

g = nx.Graph()
g.add_nodes_from(Pos.keys())
g.add_edges_from(E)
# pickle.dump(g,open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(20,20).txt','wb'))
# pickle.dump(Pos,open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(20,20).txt','wb'))

nx.draw(G,pos=pos,node_size=5)
nx.draw(g,pos=Pos,node_size=5,node_color='red')
plt.show()

retyreyt

import random
from networkx.algorithms import bipartite
from hopcroftkarp import HopcroftKarp
def MaxMatch(G,n,Nodes,x):
    #NO = Nodes.copy()
    def mapping():
        #random.shuffle(NO)
		#map1 = dict(map(lambda i,j : (i,j), Nodes, NO))
        map = {}
        for j in Nodes:
            flag = True
            while flag == True:
                r = random.randint(0,x)
                if(r not in map.values()):
                    map[j]=r
                    flag = False
        # for i in G.edges():  #Edges:
		# 	map[i] = random.uniform(0.1,1.1)
		# nx.set_edge_attributes(G, values=map, name='weight')
        return map
    MM=[]
    for j in range(n):
        #mapping()
        mapp = mapping()
        H = nx.relabel_nodes(G, mapp)
        X, Y = bipartite.sets(H)
        print(len(X),len(Y))
        BP = dict(map(lambda i : (i,0), X)) | dict(map(lambda i : (i,1), Y))
        BP = {}
        for i in X:
            BP[i]=0
        for i in Y:
            BP[i]=1
        nx.set_node_attributes(H, BP, name="bipartite")
        GD={}
        for line in bipartite.generate_edgelist(H, data=False):
            l=line.split(' ')
            if(int(l[0]) not in list(GD.keys())):
                GD[int(l[0])] = {int(l[1])}
            else:
                GD[int(l[0])].add(int(l[1]))
        M = HopcroftKarp(GD).maximum_matching(keys_only=True)
        # M = nx.bipartite.maximum_matching(H)
        # M = nx.bipartite.hopcroft_karp_matching(H)
        # M = nx.max_weight_matching(G,maxcardinality=True,weight='weight')
        rev_map = dict((v,k) for k,v in mapp.items())
        ME = []
        for u,v in M.items():
          ME.append((rev_map[u],rev_map[v]))
        # MM.append(ME)
        # pickle.dump(ME, open('Delta/MaxMatches_3out_NB/MaxMatch_%i.txt'%j,'wb'))
    return ME

# g = pickle.load(open('../ModifiedTiling/Mod_Tiling_(k,R)=(10,10).txt','rb'))
# Pos = pickle.load(open('../ModifiedTiling/Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
#ME = [(u,v) for u,v in nx.bipartite.hopcroft_karp_matching(g).items()]

ME = MaxMatch(g,1,list(g.nodes()),50000)

EC = []
W = []
for u,v in g.edges():
    if((u,v) in ME or (v,u) in ME):
        EC.append('red')
        W.append(1)
    else:
        EC.append('black')
        W.append(0.5)

plt.figure(dpi=200)
nx.draw(g,pos=Pos,node_size=0,node_color='blue',edge_color=EC,width=W) #,with_labels=True)
#nx.draw(G,pos=pos,node_size=3,edge_color='gray',width=0.7) #,with_labels=True,font_size=8)
plt.savefig('../ModifiedTiling/Fig_(k,R)=(20,20)_1.pdf',dpi=200)
plt.show()




