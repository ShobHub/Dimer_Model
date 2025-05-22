import Heights
import matplotlib.pylab as plt
import networkx as nx
import pickle
from collections import Counter as Ct

G = Heights.G
pos = Heights.pos
PL = Heights.PL
H = Heights.H
H = {u:v%3 for u,v in H.items()}

# obj = Heights.Heights_Penrose()
# obj.plot(H)

def center(pl):
    cx = 0
    cy = 0
    for i in pl:
        cx += pos[i[0]][0]
        cy += pos[i[0]][1]
        cx += pos[i[1]][0]
        cy += pos[i[1]][1]
    return (cx / (2 * len(pl)), cy / (2 * len(pl)))

C = Ct(H.values())
print(C)
C = [C[0]+C[1],C[1]+C[2],C[2]+C[0]]
print(C)
ind = C.index(max(C))
print(ind)

plc = []
plc_ = []
if ind == 0:
    plc = [u for u,v in H.items() if v==0]
    plc_ = [u for u,v in H.items() if v==1]
elif ind == 1:
    plc = [u for u,v in H.items() if v==1]
    plc_ = [u for u,v in H.items() if v==2]
elif ind == 2:
    plc = [u for u,v in H.items() if v==2]
    plc_ = [u for u,v in H.items() if v==0]

print(H[plc[1]],H[plc_[1]])

E = []
E_ =[]
for i in PL:
    cen = center(i)
    if cen in plc:
        E.extend(i)
    elif cen in plc_:
        E_.extend(i)

set1 = set(frozenset(edge) for edge in E)
set2 = set(frozenset(edge) for edge in E_)
CE = set1 & set2

M = []
for u,v in CE:
    M.append((u,v))

#pickle.dump(M,open('../Mij_correF/Col_RandomG_L=80.txt','wb'))
obj = Heights.Heights_Penrose()
obj.plot(H) #,CE)

dfdsfds

EC = []
W = []
for u, v in G.edges():
    if frozenset([u, v]) in CE:
        EC.append('purple')
        W.append(2)
    else:
        EC.append('gray')
        W.append(0.8)
nx.draw(G, pos, node_size=0, edge_color=EC, width=W)
plt.show()