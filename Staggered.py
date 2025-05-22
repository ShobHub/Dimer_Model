import pickle
import networkx as nx
import numpy as np
from math import degrees,atan2
import matplotlib.pylab as plt
import matplotlib.path as mpltPath
from shapely.geometry import Point, Polygon

g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))

# g = pickle.load(open('../RandomG/G_L=80_tutte.txt','rb'))
# pos = pickle.load(open('../RandomG/Pos_L=80_tutte.txt','rb'))
# PL_HC = pickle.load(open('../RandomG/Plaqs_L=80_tutte.txt','rb'))

def arrange_cycle_in_order(cycle_edges):
    current_edge = cycle_edges.pop(0)[::-1]
    ordered_edges = [current_edge]
    while cycle_edges:
        last_node = ordered_edges[-1][1]
        for i, edge in enumerate(cycle_edges):
            if last_node in edge:
                if edge[1] == last_node:
                    ordered_edges.append(edge[::-1])
                else:
                    ordered_edges.append(edge)
                cycle_edges.pop(i)
                break
    return ordered_edges

def Cen_Plaq(str,thresh):
    C = {'Penrose':np.array([0,0]),'Random':np.array([0.5,0.5])}
    for i in PL_HC:
        if len(i)>4:
            cen = 0
            visit = []
            for j in i:
                visit.append(j[0])
                visit.append(j[1])
                cen = cen + pos[j[0]]
            if np.linalg.norm(cen/len(i) - C[str]) <= thresh:
                cyc = i
                break
    if str == 'Penrose':
        cyc = arrange_cycle_in_order(cyc)
    return cyc,list(set(visit))

def cross_product(A,B,C,D):
    vec_ab = A - B
    vec_bd = B - D
    vec_cb = C - B
    cross_ab_bd = np.cross(vec_ab, vec_bd)
    cross_cb_bd = np.cross(vec_cb, vec_bd)
    return cross_ab_bd, cross_cb_bd

def is_node_in_cycle(node, cycle_edges):
    cycle_nodes = set()
    for u, v in cycle_edges:
        cycle_nodes.add(u)
        cycle_nodes.add(v)
    polygon_nodes = [pos[n] for n in cycle_nodes]
    polygon_path = mpltPath.Path(polygon_nodes)
    point = pos[node]
    is_inside = polygon_path.contains_point(point)
    return is_inside

def initialV(visit,cyc): #AltCs,index):  #
    #cyc = AltCs[index] #remove
    for u,v in cyc:
        for i in g.neighbors(u):
            if i not in visit:
                sv = i
                past = u
                return past,sv
        for i in g.neighbors(v):
            if i not in visit:
                sv = i
                past = v
                return past,sv
    # index -= 1       #remove
    # return initialV(visit,AltCs,index)    #remove

def StagSt(str,thresh):
    cyc,visit = Cen_Plaq(str,thresh)
    gc = g.copy()
    gc.remove_nodes_from(visit)
    AltCs = [cyc]
    while True:
        past,sv = initialV(visit,cyc)  #,AltCs,index=-1)
        visit.append(sv)
        gc.remove_node(sv)
        cyc = []
        while True:
            nn = []
            for i in g.neighbors(sv):
                if i not in visit:
                    nn.append(i)
            if nn==[]:
                break
            if len(nn)==1:
                if sv!=cyc[0][1] and cyc[0][0] in g.neighbors(sv):
                    break
                else:
                    cyc.append((sv,nn[0]))
                    past = sv
                    sv = nn[0]
                    gc.remove_node(sv)
                    visit.append(nn[0])
            else:
                cpA,cpC = cross_product(pos[nn[0]],pos[sv],pos[nn[1]],pos[past])
                if cpA < 0 and cpC > 0:
                    pick = 0
                elif cpA > 0 and cpC < 0:
                    pick = 1
                elif cpA < 0 and cpC < 0:
                    if cpA < cpC:
                        pick = 0
                    else:
                        pick = 1
                else:
                    if cpA < cpC:
                        pick = 0
                    else:
                        pick = 1
                cyc.append((sv,nn[pick]))
                past = sv
                sv = nn[pick]
                gc.remove_node(sv)
                visit.append(nn[pick])
        cyc.append((cyc[-1][1],cyc[0][0]))
        if g.degree(sv) == 2 or (str=='Random' and is_node_in_cycle(AltCs[0][2][0],cyc)==False):
            print(AltCs[0][3][0])
            break
        else:
            AltCs.append(cyc)
    return AltCs

AltCs = StagSt('Penrose',0.5)    #'Random',0.08
print(len(AltCs))

StagI = []
StagI.extend(AltCs[0][::2])
Stag = StagI
for i in AltCs[1:]:
    c = set()
    for j in i:
        c.add(j[0])
        c.add(j[1])
    for u,v in StagI:
        for k in g.neighbors(v):
            if k in list(c):
                print(u, v)
                break
        break
    for j in range(len(i)):
        if i[j][0] == k:
            break
    if j%2==0:
        StagI = []
        for p in i[::2]:
            StagI.append(p)
    else:
        StagI = []
        for p in i[1::2]:
            StagI.append(p)
    Stag.extend(StagI)
# pickle.dump(Stag,open('../Mij_correF/Stag_RandomG_L=80.txt','wb'))

nodes = set()
for u,v in Stag:
    nodes.add(u)
    nodes.add(v)
rn = []
for i in g.nodes():
    if i not in nodes:
        rn.append(i)
g.remove_nodes_from(rn)
# pickle.dump(g,open('../Mij_correF/Gstag_RandomG_L=80.txt','wb'))

E = []
for i in AltCs:
    E.extend(i)
EC = []
W = []
for u,v in g.edges():
    if (u,v) in Stag or (v,u) in Stag:
        EC.append('green')
        W.append(2)
    # elif (u,v) in E or (v,u) in E:
    #     EC.append('red')
    #     W.append(1.5)
    else:
        EC.append('black')
        W.append(0.7)
nx.draw(g,pos,node_size = 0,edge_color = EC,width=W) #,with_labels=True)
plt.show()


'''
Stag = []
sv = AltCs[0][0][0]
for cyc in AltCs:
    print(sv,cyc)
    for i in range(len(cyc)):
        if cyc[i][0] == sv:
            break
    if i%2==0:
        for j in cyc[::2]:
            Stag.append(j)
            x = cyc[0][0]      #== columnar ??
    else:
        for j in cyc[1::2]:
            Stag.append(j)
            x = cyc[1][0]      #== columnar ??
    #x = cyc[-1][0]            #Stag[-1][0] == columnar ??
    for k in g.neighbors(x):
        if (k,x) not in cyc and (x,k) not in cyc:
            sv=k
            break
'''