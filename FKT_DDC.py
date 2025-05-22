from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import networkx as nx
import numpy as np
import pickle
import sys
import os
sys.path.append('../../HatTile/Spectre')
from FKT import FKT

# g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))

# g = pickle.load(open('../all_graphs_plantri/G_out_20_15k.txt','rb'))
# pos = pickle.load(open('../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
# PL_HC = pickle.load(open('../all_graphs_plantri/PL_out_20_15k.txt','rb'))

g = pickle.load(open('../RandomG/G_L=80_tutte.txt','rb'))
pos = pickle.load(open('../RandomG/Pos_L=80_tutte.txt','rb'))
PL_HC = pickle.load(open('../RandomG/Plaqs_L=80_tutte.txt','rb'))
print(len(g.edges()))

def center(pl):
    cen = 0
    for u,v in pl:
        cen += pos[u]
    return cen/len(pl)

def cross_product(A,B,C,D):
    vec_ab = A - B
    vec_bd = B - D
    vec_cb = C - B
    cross_ab_bd = np.cross(vec_ab, vec_bd)
    cross_cb_bd = np.cross(vec_cb, vec_bd)
    return cross_ab_bd, cross_cb_bd

def BD_Plaq():
    past = 7                 # Plantri_20_15K: (past,sv) = (27,20) ; Penrose_10x10: (past,sv) = (7,3) ; RandomG_L=80: (past,sv) = (261,260)
    sv = 3
    visit = [7,3]
    cyc = [7,3]
    while True:
        nn = []
        for i in g.neighbors(sv):
            if i not in visit:
                nn.append(i)
        if nn == []:
            break
        if len(nn) == 1:
            if False:   #cyc[0] in g.neighbors(sv):   #-Plantri   #False-Penrose
                break
            else:
                cyc.append(nn[0])
                past = sv
                sv = nn[0]
                visit.append(nn[0])
        else:
            cpA, cpC = cross_product(pos[nn[0]], pos[sv], pos[nn[1]], pos[past])
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
            cyc.append(nn[pick])
            past = sv
            sv = nn[pick]
            visit.append(nn[pick])
    #cyc.append((sv,cyc[0][0]))
    return cyc

def distE():
    DistE = {}
    for i in list(g.edges())[::30]:                           #[:]-Plantri ; [:30]-Penrose ; [:17]-RandomG
        cen1 = (pos[i[0]] + pos[i[1]]) / 2
        for j in list(g.edges())[::30]:
            if i[0] not in BDL and i[1] not in BDL and j[0] not in BDL and j[1] not in BDL and not set(i) & set(j):
                cen2 = (pos[j[0]] + pos[j[1]]) / 2
                d = round(np.linalg.norm(cen1 - cen2),1)          #2-Plantri ; 1-Penrose ; 2-RandomG
                if d != 0 and d not in DistE.keys():
                    DistE[d] = [(i, j)]
                elif d in DistE.keys() and (j, i) not in DistE[d]:
                    DistE[d].append((i, j))
    return DistE

def dualP(l):
    d_pos = {}
    PL = {}
    c = 0
    for i in l:
        cen = center(i)
        d_pos[c] = cen
        PL[c] = i
        c += 1
    return d_pos,PL

def updatePL(l,a):
    lc = l.copy()
    rem = []
    for i in lc:
        for tup in i:
            if frozenset(tup) == frozenset(a):
                rem.extend(i)
                l.remove(i)
                break
    pl = []
    for tup in rem:
        if frozenset(tup) != frozenset(a):
            pl.append(tup)
    l.append(pl)
    return l

def fktcolormap(se):
    ec = {}
    sed = list(g.edges())[se]
    ec[sed] = 0
    PLo = PL_HC.copy()
    g.remove_edge(*sed)
    gc = g.copy()
    PLc = updatePL(PL_HC, sed)
    PLc_st = PLc.copy()
    d_pos, PL = dualP(PLc)
    d_i = FKT(gc, pos, d_pos, PL, BDL)
    for i in g.edges():
        if (i[0] not in BDL and i[1] not in BDL) and (i[0] not in sed and i[1] not in sed):
            print(i)
            gc.remove_edge(*i)
            PLc = updatePL(PLc, i)
            d_pos, PL = dualP(PLc)
            d_ij = FKT(gc, pos, d_pos, PL, BDL)

            gc = g.copy()
            gc.add_edge(*sed)
            PLc = PLo.copy()

            gc.remove_edge(*i)
            PLc = updatePL(PLc, i)
            d_pos, PL = dualP(PLc)
            d_j = FKT(gc, pos, d_pos, PL, BDL)

            ec[i] = d_ij - d_i*d_j
        else:
            ec[i] = 0
        gc = g.copy()
        PLc = PLc_st.copy()
    return ec

# BDL = BD_Plaq()
# BDL = [261, 260, 46, 47, 125, 124, 109, 108, 293, 292, 105, 104, 146, 147, 101, 100, 165, 164, 337, 336, 161, 160, 333, 332, 25, 24, 77, 76]   #RG_L=80
# DistV = distE()

se = 218
sed = list(g.edges())[se]
# ec = fktcolormap(se)
# print(ec)
# pickle.dump([se,ec],open('../CorreFns/Data_RandG/FKTcmap_Dr_L=80_se=%i.txt'%se,'wb'))
se, ec = pickle.load(open('../CorreFns/Data_RandG/FKTcmap_Dr_L=80_se=%i.txt'%se,'rb'))

for u in g.edges():
    if ec[u]!=0:
        ec[u] = ec[u]/7.43350043902808e+29    #7.43350043902808e+29 (RG)    2.095608431208898e+143 (Penrose)  #np.sign(ec[u])*np.log10(abs(ec[u]))
    else:
        del ec[u]

norm = Normalize(vmin=min(ec.values()), vmax=max(ec.values()))
cmap = plt.cm.coolwarm
color_map = [cmap(norm(ec[edge])) if (edge in ec.keys()) else 'black' for edge in g.edges()]    #edge!=sed and (edge[0] not in BDL and edge[1] not in BDL)
W = [2 if (edge in ec.keys()) else 1 for edge in g.edges()]
nx.draw(g, pos, node_size=0, edge_color=color_map,width=W)  #, with_labels=True, font_size=8)

sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for colorbar
plt.colorbar(sm, label="Exact Dimer-Dimer Correlations (FKT)",shrink=0.6)
plt.savefig('../../../Downloads/RG_DD_1.pdf')
plt.show()

gtr

x = []
y = []
for u,v in ec.items():  #nc.items():
    cen1 = (pos[u[0]] + pos[u[1]]) / 2
    cen2 = (pos[sed[0]] + pos[sed[1]]) / 2
    x.append(np.linalg.norm(cen1-cen2))
    y.append(v)
print(x)
print(y)
plt.scatter(x,y,s=5,label = 'fixed_edge = (%i,%i)'%sed)
#plt.yscale('log')
plt.ylabel('FKT output')
plt.xlabel('r')
plt.legend()
plt.show()

dfg

# gc = g.copy()
# PLc = PL_HC.copy()
# D = {}
# print(DistE.keys())
for r in []: #list(DistE.keys())[::5]:         #[::5]-Penrose ; [:]-Plantri,RandG
    print(r,len(DistE[r]))
    si = []
    sj = []
    sij = []
    for i in DistE[r][:]:             #[:]-Penrose,Plantri,RandG
        gc.remove_edge(*i[0])
        PLc = updatePL(PLc,i[0])
        d_pos, PL = dualP(PLc)
        si.append(FKT(gc,pos,d_pos,PL,BDL))

        gc.remove_edge(*i[1])
        PLc = updatePL(PLc,i[1])
        d_pos,PL = dualP(PLc)
        sij.append(FKT(gc,pos,d_pos,PL,BDL))

        gc = g.copy()
        PLc = PL_HC.copy()

        gc.remove_edge(*i[1])
        PLc = updatePL(PLc, i[1])
        d_pos, PL = dualP(PLc)
        sj.append(FKT(gc, pos, d_pos, PL, BDL))

        gc = g.copy()
        PLc = PL_HC.copy()
    D[r] = (np.mean(sij) - np.mean(si)*np.mean(sj))

# print(D)
#
# filter_D = {k: v for k, v in D.items() if v != 0}
# sort_D = dict(sorted(filter_D.items()))
# print(sort_D)

#pickle.dump(sort_D,open('../CorreFns/Data_ModP/FKT_Dr_(k,R)=(10,10).txt','wb'))
sort_D = pickle.load(open('../CorreFns/Data_RandG/FKT_Dr_L=80.txt','rb'))
sort_D = {k: v*(-1**k) for k, v in sort_D.items()}

x = np.array(list(sort_D.keys()))
y = np.array(list(sort_D.values()))

plt.scatter(sort_D.keys(),sort_D.values(),s=5)
plt.plot(sort_D.keys(),sort_D.values())
# plt.xscale('log')
# plt.yscale('log')

# def model_function(x, a, b):
#     return a*x+b
# popt, pcov = curve_fit(model_function, np.log(x)[4:-17], np.log(y)[4:-17])
# a_f, b_f = popt
# print(popt)
# y_fitted = model_function(np.log(x)[4:-17], a_f, b_f)
# plt.plot(x[4:-17], np.exp(y_fitted), label='Fitted Curve', color='red')

plt.show()