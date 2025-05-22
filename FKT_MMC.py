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

g = pickle.load(open('../3DegLattices/(4,6)_OBC.txt','rb'))
pos = pickle.load(open('../3DegLattices/Pos_(4,6)_OBC.txt','rb'))
pos = {u:np.array(v) for u,v in pos.items()}
PL_HC = pickle.load(open('../3DegLattices/Plaqs_(4,6)_OBC.txt','rb'))
print(len(PL_HC))

# g = pickle.load(open('../3DegLattices/SqOct/SqOct12_OBC.txt','rb'))
# pos = pickle.load(open('../3DegLattices/SqOct/Pos_SqOct12_OBC.txt','rb'))
# pos = {u:np.array(v) for u,v in pos.items()}
# PL_HC = pickle.load(open('../3DegLattices/SqOct/Plaqs_SqOct12_OBC.txt','rb'))
#
# g = pickle.load(open('../3DegLattices/Hexa/Hexa30_OBC.txt','rb'))
# pos = pickle.load(open('../3DegLattices/Hexa/Pos_Hexa30_OBC.txt','rb'))
# pos = {u:np.array(v) for u,v in pos.items()}
# PL_HC = pickle.load(open('../3DegLattices/Hexa/Plaqs_Hexa30_OBC.txt','rb'))

# g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))

# g = pickle.load(open('../all_graphs_plantri/G_out_20_15k.txt','rb'))
# pos = pickle.load(open('../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
# PL_HC = pickle.load(open('../all_graphs_plantri/PL_out_20_15k.txt','rb'))

# g = pickle.load(open('../RandomG/G_L=200_tutte.txt','rb'))
# pos = pickle.load(open('../RandomG/Pos_L=200_tutte.txt','rb'))
# PL_HC = pickle.load(open('../RandomG/Plaqs_L=200_tutte.txt','rb'))
# nx.draw(g,pos,node_size=5,with_labels=True)
# plt.show()
# dgf

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
    past = 7                # Plantri_20_15K: (past,sv) = (27,20) ; Penrose_10x10: (past,sv) = (7,3) ; RandomG_L=80: (past,sv) = (261,260)
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
            if False:   #cyc[0] in g.neighbors(sv):   #-Plantri    #False-Penrose
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

def distV():
    DistV = {}
    for i in list(g.nodes())[::30]:                       #[:]-Plantri ; [:30]-Penrose ; [:17]-RandomG
        for j in list(g.nodes())[::30]:
            d = round(np.linalg.norm(pos[i]-pos[j]),1)    #6-Plantri ; 1-Penrose ; 2-RandomG
            if i!=j and d not in DistV.keys() and i not in BDL and j not in BDL:
                DistV[d] = [(i,j)]
            elif d in DistV.keys() and (j,i) not in DistV[d] and i not in BDL and j not in BDL:
                DistV[d].append((i,j))
    return DistV

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

def updatePL(l,u):
    Eu = g.edges(u)
    lc = l.copy()
    rem = []
    for i in lc:
        for a,b in i:
            if (a,b) in Eu or (b,a) in Eu:
                rem.extend(i)
                l.remove(i)
                break
    pl = []
    for i,j in rem:
        if (i,j) not in Eu and (j,i) not in Eu:
            pl.append((i,j))
    l.append(pl)
    return l

def fktcolormap(st):
    nc = {}
    gc = g.copy()
    PLc = updatePL(PL_HC, st)
    PLc_st = PLc.copy()
    for i in g.nodes():
        print(i)
        gc.remove_node(st)
        if i!=st and i not in BDL:
            gc.remove_node(i)
            PLc = updatePL(PLc, i)
            d_pos, PL = dualP(PLc)
            nc[i] = FKT(gc, pos, d_pos, PL, BDL)
        else:
            nc[i] = 0
        gc = g.copy()
        PLc = PLc_st.copy()
    return nc

# BDL = BD_Plaq()
BDL = pickle.load(open('../3DegLattices/BDL_(4,6)_OBC.txt','rb'))

#[850,851,864,865,866,867,297,296,255,254,252,253,260,261,280,281,313,312,235,234,193,192,196,197,383,382,384,385,465,464,466,467,881,880,877,876,45,44,39,38,35,34,5,4,2,3,10,11,25,24,22,23,32,33,17,16,9,8,1,0,6,7,36,37,40,41,42,43,47,46,50,51,60,61,693,692,683,682,679,678,673,672,669,668,670,671,821,820,655,654,645,644,646,647,799,798,794,795,800,801,856,857,853,852]
#[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 1020, 1079, 1080, 2639, 2640, 1619, 1620, 599, 600, 2039, 3179, 3180, 2159, 2160, 3584, 1139, 1140, 3585, 119, 120, 2040, 3586, 3587, 3588, 3589, 2699, 2700, 3590, 1679, 1680, 659, 3591, 660, 539, 3592, 540, 3593, 3059, 3594, 1559, 3239, 3595, 1560, 3240, 2219, 2220, 3596, 1199, 1200, 2579, 3597, 179, 180, 1019, 2580, 3598, 3060, 3599, 2759, 2760, 1739, 1740, 719, 720, 3299, 3300, 2279, 2280, 1259, 1260, 239, 240, 2819, 2820, 1799, 1800, 779, 780, 3359, 3360, 2339, 2340, 1319, 1320, 299, 300, 2879, 2880, 2099, 1859, 1860, 2100, 839, 840, 3119, 3120, 3419, 3420, 2399, 2400, 1379, 1380, 359, 360, 479, 480, 1499, 1500, 2939, 2940, 1919, 1920, 2519, 899, 900, 2520, 3479, 3480, 2459, 2460, 1439, 1440, 419, 420, 2999, 3000, 1979, 1980, 959, 960, 3539, 3540, 3541, 3542, 3543, 3544, 3545, 3546, 3547, 3548, 3549, 3550, 3551, 3552, 3553, 3554, 3555, 3556, 3557, 3558, 3559, 3560, 3561, 3562, 3563, 3564, 3565, 3566, 3567, 3568, 3569, 3570, 3571, 3572, 3573, 3574, 3575, 3576, 3577, 3578, 3579, 3580, 3581, 3582, 3583]
#[261, 260, 46, 47, 125, 124, 109, 108, 293, 292, 105, 104, 146, 147, 101, 100, 165, 164, 337, 336, 161, 160, 333, 332, 25, 24, 77, 76]   #??

DistV = distV()
'''
# st = 79   #520(sqoct12)   #1524(sqoct20)  #216(sqoct8)
# nc = fktcolormap(st)
# print(nc)
# pickle.dump([st,nc],open('../3DegLattices/FKTcmap_Mr_hexa8_obc.txt','wb'))  #../CorreFns/Data_RandG/FKTcmap_Mr_L=200.txt','wb'))

st, nc = pickle.load(open('../3DegLattices/SqOct/FKTcmap_Mr_sqoct12_obc.txt','rb'))
print(st,nc)
NS = []
for u in g.nodes():
    if nc[u]!=0:
        nc[u] = np.log(nc[u])   #nc[u]   #
        NS.append(20)
    elif u==st:
        del nc[u]
        NS.append(20)
    else:
        del nc[u]
        NS.append(0)

norm = Normalize(vmin=min(nc.values()), vmax=max(nc.values()))
cmap = plt.cm.viridis
color_map = [cmap(norm(nc[node])) if (node in nc.keys()) else 'red' for node in g.nodes()]   #node!=st and node not in BDL and
nx.draw(g, pos, node_size=NS, node_color=color_map)  #, with_labels=True, font_size=8)
sm = ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # Only needed for colorbar
plt.colorbar(sm, shrink=0.7,pad=-0.04)  #label="FKT output")
plt.show()


gc = g.copy()
PLc = PL_HC.copy()
M = {}
print(DistV.keys())
for r in DistV.keys():
    print(r,len(DistV[r]))
    s = 0
    for u,v in DistV[r][::5]:    #:5-Penrose ; [:]-Plantri,RandG
        gc.remove_node(u)
        PLc = updatePL(PLc,u)
        gc.remove_node(v)
        PLc = updatePL(PLc,v)
        d_pos,PL = dualP(PLc)
        s += FKT(gc,pos,d_pos,PL,BDL)
        gc = g.copy()
        PLc = PL_HC.copy()
    print(s)
    M[r] = s/(len(DistV)/5)        #/5-Penrose

print(M)
filter_M = {k: v for k, v in M.items() if v != 0}
sort_M = dict(sorted(filter_M.items()))
pickle.dump(sort_M,open('../3DegLattices/Hexa/FKT_Mr_hexa30_obc_1.txt','wb'))   #../CorreFns/Data_RandG/FKT_Mr_L=200.txt','wb'))
'''
import matplotlib.pyplot as plt

# Data
x_values = list(range(1, 9))
y_values = [1,0.5,0.33,0.25,0.2,0.17,0.14,0.125]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, y_values, marker='o', linestyle='-', color='b')
plt.title("Plot of Given Data")
plt.yscale('log')
plt.xlabel("X Values")
plt.ylabel("Y Values")
plt.grid(True)
plt.show()

tyt

sort_M = pickle.load(open('../3DegLattices/Hexa/FKT_Mr_hexa30_obc_1.txt','rb'))   #../CorreFns/Data_RandG/FKT_Mr_L=200.txt','rb'))
print(sort_M)

x = np.array(list(sort_M.keys()))
y = np.array(list(sort_M.values()))

plt.scatter(x,y,s=5)
plt.plot(x,y)
# plt.xscale('log')
# plt.yscale('log')

plt.show()

# def model_function(x, a, b):
#     return b * np.exp(-a * x)
# popt, pcov = curve_fit(model_function, x, y)
# a_f, b_f = popt
# y_fitted = model_function(x, a_f, b_f)
# plt.plot(x, y_fitted, label='Fitted Curve', color='red')
# plt.show()