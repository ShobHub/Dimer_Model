from matplotlib.patches import Polygon
import colorsys
import math
from matplotlib.pyplot import pause
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from collections import deque
import matplotlib.pyplot as plt   #import matplotlib.pylab as plt
import itertools as it
import networkx as nx
import pandas as pd
import numpy as np
import random
import pickle
import pylab

class Heights_Penrose():
    def __init__(self):
        self.PL = PL
        self.PL1 = PL1
        self.M = M
        self.BN,_ = nx.bipartite.sets(G)
    def center(self,pl):
        cx = 0
        cy = 0
        for i in pl:
            cx += pos[i[0]][0]
            cy += pos[i[0]][1]
            cx += pos[i[1]][0]
            cy += pos[i[1]][1]
        return (cx/(2*len(pl)),cy/(2*len(pl)))
    def Plqs(self,node):
        plq = []
        for i in self.PL:
            if node in set(it.chain.from_iterable(i)):
                plq.append(i)
        return plq
    def Dir(self,c,plqs):
        c1 = self.center(plqs[0])
        c2 = self.center(plqs[1])
        a = np.array(c1)-np.array(c)
        b = np.array(c2)-np.array(c)
        if np.cross(a,b)<0:
            return plqs,c1,c2
        else:
            return plqs[::-1],c2,c1
    def dm(self,p1,p2):
        set1 = {frozenset(item) for item in p1}
        set2 = {frozenset(item) for item in p2}
        d1,d2 = tuple(list(set1.intersection(set2))[0])
        if (d1,d2) in self.M or (d2,d1) in self.M:
            return True
        else:
            return False
    def AssignH(self,H,l,C):
        for i in range(len(l)-1):
            flag = self.dm(l[i],l[i+1])
            if flag == True:
                H[C[i+1]] = H[C[i]]-2
            else:
                H[C[i+1]] = H[C[i]]+1
        return H
    def plot(self,H): #,CE):
        #pylab.ion()
        #plt.clf()
        # plt.style.use('dark_background')
        fig, ax = plt.subplots()

        X, Y = zip(*list(H.keys()))
        # ax.scatter(X, Y, s=1, c='gold')
        for i, label in enumerate(list(H.values())):
            ax.annotate(np.round(label,1), list(H.keys())[i], ha='center',c='gold')

        EC = []
        W = []
        for u, v in G.edges():
            if ((u,v) in self.M or (v,u) in self.M):   #frozenset([u, v]) in CE:     #
                EC.append('purple')
                W.append(3)
            else:
                EC.append('white')
                W.append(0.8)
        nx.draw(G, pos, node_size=0, edge_color=EC, width=W) #,with_labels = True,font_size=8)
        #pause(5)
        #pylab.show()
        ax.set_facecolor('black')
        ax.axis('off')
        fig.set_facecolor('black')  # Set the figure background color
        #plt.savefig('RG_H.pdf')
        plt.show()

    def Heights(self):
        H = {}
        pl = self.PL1[random.choice(range(len(self.PL1)))]
        self.PL1.remove(pl)
        c = self.center(pl)
        H[c] = 0
        dq = []
        while len(self.PL1) > 0:
            # print(len(self.PL1))
            self.PL.remove(pl)
            for i in set(it.chain.from_iterable(pl)):
                if i in self.BN:
                    node = i
                    plqs = self.Plqs(node)
                    if len(plqs) == 0:
                        continue
                    elif len(plqs) == 1:
                        continue
                    else:
                        Dplq,c1,c2 = self.Dir(c,plqs)
                        H = self.AssignH(H, [pl,Dplq[0],Dplq[1]], [c,c1,c2])
                        if plqs[0] in self.PL1:
                            dq.append(plqs[0])
                            self.PL1.remove(plqs[0])
                        if plqs[1] in self.PL1:
                            dq.append(plqs[1])
                            self.PL1.remove(plqs[1])
            self.PL.append(pl)
            pl = dq[random.choice(range(len(dq)))]
            c = self.center(pl)
            #self.plot(H)
        return H

    def plotColorH(self, H, vp, t):
        hs=list(set(H.values()))
        coms=sorted(hs)
        # hs = list(range(-16, 66))  # (-2,2)
        # coms = hs
        mycols = {}
        value_range = max(hs) - min(hs)
        for k in range(len(coms)):
            # Normalize the value between 0.2 (dark) and 1.0 (bright)
            normalized_value = 0.1 + (coms[k] - min(hs)) / value_range * 0.9
            rgb_color = colorsys.hsv_to_rgb(39/360,0.9,normalized_value)   #1, 0.5, normalized_value)
            mycols[coms[k]] = rgb_color

        fig, ax = plt.subplots()
        color_list = [mycols[val] for val in coms]
        norm = mcolors.Normalize(vmin=min(hs), vmax=max(hs))
        cmap = cm.ScalarMappable(norm=norm, cmap=plt.cm.colors.ListedColormap(color_list))

        ax = plt.gca()
        for i in self.PL:
            xc = 0
            yc = 0
            pp = []
            for u, v in i:
                xc = xc + pos[u][0]
                yc = yc + pos[u][1]
                xc = xc + pos[v][0]
                yc = yc + pos[v][1]
                pp.extend([u, v])
            x, y = (xc / (2 * len(i)), yc / (2 * len(i)))
            pp = set(pp)
            pp1 = []
            for j in pp:
                pp1.append(pos[j])
            pp1.sort(key=lambda a: math.atan2(a[1] - y, a[0] - x))
            if ((x, y) in H.keys()):
                ax.add_patch(Polygon(pp1, closed=True, fill=True, facecolor=mycols[H[(x, y)]], edgecolor='white',linewidth=0.2))
            else:
                ax.add_patch(Polygon(pp1, closed=True, fill=False))

        cbar = fig.colorbar(cmap, ax=ax, shrink=0.6, c='white')
        cbar.ax.tick_params(axis="y", color='white', labelcolor='white')
        cbar.set_label('Height Values',c='white')
        ax.autoscale_view()
        plt.gca().set_aspect('equal')
        plt.title('V=%.2f ; T=%.2f' % (vp,t), c='white')  #(k,R)=(10,10)
        plt.axis('off')
        ax.set_facecolor('black')
        ax.axis('off')
        fig.set_facecolor('black')
        #plt.savefig('RG_H1.pdf')
        plt.show()

G = pickle.load(open('../RandomG/G_L=80_tutte.txt','rb'))   #ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))         #../3DegLattices/SqOct/SqOct12_OBC.txt','rb'))              #../all_graphs_plantri/G_out_20_15k.txt','rb'))
pos = pickle.load(open('../RandomG/Pos_L=80_tutte.txt','rb'))   #ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))   #../3DegLattices/SqOct/Pos_SqOct12_OBC.txt','rb'))          #../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
# nx.draw(G,pos,node_size = 5)
# plt.show()

V = -10
T,Cv = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (-10),'rb'))

# M = nx.bipartite.hopcroft_karp_matching(G)
# ME = []
# for u,v in M.items():
# 	if((u,v) not in ME and (v,u) not in ME):
# 		ME.append((u,v))

for i in T[11::5]:
    Mlist = pickle.load(open('../CorreFns/Data_RandG/Samples/Mlist_L=80_T=%.2f_V=%.2f.txt'%(i,V),'rb'))   #[Stag_RandomG_L=80.txt','rb'))]   #Data_ModP/Samples/Mlist_(k,R)=(10,10)_T=%.2f_V=%.2f.txt'%(i,V),'rb'))   #../3DegLattices/SqOct/FS/Mlist_sqoct12_T=%.2f_V=%.2f.txt'%(i,V),'rb'))
    print(len(Mlist))
    H_ = []
    for j in Mlist[:1]:  #2]:
        M = j.copy()
        PL = pickle.load(open('../RandomG/Plaqs_L=80_tutte.txt','rb'))    #../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt', 'rb'))    #../3DegLattices/SqOct/Plaqs_SqOct12_OBC.txt','rb'))       #../all_graphs_plantri/PL_out_20_15k.txt','rb'))
        PL1 = PL.copy()   #for random PL1==PL ; [] for Penrose PL1 = PL>4
        # for j in PL:
        #     if len(j) > 4:
        #         PL1.append(j)
        obj = Heights_Penrose()
        H = obj.Heights()
        H_.append(H)
        print(len(PL)-len(H.keys()))
        # plt.scatter([-22.343434343434343],[-28.0],c='red',s=2)
        # pickle.dump(H, open('../CorreFns/Data_ModP/Samples/StagH_(k,R)=(10,10).txt', 'wb'))
        obj.plot(H)
    # df = pd.DataFrame(H_)
    # var_H = df.var(ddof=0).to_dict()
    # obj.plot(var_H)
    obj.plotColorH(H,V,i)
    # pickle.dump(H_, open('../CorreFns/Data_ModP/Samples/H_(k,R)=(10,10)_T=%.2f_V=%.2f.txt'%(i,V), 'wb'))
    # pickle.dump(var_H,open('../3DegLattices/SqOct/Heights/VarH_sqoct12_T=%.2f_V=%.2f.txt'%(i,V),'wb'))

fdhg

print(len(PL))
cent = np.array((0,0))   #0.5,0.5))-Plantri/RG  (0,0)-Penrose
DistC = {}
for c in H.keys():
    #c = obj.center(i)
    d = np.round(np.linalg.norm(np.array(c) - cent),2)
    if d not in DistC.keys():
        DistC[d] = [c]
    else:
        DistC[d].append(c)
print(len(DistC.keys()),max(DistC),DistC)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(1, 1, figsize=(79, 74)) #, sharex=True, sharey=True)
norm = mcolors.Normalize(vmin=min(T[10:-30:2]), vmax=max(T[10:-30:2]))   #[25:-25:2]  [0:-30:3]
cmap = cm.get_cmap('rainbow')

t = T[0:-30:3][2:]
c = 0
L = np.linspace(0.1,35,60)     #(0.1,44,60)-Penrose  #(0.15,0.5,12)-Random  #(0.1,0.38,5)-Plantri  (0.1,35,60)-SqOct
# F = np.log(L/90)
for i in range(0,len(T)-30,3)[2:]:  #[25:-25:2]
    var_H = pickle.load(open('../3DegLattices/SqOct/Heights/VarH_sqoct12_T=%.2f_V=%.2f.txt'%(T[i],V),'rb'))   #../CorreFns/Data_RandG/Samples/VarH_L=80_T=%.2f_V=%.2f.txt' % (T[i], V), 'rb'))
    # obj.plotColorH(var_H, V, i)
    W = []
    for j in L:
        W.append(np.mean([var_H[q] for k, c in DistC.items() if k <= j for q in c]))
    plt.scatter(L/40,W,s=3000,label='T=%.2f'%(T[i]),color=cmap(norm(T[i])))   #L/90
    plt.plot(L/40,W,color=cmap(norm(T[i])),linewidth=12)
    # plt.plot(L,F,c='r')
    c += 1

plt.xscale('log')
plt.xlabel(r"System Size ($L$)",fontsize = 270)
plt.ylabel(r"Height Variance ($W^2$)",fontsize = 270)
axs.tick_params(axis='x', labelsize=230, length=60, width=15)  # Increase x-axis numbers size
plt.tick_params(axis='y', labelsize=230, length=60, width=15)
plt.legend(fontsize = 210,loc='lower right',ncol=2)
plt.title(r"Sq-Oct ($v = -1$)",fontsize = 270)
plt.xlim(left = 0)  #0.14,right=0.55)  #0.04)
plt.ylim(bottom = -2.5)
plt.savefig('../../papers/RandomFieldTheory/varH_sqoct_thesis.pdf')
# plt.show()












