from NX_to_koala import nx_to_koala
from koala_to_NX import koala_to_nx
from koala import plotting
from koala.graph_utils import clockwise_about
from koala import graph_utils as gu
from koala.functions.koala_plantri import tutte_embedding
from Dimerisations import count_dimers,plot_dimers
from koala.lattice import Lattice

import networkx as nx
import shapely as shp
import itertools as it
from math import atan2, degrees
from shapely.geometry import LineString,Point
from scipy.spatial import distance
from shapely.ops import polygonize
from itertools import combinations
import matplotlib.pylab as plt
from matplotlib.pyplot import pause
from mpl_toolkits.mplot3d import Axes3D
from skspatial.objects import Line
import numpy as np
import pickle
import random
import collections
import pylab
import math

def Square(le):
    Nodes = [0,1,2,3]
    pos = {0:(-le/2,le/2),1:(le/2,le/2),2:(le/2,-le/2),3:(-le/2,-le/2)}
    BE = [(0,1),(1,2),(2,3),(3,0)]
    return BE,pos

def checkOnBE(x,y,l):
    BE,pos = Square(le)
    if l != []:
        x1,y1 = l[-1]
    for u,v in BE:
        line = Line.from_points(pos[u],pos[v])
        if line.side_point((x,y)) == 0:
            return (u,v),'on'
        elif line.side_point((x,y)) == -1:
            line_a = Line.from_points([x,y],[x1,y1])
            line_b = Line.from_points(pos[u], pos[v])
            x2,y2 = tuple(line_a.intersect_line(line_b))
            if abs(x2)<=le/2 and abs(y2)<=le/2 :
                return (u,v),'out'
            else:
                continue
    return (-1,-1),'in'

def Bpoints(a,x,y):
    BE, pos = Square(le)
    ind = BE.index(a)
    if ind==0:
        return x,y-le
    elif ind==1:
        return x-le,y
    elif ind==2:
        return x,y+le
    elif ind==3:
        return x+le,y

def sips(p1,p2,Pos):
    pts = {}
    line = np.array(p2)-np.array(p1)
    for u,v in Pos.items():
        pt = np.array(v)-np.array(p1)
        pt1 = np.array(v)-np.array(p2)
        cosang = abs(round(np.dot(line,pt)/(np.linalg.norm(line)*np.linalg.norm(pt)),1))
        if abs(cosang-1)<=0.000001 and abs(np.linalg.norm(pt)+np.linalg.norm(pt1)-np.linalg.norm(line))<=0.000001:
            pts[np.linalg.norm(np.array(v)-np.array(p1))] = u
    return pts

def Random_walk(N,le):
    x = random.uniform(-le/2+0.1,le/2-0.1)
    y = random.uniform(-le/2+0.1,le/2-0.1)
    r = 3
    L = []
    l = []
    while N>0:
        a,str = checkOnBE(x,y,l)
        if str == 'in':
            p = 0
            l.append((x,y))
            theta = 2*np.pi*random.uniform(0,1)
            x,y = x + r*np.cos(theta),y + r*np.sin(theta)
        elif str == 'on':
            p = 1
            l.append((x,y))
            theta = 2 * np.pi * random.uniform(0, 1)
            x, y = x + r * np.cos(theta), y + r * np.sin(theta)
        elif str == 'out':
            if p==0:
                line_a = Line.from_points([x,y], [l[-1][0],l[-1][1]])
                line_b = Line.from_points(pos[a[0]],pos[a[1]])
                l.append(tuple(line_a.intersect_line(line_b)))
                L.append(l)
                x,y = Bpoints(a, x, y)
                l = [Bpoints(a, *l[-1])]   #, (x,y)]
                #theta = 2 * np.pi * random.uniform(0, 1)
                #x, y = x + r * np.cos(theta), y + r * np.sin(theta)
            elif p==1:
                L.append(l)
                x,y = Bpoints(a, x, y)
                l = [Bpoints(a,*l[-1])]    #,(x,y)]
                #theta = 2 * np.pi * random.uniform(0, 1)
                #x, y = x + r * np.cos(theta), y + r * np.sin(theta)
        N -= 1
    return L

le = 5
BE,pos= Square(le)
gsq = nx.Graph()
gsq.add_edges_from(BE)
#nx.draw(gsq,pos,node_size=1,width=0.5,font_size=0.5)

L = Random_walk(80,le)
L[-1].append(L[0][0])

# for i in L[:]:
#     for j in range(len(i)-1):
#         #pylab.ion()
#         plt.plot(*zip(i[j],i[j+1]),c=(0,0,1,0.5))
#         plt.scatter(*zip(i[j], i[j + 1]),s=10,c='b')
#         plt.annotate((round(i[j][0],1),round(i[j][1],1)),(i[j][0],i[j][1]))
#         plt.annotate((round(i[j+1][0],1), round(i[j+1][1], 1)), (i[j+1][0], i[j+1][1]))
        # pause(0.5)
        # pylab.show()

mls = shp.MultiLineString(L)
mls = shp.unary_union(mls)
coords_list = []
for non_itersecting_ls in mls.geoms:
    coords_list.extend(non_itersecting_ls.coords)

IntPt = L[0][0]
SIPs = []
Pos = {}
n = 0
for item, count in collections.Counter(coords_list).items():
    if count > 1 and item!=IntPt:
        #plt.scatter(*item,s=20,c='red')
        SIPs.append(item)
        Pos[n] = np.array(item)
        n += 1

Lst = []
for i in L:
    for j in range(len(i)-1):
        pts = sips(i[j], i[j+1], Pos)
        if len(pts.keys())==0:
            continue
        elif len(pts.keys())==1:
            Lst.append([i[j],list(pts.values())[0],i[j+1]])
        else:
            pts = dict(sorted(pts.items()))
            pv = list(pts.values())
            Len = len(pv)
            Lst.append([i[j],pv[0],pv[1]])
            for k in range(1,Len-1):
                Lst.append([pv[k-1],pv[k],pv[k+1]])
            Lst.append([pv[-2], pv[-1], i[j+1]])

while True:
    Ls = Lst.copy()
    c=0
    for i in range(len(Lst)-1):
        if Lst[i][1] == Lst[i+1][1]:
            del Ls[i:i+2]
            c = 1
    if c==1:
        Lst = Ls
    else:
        while Lst[0][1] == Lst[-1][1]:
            del Lst[0]
            del Lst[-1]
        break
Lst.append(Lst[0])

g = nx.Graph()
pos = {}
n = max(Pos.keys())+1
SqN = {}
for i in range(len(Lst)-1):
    if isinstance(Lst[i][2],int):
        vec1 = Pos[Lst[i][2]] - Pos[Lst[i][1]]
    else:
        vec1 = np.array(Lst[i][2]) - Pos[Lst[i][1]]
    if isinstance(Lst[i+1][0],int):
        vec2 = Pos[Lst[i+1][0]] - Pos[Lst[i+1][1]]
    else:
        vec2 = np.array(Lst[i+1][0]) - Pos[Lst[i+1][1]]
    pos[n] = Pos[Lst[i][1]]+0.05*(vec1/np.linalg.norm(vec1))
    if Lst[i][1] in SqN.keys():
        SqN[Lst[i][1]].append(n)
    else:
        SqN[Lst[i][1]] = [n]
    pos[n+1] = Pos[Lst[i+1][1]] + 0.05*(vec2/np.linalg.norm(vec2))
    if Lst[i+1][1] in SqN.keys():
        SqN[Lst[i+1][1]].append(n+1)
    else:
        SqN[Lst[i+1][1]] = [n+1]
    g.add_edge(n,n+1)
    n = n+2
#print(SqN)

pl_sq = []
for i in SqN.keys():
    ord = {}
    pl = []
    for j in SqN[i]:
        angle = degrees(atan2(pos[j][1]-Pos[i][1],pos[j][0]-Pos[i][0]))
        ord[(angle + 360) % 360] = j
    ord_s = sorted(ord)
    for j in range(len(ord_s)-1):
        g.add_edge(ord[ord_s[j]],ord[ord_s[j+1]])
        pl.append((ord[ord_s[j]],ord[ord_s[j+1]]))
    g.add_edge(ord[ord_s[-1]], ord[ord_s[0]])
    pl.append((ord[ord_s[-1]], ord[ord_s[0]]))
    pl_sq.append(pl)
g.add_nodes_from(pos.keys())
#print(len(pl_sq))

pickle.dump(g,open('../RandomG/G_L=80_3.txt','wb'))
pickle.dump(pos,open('../RandomG/Pos_L=80_3.txt','wb'))
pickle.dump(pl_sq,open('../RandomG/PlaqSQ_L=80_3.txt','wb'))
#pickle.dump(EB1,open('../RandomG/PbcE_L=80_3.txt','wb'))

# g = pickle.load(open('../RandomG/G_L=80.txt','rb'))
# pos = pickle.load(open('../RandomG/Pos_L=80.txt','rb'))
# PL = pickle.load(open('../RandomG/Plaqs_L=80.txt','rb'))
# pos = nx.planar_layout(g)

for u,v in g.degree():
    if v!=3:
        print(u,v,'##')

print(nx.is_bipartite(g))
print(nx.is_planar(g))
#nx.draw(g,pos,node_size=5) #,with_labels=True)
print(len(g.nodes),len(g.edges))  #,len(PL))

# nx.draw(g,pos,node_size=5)
# plt.show()

lattice = nx_to_koala(pos,g)
print(lattice,lattice.n_plaquettes)

lat_list = []
for i in range(len(g.nodes)):
    l = clockwise_about(i,lattice)[0]
    lat_list.append(l[::-1])

PosT = tutte_embedding(lat_list)
E = lattice.edges.indices
crossing = np.zeros((len(E),2))
lattice = Lattice(PosT, E, crossing)
print(lattice,lattice.n_plaquettes)

# dimers = gu.dimerise(lattice,500)
#
# PLSC,EXSC,MFI,LFI = count_dimers(dimers,lattice)
# plot_dimers(lattice,dimers,PLSC,EXSC,MFI,LFI)
#
# plotting.plot_edges(lattice,dimers[LFI[2]])
# plt.show()

pos,E,PL = koala_to_nx(lattice)
print(len(PL))

G = nx.Graph()
G.add_edges_from(E)
nx.draw(G,pos,node_size=0) #,with_labels=True)
plt.show()

# pickle.dump(G,open('../RandomG/G_L=150_tutte.txt','wb'))
# pickle.dump(pos,open('../RandomG/Pos_L=150_tutte.txt','wb'))
# pickle.dump(PL,open('../RandomG/Plaqs_L=150_tutte.txt','wb'))





'''
gex = g.copy()
pos_ex = pos.copy()
Nex = len(g.nodes())
c_ex = ['b']*Nex
R_ex = {}
for u,v in EB1.keys():
    pos_ex[Nex+1] = L[EB1[(u,v)][0]][-1]
    pos_ex[Nex+2] = L[EB1[(u,v)][1]][0]
    gex.add_edges_from([(u,Nex+1),(v,Nex+2)])
    gex.remove_edge(u,v)
    gex.add_nodes_from([Nex+1,Nex+2])
    colt = (random.uniform(0,1),random.uniform(0,1),random.uniform(0,1),1)
    c_ex.append(colt)
    c_ex.append(colt)
    R_ex[(u,v)] = [(u,Nex+1),(v,Nex+2)]
    Nex += 2
nx.draw(gex,pos_ex,node_size=10,node_color=c_ex,with_labels=True)
plt.show()

ME = [(u,v) for u,v in nx.bipartite.hopcroft_karp_matching(gex).items()]

# nx.draw_networkx_nodes(g, pos, node_color = 'r', node_size = 10, alpha = 1)
# ax = plt.gca()
# visit = []
# for e in g.edges:
#     u,v = (e[0],e[1])
#     if(((u,v) in ME or (v,u) in ME) and (u,v) not in visit):
#         visit.append((u,v))
#         ax.annotate("",
#                     xy=pos[e[0]], xycoords='data',
#                     xytext=pos[e[1]], textcoords='data',
#                     arrowprops=dict(arrowstyle="-", color="purple",lw=2,
#                                     shrinkA=5, shrinkB=5,
#                                     patchA=None, patchB=None,
#                                     connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
#                                     ),
#                                     ),
#                     )
#     else:
#         ax.annotate("",
#                     xy=pos[e[0]], xycoords='data',
#                     xytext=pos[e[1]], textcoords='data',
#                     arrowprops=dict(arrowstyle="-", color="0",
#                                     shrinkA=5, shrinkB=5,
#                                     patchA=None, patchB=None,
#                                     connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
#                                                                            ),
#                                     ),
#                     )
# plt.axis('off')
# plt.show()

# pos_ex = nx.spring_layout(gex)
# pos = nx.spring_layout(g)
nx.draw_networkx_nodes(gex, pos_ex, node_color = c_ex, node_size = 10, alpha = 1)
ax = plt.gca()
visit = []
for e in g.edges:
    u,v = (e[0],e[1])
    if((u,v) not in R_ex.keys() and (v,u) not in R_ex.keys()):
        if(((u,v) in ME or (v,u) in ME) and (u,v) not in visit):
            visit.append((u,v))
            ax.annotate("",
                        xy=pos_ex[e[0]], xycoords='data',
                        xytext=pos_ex[e[1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="red",lw=2,
                                        shrinkA=2, shrinkB=2,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
                                        ),
                                        ),
                        )
        else:
            ax.annotate("",
                        xy=pos_ex[e[0]], xycoords='data',
                        xytext=pos_ex[e[1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="0",
                                        shrinkA=2, shrinkB=2,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
                                                                               ),
                                        ),
                        )
    else:
        if(((u, v) in ME or (v, u) in ME) and (u,v) not in visit):
            visit.append((u,v))
            if((u,v) in R_ex.keys()):
                a = R_ex[(u,v)]
            else:
                a = R_ex[(v,u)]
            ax.annotate("",
                        xy=pos_ex[a[0][0]], xycoords='data',
                        xytext=pos_ex[a[1][1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="red", lw=2,
                                        shrinkA=2, shrinkB=2,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
                                                                               ),
                                        ),
                        )
            ax.annotate("",
                        xy=pos_ex[a[1][0]], xycoords='data',
                        xytext=pos_ex[a[1][1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="red", lw=2,
                                        shrinkA=2, shrinkB=2,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
                                                                               ),
                                        ),
                        )
        else:
            if ((u, v) in R_ex.keys()):
                a = R_ex[(u, v)]
            else:
                a = R_ex[(v, u)]
            ax.annotate("",
                        xy=pos_ex[a[0][0]], xycoords='data',
                        xytext=pos_ex[a[0][1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="0",
                                        shrinkA=5, shrinkB=5,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
                                                                               ),
                                        ),
                        )
            ax.annotate("",
                        xy=pos_ex[a[1][0]], xycoords='data',
                        xytext=pos_ex[a[1][1]], textcoords='data',
                        arrowprops=dict(arrowstyle="-", color="0",
                                        shrinkA=5, shrinkB=5,
                                        patchA=None, patchB=None,
                                        connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
                                                                               ),
                                        ),
                        )
plt.axis('off')
plt.show()
'''























#------------------------------*****---------------------------------

# def plot_2torus(R,r):
#     angle = np.linspace(0, 2 * np.pi, 32)
#     u,v = np.meshgrid(angle, angle)
#     X = (R + r * np.cos(v)) * np.cos(u)
#     Y = (R + r * np.cos(v)) * np.sin(u)
#     Z = r * np.sin(v)
#     ax.plot_surface(X, Y, Z, color=[0,0,1,0.3])
#
# def random_angles_2torus(u1,u2,v1,v2):
#     u = np.random.uniform(u1,u2)
#     v = np.random.uniform(v1,v2)
#     return u,v
#
# def random_points_2torus(R,r,u,v):
#     x = (R + r * np.cos(v)) * np.cos(u)
#     y = (R + r * np.cos(v)) * np.sin(u)
#     z = r * np.sin(v)
#     # ax.scatter(x, y, z, s=1, c='red')
#     # pause(0.5)
#     return x,y,z
#
# def random_walk(N,R,r):
#     su,sv = random_angles_2torus(0,2*np.pi,0,2*np.pi)
#     SU,SV = su,sv
#     L = [random_points_2torus(R,r,su,sv)]
#     du = 0.1
#     dv = 0.1
#     while N>=0:
#         u,v = random_angles_2torus(su-du,su+du,sv-dv,sv+dv)
#         L.append(random_points_2torus(R,r,u,v))
#         su = u
#         sv = v
#         N = N-1
#     u = np.linspace(su,SU,10)
#     for i in u:
#         v = np.random.uniform(sv,SV)
#         L.append(random_points_2torus(R,r,i,v))
#     L.append(L[0])
#     return L
#
# #ax = plt.axes(projection='3d')
# #plot_2torus(2,1)
# L = random_walk(5,2,1)
# IntPt = L[0]
# L = shp.LineString(L)
#
# coords = L.coords
# x, y, z = zip(*coords)
# L_2d = shp.LineString(zip(x,y))
# #plt.show()
# #ax.plot(x,y,z, marker='o', markersize=2, linestyle='-', color='b')
# #plt.show()
#
# mls = shp.unary_union(L)
# coords_list = []
# for non_itersecting_ls in mls.geoms:
#     coords_list.extend(non_itersecting_ls.coords)
#
# SIPs = []
# Pos = {}
# n = 0
# for item, count in collections.Counter(coords_list).items():
#     if count > 1 and item!=IntPt:
#         #ax.scatter(*item,s=10,c='yellow')
#         SIPs.append(item)
#         Pos[n] = item
#         n += 1
#     elif count > 1 and item==IntPt:
#         print('*')
#
# #plt.show()
#
# G = nx.MultiGraph()
# G.add_nodes_from(Pos.keys())
# E = []
# mls_2d = shp.unary_union(L_2d)
# polygons = list(polygonize(mls_2d))
# for p in polygons:
#     c = 0
#     C = p.exterior.coords
#     I = []
#     for u,v,w in SIPs:
#         if((u,v) in C or (v,u) in C):
#             c += 1
#             I.append((u,v,w))
#     if c==1:
#         a = list(Pos.keys())[list(Pos.values()).index(I[0])]
#         E.append((a,a))
#     elif c==2:
#         a = list(Pos.keys())[list(Pos.values()).index(I[0])]
#         b = list(Pos.keys())[list(Pos.values()).index(I[1])]
#         E.append((a,b))
#         E.append((a,b))
#         # xs,ys,zs = list(zip(Pos[a],Pos[b]))
#         # ax.plot(xs,ys,zs)
#     elif c==3:
#         a = list(Pos.keys())[list(Pos.values()).index(I[0])]
#         b = list(Pos.keys())[list(Pos.values()).index(I[1])]
#         c = list(Pos.keys())[list(Pos.values()).index(I[2])]
#         E.append((a,b))
#         E.append((b,c))
#         E.append((c,a))
#         # xs,ys,zs = list(zip(Pos[a],Pos[b]))
#         # ax.plot(xs,ys,zs)
#
# Pos = {u:(v[0],v[1]) for u,v in Pos.items()}
# G.add_edges_from(E)
# print(G.edges())
# for u,v in G.degree():
#     if v!=4:
#         print(u,v,'#')
#
# EC = []
# W = []
# visit = []
# for u,v in G.edges():
#     if((u,v) not in visit and (v,u) not in visit):
#         EC.append('red')
#         W.append(4)
#         visit.append((u,v))
#     else:
#         EC.append('black')
#         W.append(1)
#         visit.append((u,v))
#
# nx.draw(G,pos=Pos,node_size=10,edge_color=EC,width=W,alpha=0.5) #with_labels=True)
# plt.show()
# # hgf
#
# #SIPs
# #edges 4-regular b/w SIPs using polygons - doesn't consider all cases of polygons
#
# g = nx.MultiGraph()
# pos = {}
# e = []
# n=0
# vec1 = np.array((1/np.sqrt(2),1/np.sqrt(2)))  #,0))
# vec2 = np.array((-1/np.sqrt(2),1/np.sqrt(2)))  #,0))
# Et = {i:[] for i in G.edges() if(i[0]!=i[1])}
# for i in G.nodes():
#     cen = Pos[i]
#     dir = [vec1, -1*vec2, -1*vec1, vec2]
#     e.extend([(n,n+1),(n+1,n+2),(n+2,n+3),(n+3,n)])
#     for j in G.edges(i):
#         if(j[0]==j[1]):
#             pos[n] = cen + 0.01*dir[0]
#             pos[n+1] = cen + 0.01*dir[1]
#             e.append((n,n+1))
#             dir.pop(0)
#             dir.pop(0)
#             n += 2
#         else:
#             pos[n] = cen + 0.01*dir[0]
#             if((j[0],j[1]) in Et.keys()):
#                 Et[j].append(n)
#             elif((j[1],j[0]) in Et.keys()):
#                 Et[(j[1],j[0])].append(n)
#             dir.pop(0)
#             n += 1
#
# for i in Et.keys():
#     l = list(zip(Et[i][:int(len(Et[i])/2)],Et[i][int(len(Et[i])/2):][::-1]))
#     e.extend(l)
#
# pos = {u:(v[0],v[1]) for u,v in pos.items()}
# g.add_nodes_from(pos.keys())
# g.add_edges_from(e)
# print(g.edges())
# for u,v in g.degree():
#     if v!=3:
#         print(u,v,'##')
# # nx.draw(g,pos,node_size=5,with_labels=True)
# # plt.show()
#
# ME = [(u,v) for u,v in nx.bipartite.hopcroft_karp_matching(g).items()]
# print(ME)
#
# nx.draw_networkx_nodes(g, pos, node_color = 'r', node_size = 10, alpha = 1)
# ax = plt.gca()
# visit = []
# for e in g.edges:
#     u,v = (e[0],e[1])
#     if(((u,v) in ME or (v,u) in ME) and (u,v) not in visit):
#         visit.append((u,v))
#         ax.annotate("",
#                     xy=pos[e[0]], xycoords='data',
#                     xytext=pos[e[1]], textcoords='data',
#                     arrowprops=dict(arrowstyle="-", color="purple",lw=2,
#                                     shrinkA=5, shrinkB=5,
#                                     patchA=None, patchB=None,
#                                     connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
#                                     ),
#                                     ),
#                     )
#     else:
#         ax.annotate("",
#                     xy=pos[e[0]], xycoords='data',
#                     xytext=pos[e[1]], textcoords='data',
#                     arrowprops=dict(arrowstyle="-", color="0",
#                                     shrinkA=5, shrinkB=5,
#                                     patchA=None, patchB=None,
#                                     connectionstyle="arc3,rad=rrr".replace('rrr', str(0.3 * e[2])
#                                                                            ),
#                                     ),
#                     )
# plt.axis('off')
# plt.show()

#------------------------------*****---------------------------------

# def generate_random_line():
#     x = np.random.rand(2)
#     y = np.random.rand(2)
#     return LineString(np.column_stack((x, y)))
#
# def plot_lines(lines):
#     for line in lines:
#         plt.plot(line.xy[0], line.xy[1], c='gray')
#     plt.show()
#
# def find_intersection_points(lines):
#     intersections = np.zeros((1,6))
#     line_combinations = combinations(lines, 2)
#     for line1, line2 in line_combinations:
#         l1c = np.array(line1.coords)
#         l2c = np.array(line2.coords)
#         uL1 = (l1c[0]-l1c[1])/np.linalg.norm(l1c[0]-l1c[1])
#         uL2 = (l2c[0]-l2c[1])/np.linalg.norm(l2c[0]-l2c[1])
#         if line1.intersects(line2):
#             intersection = line1.intersection(line2)
#             if intersection.is_empty:
#                 print('**')
#                 continue
#             if intersection.geom_type == 'Point':
#                 intersections = np.vstack([intersections, [intersection.x, intersection.y,uL1[0],uL1[1],uL2[0],uL2[1]]])
#             elif intersection.geom_type == 'MultiPoint':
#                 print('*')
#                 intersections.extend([(point.x, point.y) for point in intersection])
#     return intersections
#
# n = 15
# random_lines = [generate_random_line() for _ in range(n)]
# IP = find_intersection_points(random_lines)
# IP = np.delete(IP,0,axis=0)
#
# n = 0
# Pos = {}
# E = []
# for i in range(len(IP)):
#     cen = np.array([IP[i][0],IP[i][1]])
#     vec1 = np.array([IP[i][2],IP[i][3]])
#     vec2 = np.array([IP[i][4],IP[i][5]])
#     Pos[4*n] = cen + 0.01*vec1
#     Pos[4*n+1] = cen - 0.01*vec1
#     Pos[4*n+2] = cen + 0.01*vec2
#     Pos[4*n+3] = cen - 0.01*vec2
#     E.extend([(4*n,4*n+2),(4*n,4*n+3),(4*n+1,4*n+2),(4*n+1,4*n+3)])
#     n += 1
#
# g = nx.Graph()
# g.add_nodes_from(Pos.keys())
# g.add_edges_from(E)
# nx.draw(g,pos=Pos,node_size=10,node_color='red',edge_color='red')
#
# plt.scatter(IP[:,0],IP[:,1])
# plot_lines(random_lines)






