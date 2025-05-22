import matplotlib.pylab as plt
from scipy.spatial import KDTree
from shapely.geometry import Point, Polygon
import networkx as nx
import numpy as np
import pickle
import random

G = pickle.load(open('../Cen_Tiling_(k,R)=(20,20).txt','rb'))
Pos = pickle.load(open('../Cen_Pos_(k,R)=(20,20).txt','rb'))
Pos = {u:v*4 for u,v in Pos.items()}
PL = pickle.load(open('../Cen_Plaqs_(k,R)=(20,20).txt','rb'))
g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(20,20).txt','rb'))
pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(20,20).txt','rb'))
PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(20,20).txt','rb'))
PL_BP = pickle.load(open('../ModifiedTiling/Cen_PlaqBPCharge_(k,R)=(20,20).txt','rb'))
rev_pos = {tuple(v): k for k, v in pos.items()}
BN,TN = nx.bipartite.sets(G)

def contains_tuple(target_tuple, lst_of_lst_of_tuples):
    target_set = set(target_tuple)
    result = []
    for lst in lst_of_lst_of_tuples:
        for tpl in lst:
            if set(tpl) == target_set:
                result.append(lst)
                break
    return result
def find_opposite_vertex(vertex, edges):
    quad = nx.Graph()
    quad.add_edges_from(edges)
    quad.remove_node(vertex)
    return next(node for node, data in quad.degree() if data == 2)

def signed_angle_with_negative_y_axis(p1, p2):
    v, neg_y = np.array(p2) - np.array(p1), np.array([0, -1])
    angle = np.arccos(np.dot(v, neg_y) / np.linalg.norm(v))
    return np.degrees(angle if np.cross(neg_y, v) > 0 else -angle)

def G_lines(stp):
    visit = [2602,2723,7337,2673,2,194,8035,9017,9209]    #[442,24,8,2,1361,1579]
    #BD = [7337, 2673, 1, 32, 7049, 2720, 104, 13275, 193]
    #[3010,2746,1699,443,441,440,688,25,24,17,16,2739,697,9,8,425,1,0,705,707,2651,1819,1817,32,34,2329,682,80,1392,136,200,1377,2081,1275,2769,616,1456,1361,658]
    L = []
    for i in stp:
        st = i
        Lst = [st]
        while True:
            visit.append(st)
            neigh = G.neighbors(st)
            res = []
            for j in neigh:
                res.extend(contains_tuple((j,st),PL))
            unique_res = set(tuple(sorted(sb)) for sb in res)
            nxt = {}
            for j in [list(t) for t in unique_res]:
                a = find_opposite_vertex(st,j)
                if Pos[st][1]-Pos[a][1]>0 and a not in visit:  # and a not in BD:
                    nxt[a] = round(Pos[st][1]-Pos[a][1],2)
            if nxt == {}:
                break
            elif len(nxt)==1:
                st = list(nxt.keys())[0]
            else:
                max_value = max(nxt.values())
                max_pairs = [(key, value) for key, value in nxt.items() if value == max_value]
                if len(max_pairs) == 1:
                    st = max_pairs[0][0]
                else:
                    a = max_pairs[0][0]
                    b = max_pairs[1][0]
                    if Pos[a][0]<Pos[b][0]:
                        st = a
                    else:
                        st = b
            if abs(Pos[st][0]-Pos[Lst[-1]][0])<3:
                Lst.append(st)
            else:
                break
        if len(Lst)>1:
            #print(Lst)
            L.append(Lst)
    return L

def points_inside_polygon(polygon, points_list):
    polygon_coords = [tuple(Pos[polygon[0][0]]),tuple(Pos[polygon[1][0]]),tuple(Pos[polygon[2][0]]),tuple(Pos[polygon[3][0]])]
    polygon = Polygon(polygon_coords)
    points_array = np.array(points_list)
    shapely_points = np.apply_along_axis(lambda point: Point(point), 1, points_array)
    inside_mask = np.vectorize(polygon.contains)(shapely_points)
    inside_points = points_array[inside_mask]
    inside_nodes = []
    for i in inside_points:
        inside_nodes.append(rev_pos[tuple(i)])
    return inside_nodes

def PenPlaq(a,b):
    for i in PL:
        l = list(set(element for tup in i for element in tup))
        if a in l and b in l:
            return i

def PlaqG_lines(L):
    LP = []
    for i in L:
        l = []
        for j in range(len(i)-1):
            pl = PenPlaq(i[j],i[j+1])
            sq = points_inside_polygon(pl,list(pos.values()))
            l.extend([PL_BP[i[j]],sq])
        #l.append(PL_BP[i[-1]])
        LP.append(l)
    return LP

def g_lines(SE,LP,E,q):
    visit = []
    for i in range(len(LP)):
        se1 = SE[i][q]
        e1 = [se1]
        visit.extend(list(se1))
        j = 0
        while j != len(LP[i]):  #-1:
            if se1[1] not in LP[i][j+1]:
                for k in g.neighbors(se1[1]):
                    if (se1[1],k) in LP[i][j] and k not in visit and k not in LP[i][j-1]:
                        se1 = (se1[1],k)
                        visit.append(k)
                        break
                    elif (k,se1[1]) in LP[i][j] and k not in visit and k not in LP[i][j-1]:
                        se1 = (se1[1],k)
                        visit.append(k)
                        break
            else:
                for k in g.neighbors(se1[1]):
                    if k in LP[i][j+1] and (se1[1],k) not in LP[i][j] and (k,se1[1]) not in LP[i][j] and k not in visit:
                        se1 = (se1[1],k)
                        visit.append(k)
                        break
                j = j+2
            e1.append(se1)
        E.append(e1)
    return E
def start_points(stp,BD):
    a = stp[1]
    pt = (Pos[a][0], 0)
    if stp[0] in BN and stp[1] in BN:
        pos_bp = {u:(v[0],0) for u,v in Pos.items() if u in BN}
        tree = KDTree(list(pos_bp.values()))
    else:
        pos_bp = {u:(v[0],0) for u,v in Pos.items() if u in TN}
        tree = KDTree(list(pos_bp.values()))
    kbp = list(pos_bp.keys())
    while True:
        nn = {}
        for i in tree.query_ball_point(pt, r=9.7):   #12.4):
            #print(kbp[i],Pos[kbp[i]][0]>pt[0])
            if Pos[kbp[i]][0]>=pt[0] and kbp[i] not in BD:
                nn[kbp[i]] = Pos[kbp[i]][1]
        #print(nn)
        if nn == {}:
            break
        else:
            max_value = max(nn.values())
            max_pairs = [(key, value) for key, value in nn.items() if value == max_value]
            if len(max_pairs) == 1:
                a = max_pairs[0][0]
            else:
                a = max_pairs[0][0]
                b = max_pairs[1][0]
                if Pos[a][0] > Pos[b][0]:
                    a = b
            #print(a)
            stp.append(a)
            pt = (Pos[a][0],0)
    print(stp)
    return stp

# def BDdimers(MEn):
#     bd = list(set(g.nodes()).difference(MEn))
#     meBD = []
#     while len(bd)!=0:
#         print(len(bd))
#         nd = random.choice(bd)
#         for i in g.neighbors(nd):
#             if i in bd:
#                 meBD.append((nd,i))
#                 bd.remove(i)
#                 break
#         bd.remove(nd)
#     return meBD

BD = [2602,7506,2640,5307,7658,3688,2643,3698,2025,11656,2729,6659,1675,97,2681,7792,7337,2673,1,32,7049,2720,104,13275,193,1737,3898,2792,6843,1873,521,6907,7946,11496,6915,3160,7962,3242,4714,7955,7952,9010,3394,1347,11128,3528,7914,1531,2594,6762,1595,2634,7792,1603,2650,9209]
#stp = start_points([2659,88],BD)

stp = [2659,88,98,184,288,291,6755,1874,522,10465,2026,760,763,883,3241,1003,1121,3393,1235,1344,1339,1443,1528,3592,1592,1579]  #,1539]
      #[27,73,2027,499,192,195,571,320,610,651,1755,675,945]

L = G_lines(stp)
LP = PlaqG_lines(L)
SE = [[(6479,6125),(6478,6484)],[(364,383),(367,360)],[(403,407),(402,415)],[(716,426),(719,712)],[(1124,1135),(1127,1120)],[(1143,1121),(1142,1148)],[(6206,6534),(6205,1582)],
      [(6208,6201),(6211,1998)],[(1995,2000),(1994,2007)],[(2011,2012),(2010,2497)],[(6252,6245),(6255,6256)],[(2920,2934),(2923,2916)],[(2945,2939),(2944,2940)],[(3403,3389),(3402,3412)],
      [(6584,3413),(6587,6604)],[(3896,3891),(3899,3900)],[(4293,4299),(4292,4281)],[(4315,4300),(4314,6624)],[(4776,4771),(4779,4780)],[(5196,5207),(5199,5192)],[(5191,5183),(5190,5184)],
      [(5535,5521),(5534,5544)],[(5864,5514),(5867,5860)],[(6447,5850),(6446,5845)],[(6088,5826),(6091,6084)],[(6055,6047),(6054,6048)],[(5903,5895),(5902,5896)]]
#[[(95,89),(94,97)],[(264,271),(267,261)],[(279,288),(278,1582)],[(1576,490),(1579,1584)],[(728,503),(731,720)],[(732,717),(735,744)],[(1612,954),(1615,969)],[(967,960),(966,1224)],[(1234,1229),(1233,1632)],[(1647,1378),(1646,1384)],[(1388,1380),(1391,1530)],[(1534,1655),(1533,1508)],[(1695,1696),(1694,1708)]]
#print(len(stp),len(SE))

E1 = g_lines(SE,LP,[],0)
E2 = g_lines(SE,LP,[],1)

ME = [(5924,5927),(5900,5903),(5895,5892),(5873,5874),(5899,5898),(6055,6052),(6079,6076),(6075,6072),(6058,6057),(6068,6071),(6080,6081)]
# MEn = set()
# for edge in ME:
#     MEn.update(edge)
for i in E1:
    for j in i[1::2]:
        ME.append(j)
        # MEn.update(j)
for i in E2:
    for j in i[::2]:
        ME.append(j)
        # MEn.update(j)

# ME.extend(BDdimers(MEn))

#pickle.dump(ME,open('Stag_Cen_Mod_Tiling_(k,R)=(20,20).txt','wb'))

Col = {}
for i in range(len(L)):
    Col[i] = np.random.choice(range(256), size=3)/256
NC = []
for i in G.nodes():
    flag = False
    for j in range(len(L)):
        if i in L[j]:
            NC.append(Col[j])
            flag = True
            break
    if flag==False:
        NC.append('white')
EC = []
W = []
ed = []
for i in E1:
    ed.extend(i)
for i in E2:
    ed.extend(i)
for u,v in g.edges():
    if (u,v) in ME or (v,u) in ME:    #ed   ME
        EC.append('red')
        W.append(2)
    else:
        EC.append('black')
        W.append(1)
nx.draw(G,Pos,node_size=20,node_color=NC,width=0.7,edge_color='gray') #,with_labels=True)
nx.draw(g,pos,node_size=0,edge_color=EC,width=W) #,with_labels=True)
plt.show()

