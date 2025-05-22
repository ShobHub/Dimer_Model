import matplotlib.pylab as plt
from scipy.optimize import curve_fit
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Polygon
from matplotlib.colors import Normalize
import scipy.special as sp
import scipy.optimize as opt
import networkx as nx
import pandas as pd
import numpy as np
import colorsys
import pickle
import pylab

G = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))

# g = pickle.load(open('../all_graphs_plantri/G_out_20_15k.txt','rb'))
# pos = pickle.load(open('../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
# PL_HC = pickle.load(open('../all_graphs_plantri/PL_out_20_15k.txt','rb'))

# g = pickle.load(open('../RandomG/G_L=80_tutte.txt','rb'))
# pos = pickle.load(open('../RandomG/Pos_L=80_tutte.txt','rb'))
# PL_HC = pickle.load(open('../RandomG/Plaqs_L=80_tutte.txt','rb'))

def plot(H,t):

    # X, Y = zip(*list(H.keys()))
    # plt.scatter(X, Y, s=1, c='black')
    # for i, label in enumerate(list(H.values())):
    #     plt.annotate(np.round(label, 1), list(H.keys())[i], ha='center', c='r')
    # nx.draw(G, pos, node_size=0,label='T=%.2f'%t)
    # # pause(5)
    # # pylab.show()
    # plt.legend()
    # plt.show()

    hs = list(set(H.values()))
    coms = sorted(hs)
    mycols = {}
    part = 1.0 / len(coms)
    for k in range(len(coms)):
        mycols[coms[k]] = colorsys.hsv_to_rgb(50, 0.4, k * part)

    fig, ax = plt.subplots()
    color_list = [mycols[val] for val in coms]
    norm = Normalize(vmin=min(coms), vmax=max(coms))
    sm = ScalarMappable(cmap=plt.cm.colors.ListedColormap(color_list), norm=norm)
    sm.set_array([])

    for i in PL_HC:
        plN = []
        x = 0
        y = 0
        for j in i:
            plN.extend(list(j))
            x += pos[j[0]][0]
            y += pos[j[0]][1]
            x += pos[j[1]][0]
            y += pos[j[1]][1]
        cx = x/(2*len(i))
        cy = y/(2*len(i))
        plN = set(plN)
        plN = [pos[i] for i in plN]
        def angle_from_center(point):
            x, y = point
            return np.arctan2(y - cy, x - cx)  # Angle in radians
        plN = sorted(plN, key=angle_from_center, reverse=True)
        if (cx,cy) in H.keys():
            ax.add_patch(Polygon(plN, closed=True, fill=True, facecolor=mycols[H[(cx,cy)]], edgecolor='black',linewidth=0.15))  #H[(x/len(i),y/len(i))]
        else:
            ax.add_patch(Polygon(plN, closed=True, fill=True, facecolor='white', edgecolor='black',linewidth=0.15))  #H[(x/len(i),y/len(i))]

    cbar = fig.colorbar(sm, ax=ax, shrink=0.6)
    cbar.set_label('<HH> Values')
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.set_title('T=%.2f' % t)
    ax.axis('off')
    plt.title('T = %.2f'%t)
    plt.show()

# StagH = pickle.load(open('../CorreFns/Data_ModP/Samples/StagH_(k,R)=(10,10).txt','rb'))   #/ModifiedTiling/Heights/Stag_Cen_(k,R)=(10,10).txt','rb'))
# plot(StagH,1)
# tg
T,Cv = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (-10),'rb'))
fix_cen = (-1.0000000000000004, -3.077683537175253) #ModP(10,10)    #(0.5671648285092823, 0.20789275642336102)RG_L=80        #(1.5776939943554096, 1.5726909306322203)AB_I=4
for i in T[:11:2]+T[11::5]:
    H_ = pickle.load(open('../CorreFns/Data_ModP/Samples/H_(k,R)=(10,10)_T=%.2f_V=-10.00.txt'%(i), 'rb'))
    H0 = []
    HR0 = []
    HR = []
    for j in range(len(H_)):
        c = 0
        hr0 = {}
        hr = {}
        h0 = H_[j][fix_cen]
        # del H_[j]['out']    #AB_I=4
        for k in H_[j].keys():
            r = np.linalg.norm(np.array(k) - np.array(fix_cen))
            hr0[k] = H_[j][k] * h0   #r
            hr[k] = H_[j][k]
        HR0.append(hr0)
        HR.append(hr)
        H0.append(h0)
    df = pd.DataFrame(HR0)
    HR0 = df.mean().to_dict()
    df = pd.DataFrame(HR)
    HR = df.mean().to_dict()
    H0 = np.mean(H0)
    Hcorre = {key: HR0[key] - HR[key] * H0 for key in HR0}
    plot(Hcorre,i)
    '''
    R = {}
    for a,b in Hcorre.items():
        if round(a,3) not in R.keys():      #2-ModP,  3-RG
            R[round(a,3)] = [b]
        else:
            R[round(a,3)].append(b)
    R = {key: np.mean(R[key]) for key in R}
    R = dict(sorted(R.items()))
    x = np.array(list(R.keys()))
    y = np.array(list(R.values()))
    # print(R)

    def bessel_fit_func(x, A, B):
        return A * sp.k0(B * x)
    popt, pcov = opt.curve_fit(bessel_fit_func, x[1:], y[1:], p0=[np.max(y[1:]), 1],maxfev=50000)  # Avoid zero in x_data to prevent singularity in K_0
    A_fit, B_fit = popt
    y_fit = bessel_fit_func(x[:], A_fit, B_fit)
    # data_to_save = np.row_stack((x, y, y_fit))
    # np.savetxt("../CorreFns/Data_RandG/Samples/<HH>_T=%.2f_V=-10.00.txt"%i, data_to_save)

    # data = np.loadtxt("../CorreFns/Data_ModP/Samples/<HH>_T=%.2f_V=-10.00.txt"%i)
    # Extract x, y, and y_fit from the loaded data
    # x = data[0, :]
    # y = data[1, :]
    # y_fit = data[2, :]

    plt.scatter(x, y, color="b", s=10, label='T=%.2f'%i)
    plt.plot(x, y, color="b")
    plt.plot(x, y_fit, label=f"Fit: A={A_fit:.2f}, B={B_fit:.2f}", color="r")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.title("Fit of Modified Bessel Function K_0")
    plt.show()
    '''











