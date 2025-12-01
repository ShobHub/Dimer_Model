import collections
import networkx as nx
import matplotlib.pyplot as plt  #import matplotlib.pylab as plt
from scipy.stats import linregress
from matplotlib.pyplot import pause
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
import pickle
import random
import pylab

class EWorm():
	def __init__(self,M,PL_HC,g,e,T,V):
		self.M = M            # A perfect dimer matching
		self.PL_HC = PL_HC    # List of Plaquettes
		self.g = g            # Networkx graph for a latiice
		self.bn,self.tn = nx.bipartite.sets(g)
		self.e = e
		self.Nodes = list(g.nodes())
		self.T = T           # Temperature
		self.V = V           # Potential

	def plotG(self,s0,si,sj):
		''' For Testing: plot lattice with a given dimer matching and test monomers at site {s0,s1,sj} during the worm algorithm workings.'''
		
		EC=[]
		W=[]
		for u,v in self.g.edges():
			if(((u,v) in self.M or (v,u) in self.M) and ((u,v) in self.e or (v,u) in self.e)):
				EC.append((1,0,0,0.5))
				W.append(2.5)
			elif((u,v) in self.M or (v,u) in self.M):
				EC.append('purple')
				W.append(1.8)
			else:
				EC.append('black')
				W.append(0.8)
		NC = []
		w = []
		for u in self.g.nodes():
			if (u==si):
				NC.append('red')
				w.append(30)
			elif (u==sj):
				NC.append('green')
				w.append(30)
			elif (u==s0):    #in s0):
				NC.append('blue')
				w.append(30)
			else:
				NC.append('steelblue')
				w.append(2)
		# NC1 = []
		# w1 = []
		# for u in G.nodes():
		# 	if u in BN:
		# 		NC1.append('blue')
		# 		w1.append(30)
		# 	elif u in TN:
		# 		NC1.append('green')
		# 		w1.append(30)
		# 	else:
		# 		NC1.append('b')
		# 		w1.append(0)
		pylab.ion()
		plt.clf()
		#nx.draw(G, pos=Pos, node_size=w1,width=0.8,node_color=NC1) #,with_labels=True,font_size=8)
		nx.draw(self.g,pos=pos,node_size=w,edge_color=EC,width=W,node_color=NC) #,with_labels=True,font_size=6)
		#plt.savefig('V=5_(10,10).pdf')
		pause(0.01)
		pylab.show()
		# plt.show()

	def NoFP(self, pl, mm):
		''' 
  		mm: current dimer edge
  		Calculate number of flippable plaquttes in a given list of plaquettes 'pl', when there is a dimer on 'mm'.
  		Return number of flippable plaquttes (c) and list of flippable plaquttes (fpl)
    	'''
		c = 0
		fpl = []
		for i in pl:
			mn = mm.copy()
			for a, b in i:
				if ((a, b) in self.M or (b, a) in self.M):
					mn.extend([a, b])
			if set(collections.Counter(mn).values()) == {1} and len(set(mn)) == len(i):
				c += 1
				fpl.append(i)
		return c,fpl

	def weights(self,u,v):
		'''
  		Calculate the number of plaquettes attached to the current dimer edge (u, v) and store the result in 'pl'. Return number of flippable plaquttes in 'pl'.
    	'''
		pl = []
		for i in self.PL_HC:
			if((u,v) in i or (v,u) in i):
				pl.append(i)
		# c = 0
		# for i in pl:
		# 	ver = []
		# 	for a, b in i:
		# 		ver.extend([a, b])
		# 	a, b = set(ver).difference(set([u, v]))
		# 	if ((a, b) in self.M or (b, a) in self.M):
		# 		c += 1
		# if(self.NoFP(pl,[u,v])[0]==c):
		# 	print('Yes')
		# return c
		return self.NoFP(pl,[u,v])[0]

	def InM(self, u):
		'''Returns whether an edge contains a dimer or not.'''
		for p,q in self.M:
			if (p == u):
				return (p,q),q
			elif (q == u):
				return (p,q),p

	def check_tuples_in_lists(self,list_of_tuples, list_of_list_of_tuples):
		normalized_target = set(tuple(sorted(t)) for t in list_of_tuples)
		for sublist in list_of_list_of_tuples:
			normalized_sublist = set(tuple(sorted(t)) for t in sublist)
			if normalized_target.issubset(normalized_sublist):
				return True
		return False

	def OrdP_MP(self,PL_BP,PL_HC):
		'''
   		Only for the Trivalent Penrose case:
    		Calculates the number of flippable plaquettes corresponding to each vertex degree 
    		in the dual standard Penrose lattice separately. This helps determine which plaquette 
    		type is most prevalent in the columnar state.
    	Also calculates the number of flippable rectangular plaquettes separately.
		'''
		NA = 0
		NB = 0
		bpp = []
		tpp = []
		DV = {3:0,4:0,5:0,6:0,7:0,8:0}
		TV = {3:0,4:0,5:0,6:0,7:0,8:0}
		for i in PL_BP.keys():
			c = self.NoFP([PL_BP[i]], [])[0]
			TV[G.degree(i)] += 1
			if c==1:
				DV[G.degree(i)] += 1
			# if c==1 and i in BN:
			# 	NA += 1
			# 	bpp.append(i)
			# elif c==1 and i in TN:
			# 	NB += 1
			# 	tpp.append(i)
		ts = 0
		fs = 0
		for i in PL_HC:
			if len(i) == 4:   #self.check_tuples_in_lists(i,PL_sq):
				ts += 1
				if self.NoFP([i], [])[0] == 1:
					fs += 1
		# print(TV,DV,fs,ts) #DV[3]+fs,TV[3]+ts)        #TV,DV,
		return (DV[3]+fs)/(TV[3]+ts)    #fs/ts    #

	# def OrdP(self):
	# 	NV = 0
	# 	NH = 0
	# 	for u,v in self.M:
	# 		if (round(abs(pos[u][0] - pos[v][0]), 1) == 0.0):
	# 			NV += 1
	# 		elif(round(abs(pos[u][1] - pos[v][1]), 1) == 0.0):
	# 			NH += 1
	# 	d = abs(NH - NV)*(2/(L**2))
	# 	return d

	def worm_single(self,k,i,i0,visit):
		''' 
		Implement a single worm update until the test monomers recombine. At the end of a worm self.M (dimer state) is updated.
		'''
		Mr = []
		count = 0
		while k != i0:
			dim, j = self.InM(i)
			self.M.remove(dim)
			# self.plotG(i0,i,j)
			# print(obj.NoFP(PL_HC, []))
			empE = list(self.g.edges(j))
			Wei = {}
			for u, v in empE:      # calculate weights for each empty edge
				w = np.exp(-self.V * self.weights(u, v) / self.T)
				if ((u, v) == dim or (v, u) == dim):
					rw = w
				else:
					Wei[(u, v)] = w
			wei = list(Wei.values())
			S = sum(wei) + rw
			n = random.uniform(0, 1)
			#cn = [wei[0] / S, (wei[0] + wei[1]) / S, (wei[0] + wei[1] + wei[2]) / S]

			# Calculate probaility for each empty edge
			cn = []                     
			for kk in range(1,len(wei)+1):
				su = 0
				for kk1 in range(kk):
					su += wei[kk1]
				cn.append(su/S)
			for i in range(len(wei)):     #3):
				if (n < cn[i]):
					a, b = list(Wei.keys())[list(Wei.values()).index(wei[i])]
					flag = True
					break
				else:
					flag = False
			if (flag == False):
				count += 1
				if count > 10:                             #For Staggered
					a,b = random.choice(empE)              #For Staggered
				else:
					a,b = dim
			self.M.append((a, b))
			visit.append(j)
			k = next(iter(set((a, b)).difference(set((j, -10)))))
			i = k
			if (j in self.bn and i0 in self.tn) or (j in self.tn and i0 in self.bn):
				Mr.append(round(np.linalg.norm(pos[j] - pos[i0])/S,2))
		return visit, Mr

	def worm(self,sweeps):
		'''Implement multiple sweeps until convergence. Each sweep consists of single worm updates that continue until most (>85%) of the nodes in the lattice are visited.'''
		Mlist = []
		while sweeps>0:
			print(sweeps)
			visit = []
			#D.append(self.OrdP_MP(PL_BP,PL_HC))
			#D.append(self.OrdP())
			Nfp.append(self.NoFP(self.PL_HC,[])[0])  #/len(PL_HC))
			while len(set(self.Nodes).difference(set(visit))) >= 0.15*len(self.Nodes):
				i0 = random.choice(self.Nodes)
				i = i0
				k = -1
				visit.append(i0)
				visit,_ = self.worm_single(k,i,i0,visit)
			sweeps -= 1
			Mlist.append(self.M.copy())
			#self.plotG(0, 0, 0)
		#D.append(self.OrdP())
		#D.append(self.OrdP_MP(PL_BP,PL_HC))
		Nfp.append(self.NoFP(self.PL_HC, [])[0])   #/len(PL_HC))
		return D,Nfp,self.M,Mlist

	def MMC(self,Mm,st):
		'''Calculate monomer-monomer correlations using a single worm update. After convergence, the histogram of test monomer separation distances yields the correlation function M(r, r_0).'''
		MR = []
		VST = []
		avgWL = 0
		for i in range(30):
			# print(i)
			self.M = Mm.copy()
			i0 = st  #random.choice(self.Nodes)
			i = i0
			k = -1
			visit = [i0]
			vst,Mr = self.worm_single(k, i, i0, visit)
			MR.extend(Mr)
			VST.extend(vst)
		# MR = collections.Counter(MR)
		# MR = {u:v/30 for u,v in dict(sorted(MR.items())).items()}
		VST = collections.Counter(VST)
		VST = {u:v/30 for u,v in VST.items()}
		MR = {}
		for u in TNmp:
			R = np.round(np.linalg.norm(pos[u] - pos[st]),3)
			if u in VST.keys():
				if R not in MR.keys():
					MR[R] = [VST[u]]
				else:
					MR[R].append(VST[u])
			else:
				if R not in MR.keys():
					MR[R] = [0]
				else:
					MR[R].append(0)
		MR = {u: np.mean(v) for u, v in MR.items()}
		# MR = {np.round(np.linalg.norm(pos[u]-pos[st]),1):VST[u] if u in VST.keys() else 0 for u in TNmp}   #else 0
		MR = dict(sorted(MR.items()))
		return VST,MR

## LOAD FILES
g = pickle.load(open('../3DegLattices/SqOct/SqOct12_OBC.txt','rb'))   #Hexa/Hexa24_PBC.txt','rb'))   # NetworkX graphs
pos = pickle.load(open('../3DegLattices/SqOct/Pos_SqOct12_OBC.txt','rb'))  #Hexa/Pos_Hexa24_PBC.txt','rb'))   # Node positions Dictionary
pos = {u:np.array(v) for u,v in pos.items()}
PL_HC = pickle.load(open('../3DegLattices/SqOct/Plaqs_SqOct12_OBC.txt','rb'))  #Hexa/Plaqs_Hexa24_PBC.txt','rb'))  #Plaquette List 

# G = pickle.load(open('../Cen_Tiling_(k,R)=(10,10).txt','rb'))
# Pos = pickle.load(open('../Cen_Pos_(k,R)=(10,10).txt','rb'))
# Pos = {u:v*4 for u,v in Pos.items()}
# BN, TN = nx.bipartite.sets(G)
#
# g = pickle.load(open('../ModifiedTiling/Cen_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# pos = pickle.load(open('../ModifiedTiling/Cen_Pos_Mod_Tiling_(k,R)=(10,10).txt','rb'))
# PL_HC = pickle.load(open('../ModifiedTiling/Cen_PlaquetteCycles_(k,R)=(10,10).txt','rb'))
# PL_BP = pickle.load(open('../ModifiedTiling/Cen_PlaqBPCharge_(k,R)=(10,10).txt','rb'))

# g = pickle.load(open('../all_graphs_plantri/G_out_20_15k.txt','rb'))
# pos = pickle.load(open('../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
# PL_HC = pickle.load(open('../all_graphs_plantri/PL_out_20_15k.txt','rb'))
# print(len(g.edges))
# fdg

# g = pickle.load(open('../RandomG/G_L=200_tutte.txt','rb'))
# pos = pickle.load(open('../RandomG/Pos_L=200_tutte.txt','rb'))
# PL_HC = pickle.load(open('../RandomG/Plaqs_L=200_tutte.txt','rb'))
# PL_sq = pickle.load(open('../RandomG/PlaqSQ_L=80.txt','rb'))

# nx.draw(g,pos,node_size=0,width=1.5,with_labels=True)
# nx.draw(G,Pos,node_size=0,width=0.7,edge_color='gray') #, with_labels=True)
# plt.show()

BNmp, TNmp = nx.bipartite.sets(g)
print(len(g.nodes()), len(PL_HC), len(BNmp),len(TNmp), 622 in BNmp)  #864,64(157),18     PL=816,443,19

e = []
M = nx.bipartite.hopcroft_karp_matching(g)
ME = []
for u,v in M.items():
	if((u,v) not in ME and (v,u) not in ME):
		ME.append((u,v))

V = -10
# T = [*np.linspace(0.02,15,40),*np.linspace(16,50,20)]
# T1,Cv1 = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (0.01),'rb'))
T2,Cv2 = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (-10),'rb'))
# T3,Cv3 = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (-30),'rb'))
# T4,Cv4 = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (10),'rb'))
# T5,Cv5 = pickle.load(open('../ModifiedTiling/FinalStates/CvNvsT_(k,R)=(10,10)_V=%.2f.txt'% (30),'rb'))
T = T2

'''
# T = 0.02
# V,_ = pickle.load(open('../ModifiedTiling/FinalStates/NFPvsV_T=0.02.txt','rb'))
# print(V)
# V = [-1.0,-0.5,-0.35,-0.2,-0.16111,-0.13,-0.12222,-0.11444,-0.09889,-0.08333,-0.06778,-0.05,-0.04444,-0.03667,-0.02111,-0.00556,0.01000,0.02000,0.02500,0.04375,0.06250,0.08125,0.1,0.2,0.5,0.75,1]
CN = {}
for i in T[:0]: #V[:0]:
    # print(i,'*')
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    Mm = pickle.load(open('../3DegLattices/SqOct/FS/FS_sqoct12_T=%.2f_V=%.2f.txt'%(i,V),'rb'))  #ModifiedTiling/FinalStates/FS_(k,R)=(10,10)_T=%.2f_V=%.5f.txt' %(0.0,i), 'rb'))    #   RandomG/FS ModifiedTiling/FinalStates (k,R)=(10,10)
    #Mm = pickle.load(open('../CorreFns/Stag_Plantri_20_15k.txt', 'rb'))
    D = []
    Nfp = []
    cn = []
    for j in [622]:
        obj = EWorm(Mm,PL_HC,g,e,T=i,V=V)
        VST, MR = obj.MMC(Mm,j)
        # print(MR)

        # print(len(VST.keys())/len(g.nodes()))
        # if len(VST.keys())/len(g.nodes()) < 0.4:
        #     print('*')
        #     cn.append(j)
        pylab.ion()
        plt.clf()
        # node_colors = [VST.get(node, 0) for node in g.nodes()]
        # norm = plt.Normalize(min(node_colors), max(node_colors))
        # cmap = plt.cm.Reds
        # nx.draw(g,pos,node_size=15,node_color=node_colors,cmap=cmap,width=1,label='T=%0.2f'%i) #,with_labels=True) #,ax=ax)
        pickle.dump(MR, open('../3DegLattices/SqOct/MR/Mr_sqoct12_T=%.2f_V=%.2f.txt' %(i,V), 'wb'))  #(k,R)=(10,10)  CorreFns/Data_ModP/MR/Mr_(k,R)=(10,10)_T=%.2f_V=%.5f.txt
        # sns.histplot(sorted(MR),kde=True,label='V=%0.5f'%i)   #binwidth=10
        plt.scatter(MR.keys(),MR.values(),s=20,label='V=%0.5f'%i)
        plt.plot(MR.keys(), MR.values())
        # plt.xscale('log')
        # plt.yscale('log')
        # plt.ylim(-100,3)
        plt.legend()
        pause(0.5)
        pylab.show()
        # plt.show()
    #CN[i] = cn
'''
'''
# T = [0.002, 0.02, 0.5, 1, 2,3.0, 3.7586206896551726,10.586206896551724, 11.344827586206897, 12.86206896551724, 13.620689655172413, 14.379310344827585, 15.137931034482758, 15.89655172413793, 16.655172413793103, 17.413793103448278, 18.17241379310345, 39.33333333333333, 42.0, 44.666666666666664, 47.33333333333333]
# T = [2.7087179487179487, 3.861025641025641, 5.013333333333333, 5.397435897435897, 5.781538461538461,6.165641025641025, 6.549743589743589, 6.9338461538461535,7.317948717948718,9.238461538461538,12.31128205128205,26.736842105263158,35.684210526315795]
v1 = [-1,-0.5,0.08333333333333326,-0.5,-1] # 0.08333333333333326 0.11666666666666667
v2 = [-1,-0.5,0.016666666666666663,0.5,1] #0.05000000000000002 0.016666666666666663

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(1, 1, figsize=(48, 48), sharex=True, sharey=True)
norm = mcolors.Normalize(vmin=min(v2), vmax=max(v2))
cmap = cm.get_cmap('rainbow')

for i in range(len(v1)):  #V[19:-14]:
	color = cmap(norm(v2[i]))
	MR = pickle.load(open('../CorreFns/Data_RandG/MR/Mr_L=200_T=%.2f_V=%.5f.txt' %(T,v1[i]), 'rb'))  #(k,R)=(10,10)
	mr = {}
	for u,v in MR.items():
		if round(u,3) in mr:
			mr[round(u,3)].append(v)
		else:
			mr[round(u,3)] = [v]
	mr = {a: np.mean(b) for a,b in mr.items()}
	MAX = max(list(mr.values()))
	plt.scatter(list(mr.keys()), list(mr.values())/MAX, s=1000, color = color, label='V=%.2f' % (v2[i]))
	plt.plot(list(mr.keys()), list(mr.values())/MAX, color = color, linewidth = 5)

plt.xlabel(r'Separation $(r)$',fontsize = 200)
plt.ylabel(r'$M(r)$',fontsize = 200)
plt.tick_params(axis='x', labelsize=130,length=20, width=2)
plt.tick_params(axis='y', labelsize=130,length=20, width=2)
plt.legend(loc='upper right', fontsize = 140) #, ncol=5)
plt.title(r'Random $(T = 0.02)$',fontsize = 200)
# plt.xscale('log')
# plt.yscale('log')
plt.ylim(bottom=-0.02)
plt.xlim(left=0)
plt.savefig('../../papers/RandomFieldTheory/Mr_Random_T=0.02.pdf')
# plt.show()
'''

XD = []
Dt = []
Nfpt = []
Cv = []
EC = []
err = []
swps = [100]*200

''' ACROSS DIFFERENT TEMPERATURE RANGE - FINDING GROUND STATES'''
for k in range(0):  #len(T)-1,-1,-2):
	print(V,T[k])
	D = []    # Order Parameter
	Nfp = []  # Number of flippable plaquette
	# ME = pickle.load(open('../3DegLattices/Hexa/FS/FS_hexa24_T=%.2f_V=%.2f.txt'%(T[k],V), 'rb'))
	# ME = pickle.load(open('../RandomG/FS/FS_L=200_T=%.2f_V=%.2f.txt'%(0.53,-10), 'rb'))
	obj = EWorm(ME,PL_HC,g,e,T=T[k],V=V)
	D,Nfp,Mm,Mlist = obj.worm(swps[k])
	# obj.plotG(0, 0, 0)
	pickle.dump(Mm, open('../3DegLattices/Hexa/FS/FS_hexa24_T=%.2f_V=%.2f.txt'%(T[k],V), 'wb'))  #CorreFns/Data_ModP  (k,R)=(10,10) RandomG/FS  ModifiedTiling/FinalStates
	print(len(Mlist),Mlist[0]==Mlist[-1],Mlist[-1]==Mlist[45],Mlist[45]==Mlist[-7])
	pickle.dump(Mlist[30:], open('../3DegLattices/Hexa/FS/Mlist_hexa24_T=%.2f_V=%.2f.txt' % (T[k],V), 'wb'))
	# pickle.dump(D, open('../ModifiedTiling/FinalStates/Dconv1_(k,R)=(10,10)_T=%.2f_V=%.2f.txt' % (T[k],V), 'wb'))
	pickle.dump(Nfp, open('../3DegLattices/Hexa/FS/NFP_hexa24_T=%.2f_V=%.2f.txt' % (T[k],V),'wb'))           #../CorreFns/Data_Plantri/NFP_20_15k_T=%.2f_V=%.5f.txt'%(T[k],V), 'wb'))
	# D = pickle.load(open('../ModifiedTiling/FinalStates/Dconv_(k,R)=(10,10)_T=%.2f_V=%.2f.txt' % (T[k],V), 'rb'))
	# Nfp = pickle.load(open('../RandomG/FS/NFP1nn_L=200_T=%.2f_V=%.2f.txt' % (T[k],V), 'rb'))
	# dt = np.array(Nfp)[25:]/len(PL_HC)  #int(len(Nfp)/2)
	# XD.append((np.mean(dt**2) - np.mean(dt)**2))    #Penrose L = diameter = 80
	E = V*np.array(Nfp)[30:]
	EC.append(1 - np.mean(E**4)/(3*np.mean(E**2)**2))
	Cv.append((np.mean(E**2) - np.mean(E)**2) / (1152*T[k]**2))  # np.mean(E)/1720  EC   1720,348(884),36
	# Dt.append(np.mean(D[-int(len(D)/2):]))
	Nfpt.append(np.mean(Nfp[30:]))
	print(Nfpt)
	err.append(np.std(Nfp[30:])/np.sqrt(len(Nfp[30:])))
	pylab.ion()
	# plt.clf()
	# plt.plot(list(range(len(D))), D, label='T=%0.5f'%T[k])
	# plt.scatter(list(range(len(D))), D, s=5)
	# plt.legend()
	# plt.ylim(0, 1.1)
	# pause(0.1)
	plt.clf()
	plt.plot(list(range(len(Nfp))), Nfp, label='T=%0.2f' % T[k])
	plt.scatter(list(range(len(Nfp))), Nfp, s=5)
	plt.legend()
	pause(0.1)
	pylab.show()
	# plt.ylim(0,1) #len(PL_HC))
	# plt.show()

print(Nfpt)
print(err)
print(Cv)
print(EC)








''' PLOTTING FOR THESIS '''
# T, Nfpt,err = pickle.load(open('../3DegLattices/Hexa/NFPvsT_hexa24_V=%.2f.txt' % (V),'rb'))
# T,Cv = pickle.load(open('../3DegLattices/Hexa/CvNvsT_hexa24_V=%.2f.txt' % (V),'rb'))
# T,EC = pickle.load(open('../3DegLattices/Hexa/ECvsT_hexa24_V=%.2f.txt' % (V),'rb'))

# pickle.dump([T[::-2][::-1],Nfpt[::-1],err[::-1]],open('../3DegLattices/Hexa/NFPvsT_hexa24_V=%.2f.txt' % (V),'wb'))
# pickle.dump([T[::-2][::-1],Cv[::-1]],open('../3DegLattices/Hexa/CvNvsT_hexa24_V=%.2f.txt' % (V),'wb'))
# pickle.dump([T[::-2][::-1],EC[::-1]],open('../3DegLattices/Hexa/ECvsT_hexa24_V=%.2f.txt' % (V),'wb'))

# plt.plot(T,Cv)
# plt.show()
# plt.plot(T,EC)
# plt.show()
# plt.errorbar(T, Nfpt, yerr=err, fmt='o', capsize=5)
# plt.show()

# D = []
# ME = pickle.load(open('../RandomG/FS/FS_L=200_T=%.2f_V=%.2f.txt'%(T[11],V),'rb'))
# for i in T[:0]:
# 	mlist = pickle.load(open('../CorreFns/Data_RandG/Samples/Mlist_L=200_T=%.2f_V=%.2f.txt'%(i,V),'rb'))
# 	# print(len(mlist))
# 	cl = []
# 	for j in mlist:
# 		c = 0
# 		for u,v in j:
# 			if (u,v) in ME or (v,u) in ME:
# 				c += 1
# 			else:
# 				c -= 1
# 		cl.append(c/len(ME))
# 	D.append(np.mean(cl))

# pickle.dump([T,D],open('../RandomG/FS/OrdPvsT_L=200_V=%.2f.txt'%V,'wb'))
# plt.plot(T,D)
# plt.show()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, axs = plt.subplots(1, 1, figsize=(88,84), sharex=True, sharey=True)  #(60,60)

Str = [('../ModifiedTiling/FinalStates/1NFPvsT_(k,R)=(10,10)_V=%.2f.txt','Penrose',816,[6,0]),('../RandomG/FS/1NFPvsT_L=200_V=%.2f.txt','Height gem',443,[6,0]),('../CorreFns/Data_Plantri/1NFPvsT_20_15k_V=%.2f.txt','Random',19,[6,0]),('../3DegLattices/SqOct/NFPvsT_sqoct12_V=%.2f.txt','Sq-Oct',529,[6,0]),('../3DegLattices/Hexa/NFPvsT_hexa24_V=%.2f.txt','HC',576,[3,0])]
# Str = [('../ModifiedTiling/FinalStates/XDvsT_(k,R)=(10,10)_V=%.2f.txt','Penrose',160),('../RandomG/FS/XDvsT_L=200_V=%.2f.txt','Random',80),('../CorreFns/Data_Plantri/XDvsT_20_15k_V=%.2f.txt','Plantri',20)]
Str = [('../ModifiedTiling/FinalStates/1CvNvsT_(k,R)=(10,10)_V=%.2f.txt','Penrose'),('../RandomG/FS/1CvNvsT_L=200_V=%.2f.txt','Height gem'),('../CorreFns/Data_Plantri/1CvNvsT_20_15k_V=%.2f.txt','Random'),('../3DegLattices/SqOct/CvNvsT_sqoct12_V=%.2f.txt','Sq-Oct'),('../3DegLattices/Hexa/CvNvsT_hexa24_V=%.2f.txt','HC')]
# Str = [('../ModifiedTiling/FinalStates/1ECvsT_(k,R)=(10,10)_V=%.2f.txt','Penrose'),('../RandomG/FS/1ECvsT_L=200_V=%.2f.txt','Height gem'),('../CorreFns/Data_Plantri/1ECvsT_20_15k_V=%.2f.txt','Random'),('../3DegLattices/SqOct/ECvsT_sqoct12_V=%.2f.txt','Sq-Oct'),('../3DegLattices/Hexa/ECvsT_hexa24_V=%.2f.txt','HC')]
# Str = [('../ModifiedTiling/FinalStates/OrdPvsT_(k,R)=(10,10)_V=%.2f.txt','Penrose'),('../RandomG/FS/OrdPvsT_L=200_V=%.2f.txt','Random'),('../CorreFns/Data_Plantri/OrdPvsT_20_15k_V=%.2f.txt','Plantri')]
# Str = [('../ModifiedTiling/FinalStates/NFPvsV_T=%.2f.txt','Penrose',816),('../RandomG/FS/NFPvsV_L=200_T=%.2f.txt','Random',443),('../CorreFns/Data_Plantri/NFPvsV_20_15k_T=%.2f.txt','Plantri',19)]

V = [-10] #,10]
maxF = {'Penrose':640,'Height gem':344,'Random':14}
for i in range(len(V)):
	for s,name in Str:  #,L,lim
		T, Cv = pickle.load(open('%s'%s%V[i], 'rb'))   #yerr
		# plt.plot(np.array(T[:len(T)-lim[i]])/10, np.array(Cv[:len(T)-lim[i]])/L, linewidth = 20)    #axs[i]
		# plt.scatter(np.array(T[:len(T)-lim[i]])/10, np.array(Cv[:len(T)-lim[i]])/L, s=4000, label='V=%.2f %s' % (V[i]/10, name))
		# plt.errorbar(np.array(T[:len(T) - lim[i]])/10, np.array(Cv[:len(T)-lim[i]])/L, yerr = np.array(yerr[:len(T) - lim[i]])/L, fmt = '-o', linewidth = 5, capsize=60, markersize=40, markerfacecolor='none',label='V=%.2f %s' % (V[i]/10, name))

		plt.plot(np.array(T[1:])/10,Cv[1:],linewidth = 20)   #Cv[1:]     #CV, EC
		plt.scatter(np.array(T[1:])/10,Cv[1:],s=4000,label='V=%.2f %s' % (V[i]/10, name))

		# plt.plot(T, np.array(Cv)/L, linewidth=10)
		# plt.scatter(T, np.array(Cv)/L, s=1000, label='T=%.2f %s' % (V[i], name))
		# if i == 0 and name in ['Penrose','Height gem','Random']:
		# 	line = plt.axhline(y=maxF[name]/L, color = 'gray', linestyle='--', linewidth=25)   #, label='analytical (%s)'%name)
	# plt.legend(loc='upper center', bbox_to_anchor=(0.48, 1.17), fontsize=199, ncol=4, columnspacing=0.1, handletextpad=0.002, frameon=False)
	plt.legend(loc='upper right', fontsize=200)
	plt.ylim(bottom=0,top=0.4)   #,top=0.4)    #,top=0.3)  #0.62,top=0.67)    #0,top=0.25)   #0.07)
	plt.xlim(left=0) #,right=4)


# Str = [('../CorreFns/Data_ModP/FKT_Mr_(k,R)=(10,10).txt','Penrose'),('../CorreFns/Data_RandG/FKT_Mr_L=200.txt','Random'),('../CorreFns/Data_Plantri/FKT_Mr_20_15k.txt','Plantri'),('../3DegLattices/SqOct/FKT_Mr_sqoct12_obc.txt','Sq-Oct')] #,('../3DegLattices/Hexa/FKT_Mr_hexa30_obc.txt','HC')]
# for s,name in Str[:]:
# 	sort_M= pickle.load(open('%s'%s, 'rb'))
# 	if name == 'Penrose':
# 		x = (np.array(list(sort_M.keys())) - min(list(sort_M.keys()))) / 90
# 		y = np.array(list(sort_M.values()))/max(list(sort_M.values()))
# 	elif name == 'Sq-Oct':
# 		x = (np.array(list(sort_M.keys())) - min(list(sort_M.keys()))) / 40
# 		y = np.array(list(sort_M.values()))/max(list(sort_M.values()))
# 	else:
# 		x = np.array(list(sort_M.keys())) - min(list(sort_M.keys()))
# 		y = np.array(list(sort_M.values())) / max(list(sort_M.values()))
# 	plt.scatter(x, y, s = 3000, label='%s'%name)
# 	plt.plot(x, y,linewidth = 15)
# 	plt.xlim(left=0,right=0.95)
# 	plt.title(r'$T = \infty$',fontsize=300)
#	plt.legend(loc='upper right', fontsize = 200)

# norm = mcolors.Normalize(vmin=min(T[5:-23:2]), vmax=max(T[5:-23:2]))
# cmap = cm.get_cmap('rainbow')
#
# Str1 = [('../CorreFns/Data_ModP/MR/Mr_(k,R)=(10,10)_T=%.2f_V=%.5f.txt','Penrose'),('../CorreFns/Data_RandG/MR/Mr_L=200_T=%.2f_V=%.5f.txt','Random'),('../CorreFns/Data_Plantri/MR/Mr_20_15k_T=%.2f_V=%.5f.txt','Plantri'),('../3DegLattices/SqOct/MR/Mr_sqoct12_T=%.2f_V=%.2f.txt','Sq-Oct')]
# V = [-10] #,0.01,30]
# J = [-60,-20,-1]
# for i in range(3,4):
# 	for j in range(len(V)):
# 		T,_ = pickle.load(open('%s'%Str[i][0]%V[j], 'rb'))
# 		for k in T[0:-30:3]:  #[18:-23:2]:  #-5:J[j]]:
# 			color = cmap(norm(k))
# 			MR = pickle.load(open('%s'%Str1[i][0]%(k,V[j]), 'rb'))
# 			x = np.array(list(MR.keys())) # * factor[i] # - min(list(MR.keys()))
# 			y = np.array(list(MR.values()))
# 			plt.scatter(x/40, y, s=3000, color = color, label='T=%0.2f'%(k))  #/90
# 			plt.plot(x/40, y,color = color, linewidth = 15)
# 			plt.title(r'%s $(v=%.2f)$'%(Str[i][1],-1), fontsize=300)   #V[j]))
# 	plt.legend(loc='upper right', fontsize=240, ncol=2, columnspacing=0.1,handletextpad=0.002, frameon=False) #bbox_to_anchor=(0.48, 1.18)
# 	plt.xlim(left = 0)
# 	plt.ylim(bottom=-0.02)

plt.xlabel(r'Temperature $(T)$',fontsize = 320)  #280)    #   #fig.supylabel   Temperature $(T)$  Separation $(r)$
plt.ylabel(r'$C_v/N$',fontsize = 300) #300)   #Exact Monomer Correlations (FKT)     $\langle NFP \rangle$  $C_v/N$  $M(r)$
plt.tick_params(axis='x', labelsize=260,length=70, width=15)  #200,length=100, width=10)  # Increase x-axis numbers size
plt.tick_params(axis='y', labelsize=260,length=70, width=15)
# plt.xscale('log')
# plt.yscale('log')
plt.savefig('../../papers/RandomFieldTheory/CVvsT_thesis_1.pdf')  #Mr_SqOct_thesis.pdf')
# plt.show()







'''
T,Nfpt = pickle.load(open('../ModifiedTiling/FinalStates/NFPvsV_(k,R)=(10,10)_V=%.2f.txt'% (-0.2),'rb'))    #OPvsV_T=0.02_ModT_(k,R)=(10,10)_1.txt','rb'))


def PlaqOrdP(self, m, L):
	sum = 0
	self.M = m.copy()
	for i in PL_HC:
		cen = 0
		for j in i:
			cen += pos[j[0]]
		cen = cen / len(i)
		sum += (-1) ** int((cen[0] + cen[1])) * self.NoFP([i], [])[0]
	return abs(sum) * (2 / L ** 2)
Pavg = []
for i in T:
	Mlist  = [pickle.load(open('../ModifiedTiling/FinalStates/FS_(k,R)=(10,10)_T=%.2f_V=%.5f.txt' %(i,V), 'rb'))]
	P = 0
	for m in Mlist:
		obj = EWorm(m, PL_HC, g, e, T=i, V=V)
		P += obj.PlaqOrdP(m,80)
	Pavg.append(P/len(Mlist))


# T, Nfpt, err = pickle.load(open('../RandomG/FS/1NFPvsT_L=200_V=%.2f.txt' % (V),'rb'))
# T, Cv = pickle.load(open('../RandomG/FS/1CvNvsT_L=200_V=%.2f.txt' % (V),'rb'))
# T, EC = pickle.load(open('../RandomG/FS/1ECvsT_L=200_V=%.2f.txt' % (V),'rb'))
#
# err = [0.02634788953262613, 0.01666708705724889, 0.008200369504755885, 0.005736825579568762, 0.002814287757973408, 0.0, 0.0, 0.0, 0.0][::-1] + err[9:]
# Nfpt = [343.0969030969031, 343.72327672327674, 343.9320679320679, 343.968031968032, 343.992007992008, 344.0, 344.0, 344.0, 344.0][::-1] + Nfpt[9:]
# Cv = [0.01868371381494571, 0.009737911098210983, 0.0031966669265480648, 0.0022412387350075305, 0.0008362582563365103, 0.0, 0.0, 0.0, 0.0][::-1] + Cv[9:]
# EC = [0.6666588104938718, 0.6666635379853894, 0.6666659100828716, 0.6666662965435596, 0.6666665775887952, 0.6666666666666667, 0.6666666666666667, 0.6666666666666667,0.6666666666666667][::-1] + EC[9:]
'''
