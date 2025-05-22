import numpy as np
import pickle
import time

g = pickle.load(open('../all_graphs_plantri/G_out_20_15k.txt','rb'))
pos = pickle.load(open('../all_graphs_plantri/Pos_out_20_15k.txt','rb'))
PL_HC = pickle.load(open('../all_graphs_plantri/PL_out_20_15k.txt','rb'))
print(len(PL_HC))

V = -1
T = 1

H_local = np.matrix([[V,-T],[-T,V]])
H_global = H_local.copy()
for i in range(1,1):  #len(PL_HC)):
    H_global = np.kron(H_global,H_local)

print(H_global,H_global.shape)

st = time.time()
eigVl,eigVc = np.linalg.eig(H_global)
print(time.time()-st)

# pickle.dump(eigVl,open('../CorreFns/Data_Plantri/EigValues_20_15k.txt','wb'))
# pickle.dump(eigVc,open('../CorreFns/Data_Plantri/EigVectors_20_15k.txt','wb'))

eigVl = pickle.load(open('../CorreFns/Data_Plantri/EigValues_20_15k.txt','rb'))
eigVc = pickle.load(open('../CorreFns/Data_Plantri/EigVectors_20_15k.txt','rb'))

print(eigVc[:,1])
print(sorted(eigVl))

