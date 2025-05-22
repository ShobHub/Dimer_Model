from koala import plotting as pl
from koala import graph_utils as gu
from koala.functions.koala_plantri import plantri_to_koala, read_plantri
from Dimerisations import count_dimers,plot_dimers
from koala_to_NX import koala_to_nx

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle

filename = "../all_graphs_plantri/graphs_out_20"
lats = read_plantri(filename, verbose=True)
lattice = plantri_to_koala(lats[15300])

dimers = gu.dimerise(lattice,100)
PLSC,EXSC,MFI,LFI = count_dimers(dimers,lattice)
#plot_dimers(lattice,dimers,PLSC,EXSC,MFI,LFI)

fig, ax = plt.subplots(figsize = (10,10))
pl.plot_edges(lattice,dimers[LFI[0]],ax=ax)
plt.show()

Stag = []
stag = dimers[LFI[0]]
for i in range(len(stag)):
    if stag[i] == 1:
        Stag.append(tuple(lattice.edges.indices[i]))

pos,E,PL = koala_to_nx(lattice)
G = nx.Graph()
G.add_edges_from(E)

EC = []
W = []
for u,v in G.edges():
    if (u,v) in Stag or (v,u) in Stag:
        EC.append('g')
        W.append(1.7)
    else:
        EC.append('black')
        W.append(0.7)

print(len(PL))
nx.draw(G, pos, node_size=0, width = 2) #, edge_color = EC) #, with_labels=True)
plt.show()

# pickle.dump(G,open('../all_graphs_plantri/G_out_20_15k.txt','wb'))
# pickle.dump(pos,open('../all_graphs_plantri/Pos_out_20_15k.txt','wb'))
# pickle.dump(PL,open('../all_graphs_plantri/PL_out_20_15k.txt','wb'))

# pickle.dump(Stag,open('../Mij_correF/Stag_Plantri_20_15k.txt','wb'))







'''
import matplotlib.pylab as plt
from matplotlib.patches import Polygon
import collections
import networkx as nx
import random
import time
import pickle
import subprocess
import os

# Define the C source file and the output executable
source_file = 'plantri55/plantri.c'
output_executable = 'plantri55/hello'

# Compile the C code
compile_command = ['gcc', source_file, '-o', output_executable]
result = subprocess.run(compile_command, capture_output=True, text=True)

LV = 15
# Check if there were any errors during compilation
if result.returncode != 0:
    print("Compilation failed with the following message:")
    print(result.stderr)
else:
    print("Compilation successful.")
    # Run the compiled executable
    run_command = ['./' + output_executable,'-d','%i'%LV,'-b','-a','plantri55/outL.txt']   #'1/10000'
    result = subprocess.run(run_command, capture_output=True, text=True)
    if result.returncode != 0:
        print("Execution failed with the following message:")
        print(result.stderr)
    else:
        print("Execution output:")
        #print(result.stdout)

def read_graph_from_file(filename,Lnum):
    with open(filename, 'r') as file:
        line = file.readlines()[Lnum].strip()
        print(line)
        num_nodes, edges_data = line.split(' ', 1)
        num_nodes = int(num_nodes)
        edges_data = edges_data.split(',')
    return num_nodes, edges_data

def create_graph(num_nodes, edges_data):
    G = nx.Graph()
    for i in range(num_nodes):
        G.add_node(chr(ord('a') + i))
    for i, connections in enumerate(edges_data):
        node = chr(ord('a') + i)
        for connected_node in connections:
            G.add_edge(node, connected_node)
    return G

filename = 'plantri55/outL.txt'
Lnum = 0
num_nodes, edges_data = read_graph_from_file(filename, Lnum)
graph = create_graph(num_nodes, edges_data)

# M = nx.bipartite.hopcroft_karp_matching(graph)
# ME = []
# for u,v in M.items():
# 	if((u,v) not in ME and (v,u) not in ME):
# 		ME.append((u,v))
# print(ME)
#
# EC = []
# W = []
# for u,v in graph.edges():
#     if((u,v) in ME or (v,u) in ME):
#         EC.append('purple')
#         W.append(2)
#     else:
#         EC.append('black')
#         W.append(1)

print(nx.is_bipartite(graph),len(graph.nodes))
pos = nx.kamada_kawai_layout(graph)
nx.draw(graph,pos,node_size=1)  #,edge_color=EC,width=W,with_labels=True)
plt.show()
'''

'''
def Sread_graph_from_file(filename, Lnum):
    with open(filename, 'r') as file:
        lines = file.readlines()
        selected_line = lines[Lnum].strip()
    return selected_line

def decode_sparse6(s6):
    if not s6.startswith(':'):
        raise ValueError("Not a valid SPARSE6 format.")
    s6 = s6[1:]  # Remove the leading ':'
    s6_bytes = [ord(c) - 63 for c in s6]
    # Decode the number of vertices n
    index = 0
    if s6_bytes[0] <= 62:
        n = s6_bytes[0]
        index = 1
    else:
        n = ((s6_bytes[1] << 12) |
             (s6_bytes[2] << 6) |
             s6_bytes[3])
        index = 4
    k = (n - 1).bit_length()  # Number of bits needed to represent n-1 in binary
    G = nx.Graph()
    G.add_nodes_from(range(n))
    bitstream = []
    for byte in s6_bytes[index:]:
        bits = f"{byte:06b}"
        bitstream.extend(bits)
    v = 0
    i = 0
    while i < len(bitstream):
        if len(bitstream) - i < k + 1:
            break
        b = int(bitstream[i])
        x = int("".join(bitstream[i + 1:i + 1 + k]), 2)
        i += k + 1
        if b == 1:
            v += 1
        if x > v:
            v = x
        else:
            G.add_edge(x, v)
    return G
    
#sparse6_graph = Sread_graph_from_file(filename, Lnum)
#graph = decode_sparse6(sparse6_graph)

#___________________________________________________________________

g = pickle.load(open('../RandomG/G_L=80.txt','rb'))
pos = pickle.load(open('../RandomG/Pos_L=80.txt','rb'))
PL_HC = pickle.load(open('../RandomG/Plaqs_L=80.txt','rb'))
PL_sq = pickle.load(open('../RandomG/PlaqSQ_L=80.txt','rb'))
print(len(PL_HC),len(PL_sq))

ME = pickle.load(open('../RandomG/FS/FS_T=0.02_V=-1.00000.txt','rb'))

def NoFP(pl,M,mm):
    c = 0
    for i in pl:
        mn = mm.copy()
        for a, b in i:
            if ((a, b) in M or (b, a) in M):
                mn.extend([a, b])
        if set(collections.Counter(mn).values()) == {1} and len(set(mn)) == len(i):
            c += 1
    return c

def check_tuples_in_lists(list_of_tuples, list_of_list_of_tuples):
    normalized_target = set(tuple(sorted(t)) for t in list_of_tuples)
    for sublist in list_of_list_of_tuples:
        normalized_sublist = set(tuple(sorted(t)) for t in sublist)
        if normalized_target.issubset(normalized_sublist):
            return True
    return False

def OrdP_MP(PL_HC):
    ts = 0
    fs = 0
    for i in PL_HC:
        if check_tuples_in_lists(i,PL_sq):
            ts += 1
            if NoFP([i],ME,[]) == 1:
                fs += 1
    return fs,ts

c = NoFP(PL_HC,ME,[])
print(c)
print(OrdP_MP(PL_HC))

for i in PL_HC:
    cycle_coords = [pos[v] for u,v in i]
    if check_tuples_in_lists(i,PL_sq):
        if NoFP([i],ME,[])==1:
            polygon = Polygon(cycle_coords, closed=True, edgecolor='black',facecolor='red',alpha=0.6)
        else:
            polygon = Polygon(cycle_coords, closed=True, edgecolor='black',facecolor='white',alpha=0.3)
    else:
        print(len(i))
        if NoFP([i],ME,[])==1:
            polygon = Polygon(cycle_coords, closed=True, edgecolor='black',facecolor='blue',alpha=0.6)
        else:
            polygon = Polygon(cycle_coords, closed=True, edgecolor='black',facecolor='white',alpha=0.3)
    plt.gca().add_patch(polygon)

nx.draw(g,pos,node_size=5,node_color='orange')
plt.show()


#Ordering of node's positions in polygons are the problem, fix it (clockwise or anti) for right face coloring !!!!
'''