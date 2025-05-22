import numpy as np
import networkx as nx
from networkx import Graph, DiGraph, simple_cycles
import collections as cl
import matplotlib.pylab as plt
import pickle


# SQUARE-OCTAGON
def create_octagon_unit_cell():
	a = 1
	theta = np.linspace(np.pi / 8, 2 * np.pi + np.pi / 8, 9)[:-1]
	nodes = {}
	edges = []
	for i, angle in enumerate(theta):
		x = a * np.cos(angle)
		y = a * np.sin(angle)
		nodes[i] = (x, y)
		if i > 0:
			edges.append((i - 1, i))
	edges.append((7, 0))
	return nodes, edges

def create_square_octagon_lattice(rows, cols):
	G = nx.Graph()
	T = 1 + np.sqrt(2)
	unit_pos, unit_edges = create_octagon_unit_cell()
	pos = unit_pos.copy()
	G.add_edges_from(unit_edges)

	for i in range(1,cols):
		pos_c = {u+8*i: (v[0]+T*i , v[1]+0) for u,v in unit_pos.items()}
		edges_c = [(u+8*i,v+8*i) for u,v in unit_edges]
		G.add_edges_from(edges_c)
		G.add_edge(0+8*(i-1), 3+8*i)
		G.add_edge(7+8*(i-1), 4+8*i)
		pos.update(pos_c)

	unit_pos = pos.copy()
	unit_edges = list(G.edges()).copy()
	for i in range(1,rows):
		pos_c = {u + 8*cols*i: (v[0]+0, v[1]+T*i) for u,v in unit_pos.items()}
		edges_c = [(u + 8*cols*i, v + 8*cols*i) for u,v in unit_edges]
		G.add_edges_from(edges_c)
		for j in range(cols):
			G.add_edge(1+8*cols*(i-1)+8*j, 8*cols*i+6+8*j)
			G.add_edge(2+8*cols*(i-1)+8*j, 8*cols*i+5+8*j)
		pos.update(pos_c)
		# print(edges_c)

	for i in range(cols):                                  # |
		G.add_edge(5+8*i, 8*cols*(rows-1)+8*i+2)           # |
		G.add_edge(6+8*i, 8*cols*(rows-1)+8*i+1)           # |
															# PBCs
	for i in range(rows):                                  # |
		G.add_edge(3+8*cols*i, 8*(cols-1)+8*cols*i)        # |
		G.add_edge(4+8*cols*i, 8*cols-1+8*cols*i)          # |

	return G, pos

'''
# Parameters for the lattice
rows, cols = 12, 12

# Generate the lattice
G, pos = create_square_octagon_lattice(rows, cols)
Degrees = [d for _, d in G.degree()]
print("Are all vertices degree 3?", all(d == 3 for d in Degrees))
print(len(G.nodes()))

pickle.dump(G,open("../3DegLattices/SqOct12_PBC.txt",'wb'))
pickle.dump(pos,open("../3DegLattices/Pos_SqOct12_PBC.txt",'wb'))

# Plot the lattice
nx.draw(G, pos, node_size=5,with_labels=True)
plt.show()

BDL = []
a = [4,5,6,7]
b = [94,95,88,89]             #[62,63,56,57](8)      #[158,159,152,153] (20)
c = [1144,1145,1146,1147]     #[504,505,506,507]  #[3192,3193,3194,3195]
d = [1058,1059,1060,1061]     #[450,451,452,453]  #[3042,3043,3044,3045]
for i in range(cols):
	BDL = BDL+list(np.array(a) + 8*i)
BDL.append(88) #56) 152)
BDL.append(89) #57) 153)
for i in range(1,rows):
	BDL = BDL+list(np.array(b) + 8*cols*i)
BDL.append(1146)  #506) 3194)
BDL.append(1147)  #507) 3195)
for i in range(1,cols):
	BDL = BDL+list(np.array(c) - 8*i)
BDL.append(1060) #452) 3044)
BDL.append(1061) #453) 3045)
for i in range(1,rows):
	BDL = BDL+list(np.array(d) - 8*cols*i)

# print(BDL[:-2])
# pickle.dump(BDL[:-2],open("../3DegLattices/BDL_SqOct12_OBC.txt",'wb'))
# plt.show()
'''

# HONEYCOMB
rows = 30
cols = 30

G = nx.generators.lattice.hexagonal_lattice_graph(rows, cols, periodic=True)
G = nx.convert_node_labels_to_integers(G, first_label=0)
pos = nx.get_node_attributes(G, 'pos')
# print(len(G.nodes()),pos)

'''
pickle.dump(G,open("../3DegLattices/Hexa/Hexa30_PBC.txt",'wb'))
pickle.dump(pos,open("../3DegLattices/Hexa/Pos_Hexa30_PBC.txt",'wb'))

nx.draw(G, pos, node_size=5) #, with_labels=True)
plt.show()

BDL = list(range(rows*2+1))
a = [33,34]   #[121,122]
BDL = BDL + a
for i in range(1,cols-1):
	BDL = BDL + list(np.array(a)[::(-1)**i] + (rows*2+2)*i)
for i in range(rows*2+1):
	BDL = BDL + [159-i]   #1919
a = [126,125]   #[1798,1797]
BDL = BDL + a
for i in range(1,cols-1):
	BDL = BDL + list(np.array(a)[::(-1)**i] - (rows*2+2)*i)
print(BDL)
# pickle.dump(BDL,open("../3DegLattices/BDL_Hexa8_OBC.txt",'wb'))
# plt.show()
'''

#{p,q} HYPERBOLIC LATTICES

import matplotlib.cm as cmap
from hypertiling import HyperbolicGraph, GraphKernels
from hypertiling.graphics.plot import plot_tiling
from hypertiling.kernel.GRG_util import plot_graph
from hypertiling import HyperbolicTiling
p = 4
q = 6
nlayers = 5

# T = HyperbolicTiling(q,p,nlayers)
# print("Number of cells:", len(T))
# ffd
G = HyperbolicGraph(q, p, nlayers, kernel = "GRG")
nbrs = G.get_nbrs_list()  # get neighbors
crds = G.center_coords    # get coordinates of cell centers
p = G.p
G = plot_graph(nbrs,crds,p)
leaf_nodes = [node for node, degree in dict(G.degree()).items() if degree == 1]
G.remove_nodes_from(leaf_nodes)
G.remove_node(0)
pos = nx.get_node_attributes(G, 'pos')
print(len(G.nodes()))

print(nx.is_bipartite(G))
matching = nx.max_weight_matching(G) # Find maximum matching
print(len(matching) * 2 , len(G))

nx.draw(G,pos,node_size=0) #,with_labels = True,font_size=8)
# pickle.dump(G,open("../3DegLattices/(4,6)_OBC.txt",'wb'))
# pickle.dump(pos,open("../3DegLattices/Pos_(4,6)_OBC.txt",'wb'))

nodes = []
reverse_dict = {}
for u,v in pos.items():
	if G.degree(u)<6:
		nodes.append(v)
	reverse_dict[v] = u
cx, cy = (0,0)
def angle_from_center(point):
	x, y = point
	return np.arctan2(y - cy, x - cx)  # Angle in radians

sorted_nodes = sorted(nodes, key=angle_from_center, reverse=True)
BDL = []
for i in sorted_nodes:
	BDL.append(reverse_dict.get(i))
# print(BDL)
# pickle.dump(BDL,open("../3DegLattices/BDL_(4,6)_OBC.txt",'wb'))
plt.show()
