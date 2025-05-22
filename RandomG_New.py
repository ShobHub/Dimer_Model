import collections
from scipy.spatial import KDTree
import pickle
import networkx as nx
import matplotlib.pyplot as plt
from shapely.geometry import LineString
from math import atan2, degrees
import numpy as np
import random

from NX_to_koala import nx_to_koala
from koala_to_NX import koala_to_nx
from koala import plotting
from koala.graph_utils import clockwise_about
from koala import graph_utils as gu
from koala.functions.koala_plantri import tutte_embedding
from Dimerisations import count_dimers,plot_dimers
from koala.lattice import Lattice

def generate_polar_random_walk(n_steps):
    angles = np.random.uniform(0, 2 * np.pi, n_steps)
    step_lengths = np.random.uniform(150, 155, n_steps)  #np.ones(n_steps)

    x_steps = step_lengths * np.cos(angles)
    y_steps = step_lengths * np.sin(angles)

    x_positions = np.cumsum(x_steps)
    y_positions = np.cumsum(y_steps)
    x_positions = np.append(x_positions, x_positions[0])
    y_positions = np.append(y_positions, y_positions[0])

    return x_positions, y_positions

def find_self_intersections(x_positions, y_positions):
    intersections = {}
    num_points = len(x_positions)

    start_point = (round(x_positions[0], 5), round(y_positions[0], 5))
    end_point = (round(x_positions[-1], 5), round(y_positions[-1], 5))

    for i in range(num_points - 2):
        segment1 = LineString([(x_positions[i], y_positions[i]), (x_positions[i + 1], y_positions[i + 1])])
        for j in range(i + 2, num_points - 1):
            segment2 = LineString([(x_positions[j], y_positions[j]), (x_positions[j + 1], y_positions[j + 1])])
            if segment1.intersects(segment2):
                intersection_point = segment1.intersection(segment2)
                if intersection_point.geom_type == 'Point':
                    point = (round(intersection_point.x, 5), round(intersection_point.y, 5))
                    # Exclude the start and end points
                    if point != start_point and point != end_point:
                        if point not in intersections:
                            intersections[point] = []
                        intersections[point].append((i, j))
    return intersections

def create_multigraph(intersections, x_positions, y_positions):
    G = nx.MultiGraph()  # Create a MultiGraph
    labels = {i: point for i, point in enumerate(intersections.keys())}  # Label intersections with integers
    pos = labels.copy()
    G.add_nodes_from(labels.keys())

    # Create edges based on the random walk path
    prev_label = None  # To track the previous intersection point's label
    first_label = None  # To track the first intersection label encountered
    last_label = None  # To track the last intersection label encountered
    num_points = len(x_positions)

    # Traverse through each segment of the random walk
    for i in range(num_points - 1):
        current_segment = LineString([(x_positions[i], y_positions[i]),(x_positions[i + 1], y_positions[i + 1])])
        current_intersections = []

        # Check for intersections with existing intersection points
        for point in intersections.keys():
            # Check if the current segment intersects the point
            if current_segment.distance(LineString([point, point])) < 1e-5:  # Using distance to check proximity
                #current_intersections.append(point)
                distance = np.sqrt((point[0] - x_positions[i]) ** 2 + (point[1] - y_positions[i]) ** 2)
                current_intersections.append((distance, point))

        # Sort intersections by distance from the current position
        current_intersections.sort(key=lambda x: x[0])
        current_intersections = [p[1] for p in current_intersections]

        # If any intersections were found, create edges
        for current_point in current_intersections:
            current_label = [k for k, v in labels.items() if v == current_point][0]
            # Set the first_label if it's the first intersection encountered
            if first_label is None:
                first_label = current_label

            # Update last_label every time we encounter an intersection
            last_label = current_label
            # Create an edge from prev_label to current_label if prev_label exists
            if prev_label is not None:
                G.add_edge(prev_label, current_label)  # Create an edge from prev to current
            # Update prev_label to the current label
            prev_label = current_label

    # After traversal, add an edge between the first and last intersection points if they exist
    if first_label is not None and last_label is not None:
        G.add_edge(last_label,first_label)

    return G, pos

def remove_self_loops_until_none(G):
    # Continue to process as long as there are self-loops in the graph
    while True:
        # Find all nodes with self-loops
        nodes_with_self_loops = [node for node in G.nodes if G.has_edge(node, node)]
        # If no self-loops are found, break the loop
        if not nodes_with_self_loops:
            break
        # Process each node with a self-loop
        for node in nodes_with_self_loops:
            # Get all neighbors of the node (excluding itself)
            neighbors = list(G.neighbors(node))
            neighbors = [n for n in neighbors if n != node]
            # Add edges between all pairs of neighbors
            if len(neighbors) == 2:
                G.add_edge(neighbors[0], neighbors[1])
            else:
                G.add_edge(neighbors[0], neighbors[0])
            # Remove the node with the self-loop
            G.remove_node(node)
    return G

def transform_to_3rpb(G,pos,offset,offset1):
    G_new = nx.Graph()
    pos_new = {}
    E_new = []
    nodes = 0

    pos = {i:np.array(j) for i,j in pos.items()}
    Ecount = collections.Counter(G.edges())
    for e,c in Ecount.items():
        a,b = e
        if c == 1:
            dir = (pos[a] - pos[b])/np.linalg.norm(pos[a] - pos[b])
            pos_new[nodes] = pos[a] - offset * dir
            pos_new[nodes+1] = pos[b] + offset * dir
            E_new.append((nodes,nodes+1))
            nodes += 2
        elif c == 2:
            dir = (pos[a] - pos[b]) / np.linalg.norm(pos[a] - pos[b])
            dir_perp = np.array([dir[1],-dir[0]])
            #print(np.dot(dir,dir_perp),np.dot(dir,-dir_perp))
            pos_new[nodes] = pos[a] - offset * dir - offset1 * dir_perp
            pos_new[nodes + 1] = pos[b] + offset * dir - offset1 * dir_perp
            pos_new[nodes + 2] = pos[a] - offset * dir + offset1 * dir_perp
            pos_new[nodes + 3] = pos[b] + offset * dir + offset1 * dir_perp
            E_new.extend([(nodes, nodes + 1),(nodes + 2, nodes + 3)])
            nodes += 4
        else:
            print(c,'*')

    L = list(pos_new.values())
    kdtree = KDTree(L)
    for v in G.nodes():
        query_point = np.array(pos[v])
        _, indices = kdtree.query(query_point,k=4)
        ord = {}
        for j in indices:
            angle = degrees(atan2(L[j][1] - pos[v][1], L[j][0] - pos[v][0]))
            ord[(angle + 360) % 360] = j
        ord_s = sorted(ord)
        for j in range(3):
            E_new.append((ord[ord_s[j]],ord[ord_s[j+1]]))
        E_new.append((ord[ord_s[-1]], ord[ord_s[0]]))

    G_new.add_edges_from(E_new)
    return G_new,pos_new

x_, y_ = generate_polar_random_walk(n_steps=40)    #300)
intersections = find_self_intersections(x_,y_)
x,y = zip(*intersections.keys())
plt.plot(x_,y_,linewidth = 0.8)
plt.scatter(x_,y_,s=10)
plt.scatter(x,y,c = 'r')
plt.savefig('RGWalk.pdf')
# plt.show()
sgffdhgf
G, pos = create_multigraph(intersections,x_,y_)
#nx.draw(G, pos, node_size=5, node_color = 'yellow',width = 2,alpha = 0.5,with_labels=True)

G = remove_self_loops_until_none(G)
#nx.draw(G, pos, node_size=5, node_color = 'r',edge_color = 'r',width = 2,alpha = 0.5,with_labels=True,font_size = 8)
Degrees = [d for _, d in G.degree()]
print("Are all vertices degree 4?", all(d == 4 for d in Degrees))

G, pos = transform_to_3rpb(G,pos,0.008,0.004)   # - 300 - 0.01,0.005
# G = make_planar_and_bipartite(G)
Degrees = [d for _, d in G.degree()]
print("Are all vertices degree 3?", all(d == 3 for d in Degrees))
print("Is the graph bipartite?", nx.is_bipartite(G))
print("Is the graph planar?", nx.check_planarity(G))
nx.draw(G, pos, node_size=5) #,with_labels=True,font_size = 8)
plt.show()

g = G.copy()
lattice = nx_to_koala(pos,g)
print(len(g.nodes()),len(g.edges()))
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
# plotting.plot_edges(lattice)
# plt.show()

pos,E,PL = koala_to_nx(lattice)
print(len(PL))

G = nx.Graph()
G.add_edges_from(E)
nx.draw(G,pos,node_size=0) #,with_labels=True)
plt.show()

pickle.dump(G,open('../RandomG/G_L=100_tuttefsdfd.txt','wb'))
pickle.dump(pos,open('../RandomG/Pos_L=100_tuttedsfds.txt','wb'))
pickle.dump(PL,open('../RandomG/Plaqs_L=100_tuttesdfds.txt','wb'))



'''
def create_bipartite_planar_graph(n):
    # Create a small bipartite graph using the cycle graph as a basis
    # Cycle graph C(2n) is bipartite when n is even
    B = nx.cycle_graph(2 * n)

    # Relabel nodes to form two sets
    left_nodes = {i for i in range(2 * n) if i % 2 == 0}
    right_nodes = set(range(2 * n)) - left_nodes

    # Add additional edges while maintaining bipartiteness and planarity
    for u in left_nodes:
        for v in right_nodes:
            if B.degree[u] < 3 and B.degree[v] < 3 and not B.has_edge(u, v):
                B.add_edge(u, v)
                if not nx.check_planarity(B)[0]:
                    # If adding the edge breaks planarity, remove it
                    B.remove_edge(u, v)
    return B

# Create a 3-regular bipartite planar graph
n = 100  # Number of nodes per set (total nodes = 2 * n)
B = create_bipartite_planar_graph(n)
'''