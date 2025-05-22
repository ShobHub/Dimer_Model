from koala.lattice import Lattice
import numpy as np

def nx_to_koala(pos,G):
    x_coords = [i[0] for i in pos.values()]
    y_coords = [i[1] for i in pos.values()]

    min_x, max_x = min(x_coords), max(x_coords)
    min_y, max_y = min(y_coords), max(y_coords)

    range_x = max_x - min_x
    range_y = max_y - min_y

    Pos = np.zeros((len(G.nodes()),2))
    c = 0
    mapE = {}
    for i,j in pos.items():
        Pos[c][0] = (j[0] - min_x) / range_x
        Pos[c][1] = (j[1] - min_y) / range_y
        mapE[i] = c
        c += 1

    E = np.zeros((len(G.edges()),2),dtype = int)
    c=0
    for u,v in G.edges():
        E[c] = [mapE[u],mapE[v]]   #[int(u),int(v)]
        c += 1
    crossing = np.zeros((len(E),2))

    lattice = Lattice(Pos, E, crossing)

    return lattice