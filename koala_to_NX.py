import numpy as np

def koala_to_nx(lattice):
    pos = {}
    c = 0
    for i in lattice.vertices.positions:
        pos[c] = np.array(i)
        c += 1

    E = []
    for i in lattice.edges.indices:
        E.append(tuple(i))

    PL = []
    for i in lattice.plaquettes:
        v = i.vertices
        pl = []
        for j in range(len(v) - 1):
            pl.append((v[j], v[j + 1]))
        pl.append((v[-1], v[0]))
        PL.append(pl)

    return pos,E,PL