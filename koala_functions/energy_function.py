import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

from koala import graph_utils as gu
from koala import plotting as pl
from koala.lattice import Lattice
from dimer_models.lattice_generation import (
    bipartite_squarefull,
)


mpl.rcParams.update(
    {
        "font.size": 11,
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern"],
    }
)


def _lattice_cross_pairs(lattice: Lattice):
    # generates a list of all pairs of edges that are connected
    # diagonally across another edge

    # this creates a list of every crossed pair of bonds around each bond
    all_terms = np.zeros([lattice.n_edges, 2, 2], dtype=int)
    for n, edge in enumerate(lattice.edges.indices):

        # fint the three edges clockwise around each of the two vertices
        adj_e_1 = lattice.vertices.adjacent_edges[edge[0]]
        adj_e_2 = lattice.vertices.adjacent_edges[edge[1]]

        # roll them until they both start with edge
        a1 = np.roll(adj_e_1, -np.where(adj_e_1 == n)[0])
        a2 = np.roll(adj_e_2, -np.where(adj_e_2 == n)[0])

        # the pairs are now given by
        p1 = np.array([a1[1], a2[1]])
        p2 = np.array([a1[2], a2[2]])

        all_terms[n, 0] = p1
        all_terms[n, 1] = p2
    return all_terms


def energy_of_dimerisation(
    lattice: Lattice, dimerisation: np.ndarray, cross_pairs=None
):
    # calculates the classical energy of a consifuration, where a cross-pair
    # with both being dimers costs energy and anything else does not

    # allows you to calculate this thing once, rather than
    # recalculating the list for every energy call
    if cross_pairs is None:
        cross_pairs = _lattice_cross_pairs(lattice)

    energies_per_bond = np.prod(dimerisation[cross_pairs], axis=2)
    energies_per_bond = np.sum(energies_per_bond, axis=1)

    return energies_per_bond


def main():

    # first we generate a random bipartite lattice
    lattice = bipartite_squarefull(100)

    # next we generate a random dimerisation - if you
    # choose the [0] dimerisation it generates the flat configuration
    dimerisation = gu.dimerise(lattice, 40)[39]

    # find the energy of this dimerisation
    e_bond = energy_of_dimerisation(lattice, dimerisation)

    _, ax = plt.subplots(figsize=(5, 5))

    # plot non dimer edges in black
    pl.plot_edges(
        lattice,
        subset=np.where(dimerisation == 0),
        ax=ax,
        color_scheme=["k"],
        linewidth=2,
    )

    # plot_dimers in green
    pl.plot_edges(
        lattice,
        subset=np.where(dimerisation == 1),
        ax=ax,
        color_scheme=["g"],
        linewidth=5,
    )

    # plot the bonds with high energy in blue
    pl.plot_edges(
        lattice,
        ax=ax,
        subset=np.where(e_bond)[0],
        zorder=10,
        color_scheme="b",
        linewidth=5,
    )
    ax.set_aspect("equal")
    plt.show()


if __name__ == "__main__":
    main()
