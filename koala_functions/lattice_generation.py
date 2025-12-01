import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from koala import graph_utils as gu
from koala import plotting as pl
from dimer_models.lattice_generation import expand_edges

import pickle

mpl.rcParams.update(
    {
        "font.size": 11,
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern"],
    }
)


def load_lattice(directory, lattice_index):
    # pull the lattice_index'th lattice from the file at the directory

    with open(f"{directory}", "rb") as f:
        try:
            _ = [pickle.load(f) for _ in range(lattice_index - 1)]
            lattice_original, expandable_edges = pickle.load(f)
        except EOFError as exc:
            raise ValueError(
                f"The file does not contain {lattice_index} lattices, reduce lattice_index"
            ) from exc
    return lattice_original, expandable_edges


def reduce_lattice(lattice_original, expandable_edges, proportion):

    subset = np.random.choice(
        expandable_edges,
        np.round(proportion * len(expandable_edges)).astype(int),
        replace=False,
    )
    return expand_edges(lattice_original, subset)


if __name__ == "__main__":

    # directory for the lattice, have a look at the files I added
    directory = "koala_functions/lattices/00400.pkl"

    # Each file contains 10 lattices so this number can be 1-10
    lattice_index = 10

    # proportions of squares to remove
    rem_square_prop = 1

    # load the gem lattice
    lattice_original, expandable_edges = load_lattice(directory, lattice_index)

    # remove rem_square_prop of the squares
    reduced_lattice = reduce_lattice(lattice_original, expandable_edges, rem_square_prop)

    # clean that shiz up
    reduced_lattice = gu.com_relaxation(reduced_lattice)

    fig, ax = plt.subplots(1, 2, figsize = (10,5))
    for a in ax:
        a.set_aspect("equal")

    orig_squares = [p.n_sides == 4 for p in lattice_original.plaquettes]
    reduced_squares = [p.n_sides == 4 for p in reduced_lattice.plaquettes]

    pl.plot_edges(lattice_original, ax=ax[0])
    pl.plot_edges(reduced_lattice, ax=ax[1])

    pl.plot_plaquettes(lattice_original, orig_squares, color_scheme=['w','g'], ax = ax[0])
    pl.plot_plaquettes(reduced_lattice, reduced_squares, color_scheme=['w','g'], ax = ax[1])

    plt.show()
