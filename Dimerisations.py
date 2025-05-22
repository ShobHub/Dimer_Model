from koala.lattice import Lattice
from koala import plotting as pl
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def count_dimers(dimers,lattice):
    plaquette_scores = np.zeros([dimers.shape[0], lattice.n_plaquettes])
    for i, dimerisation in enumerate(dimers):
        plaquette_scores[i] = np.array([np.sum(1 - 2 * dimerisation[p.edges]) for p in lattice.plaquettes])
    plaquette_scores = np.array(plaquette_scores) / 2

    cut_boundaries = True
    external_scores = lattice.n_plaquettes - 1 - np.sum(plaquette_scores, axis=1)
    external_scores = external_scores.astype(int)
    if not cut_boundaries:
        external_scores *= 0

    how_flippable = np.sum(plaquette_scores == 0, axis=1) + (external_scores == 0)
    # how_flippable = np.sum(plaquette_scores == 0, axis=1)

    most_flipped = np.max(how_flippable)
    most_flippable_indices = np.where(how_flippable == most_flipped)[0]
    most_flippable_dimerisations = dimers[most_flippable_indices]
    num_most = most_flippable_indices.shape[0]
    print(f"{num_most} dimerisations have {most_flipped}/{lattice.n_plaquettes} flippable")

    least_flipped = np.min(how_flippable)
    least_flippable_indices = np.where(how_flippable == least_flipped)[0]
    least_flippable_dimerisations = dimers[least_flippable_indices]
    num_least = least_flippable_indices.shape[0]
    print(
        f"{num_least} dimerisations have {least_flipped}/{lattice.n_plaquettes} flippable"
    )
    return plaquette_scores,external_scores,most_flippable_indices,least_flippable_indices

def plot_dimers(lattice,dimers,plaquette_scores,external_scores,most_flippable_indices,least_flippable_indices):
    cmap = mpl.colormaps["Blues"]
    colors = cmap(np.linspace(0, 1, np.max([plaquette_scores.astype("int")]) + 1))

    # Take colors at regular intervals spanning the colormap.
    # colors = cmap(np.linspace(0, 1, max(least_flippable_dimerisations.max(), least_flippable_dimerisations.max()) + 1))

    num_most = most_flippable_indices.shape[0]
    num_least = least_flippable_indices.shape[0]
    x = np.max([num_least, num_most]).astype(int)
    fig, ax = plt.subplots(x, 2, figsize=(10, 5 * x))
    # ax[0].set_title("Most flippable")
    # ax[1].set_title("Least flippable")

    for n, index in enumerate(most_flippable_indices[:0]):
        dim_to_plot = dimers[index]
        score_to_plot = plaquette_scores[index]
        external = external_scores[index]

        pl.plot_plaquettes(lattice, score_to_plot, colors, ax=ax[n, 0])
        pl.plot_edges(lattice, dim_to_plot, ["k", "yellow"], ax=ax[n, 0], linewidth=2)

        for i in range(lattice.n_plaquettes):
            x, y = lattice.plaquettes[i].center % 1
            ax[n, 0].annotate(
                score_to_plot[i].astype("int"), (x, y), ha="center", va="center"
            )
            ax[n, 0].set_title(f"Most flippable, external = {external}")

    for n, index in enumerate(least_flippable_indices[:1]):
        dim_to_plot = dimers[index]
        score_to_plot = plaquette_scores[index]
        external = external_scores[index]

        pl.plot_plaquettes(lattice, score_to_plot, colors, ax=ax[n, 1])
        pl.plot_edges(lattice, dim_to_plot, ["k", "yellow"], ax=ax[n, 1], linewidth=2)

        for i in range(lattice.n_plaquettes):
            x, y = lattice.plaquettes[i].center % 1
            ax[n, 1].annotate(
                score_to_plot[i].astype("int"), (x, y), ha="center", va="center"
            )
            ax[n, 1].set_title(f"Least flippable, external = {external}")

    for a in ax.flatten():
        a.axis("off")

    plt.show()


    return lattice