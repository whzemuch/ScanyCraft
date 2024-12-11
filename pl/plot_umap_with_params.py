import matplotlib.pyplot as plt
from itertools import product
import scanpy as sc

def validate_params(min_dists, spreads):
    """Ensure parameters are lists and contain valid entries."""
    if not all(isinstance(x, (int, float)) for x in min_dists + spreads):
        raise ValueError("All entries in min_dists and spreads must be integers or floats.")
    if not min_dists or not spreads:
        raise ValueError("min_dists and spreads must be non-empty lists.")
    
    # Generate all combinations of parameters for easy iteration later
    param_combinations = list(product(enumerate(min_dists), enumerate(spreads)))
    return min_dists, spreads, param_combinations

def compute_and_plot_umap(adata, min_dist, spread, ax):
    """Compute UMAP with specified parameters and plot on the given axes."""
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)
    sc.pl.umap(
        adata,
        color=["louvain"],
        title=f"min_dist = {min_dist}, spread = {spread}",
        s=40,
        ax=ax,
        show=False
    )

def setup_plot_grid(min_dists, spreads):
    """Set up a figure with a grid of subplots."""
    fig, axes = plt.subplots(len(min_dists), len(spreads),
                             figsize=(len(spreads) * 3 + 2, len(min_dists) * 3))
    return fig, axes

def plot_umap_grid(adata, min_dists=[0.1, 1, 2], spreads=[0.5, 1, 5]):
    """
    Plot a grid of UMAP visualizations with varying min_dist and spread parameters.

    Parameters:
        adata (AnnData): An AnnData object containing the data.
        min_dists (list of float): Minimum distance values for UMAP.
        spreads (list of float): Spread values for UMAP.
    """
    min_dists, spreads, param_combinations = validate_params(min_dists, spreads)
    adata_temp = adata.copy()
    fig, axes = setup_plot_grid(min_dists, spreads)

    for (i, min_dist), (j, spread) in param_combinations:
        ax = axes[i][j] if axes.ndim > 1 else axes[max(i, j)]
        compute_and_plot_umap(adata_temp, min_dist, spread, ax)

    plt.tight_layout()
    plt.show()
    plt.close(fig)
    del adata_temp

# Example of usage
# plot_umap_grid(adata)
