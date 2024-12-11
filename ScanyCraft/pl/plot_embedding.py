import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


def plot_embeddings_with_labels(adata, basis, color_attribute, overlay_attribute, group_col=None,
                                subplot_width=5, text_props_kw=None, color_map=None, suptitle=None):
    """
    Plot UMAP embeddings with labels for an AnnData object. This function allows
    coloring by a given attribute and annotates clusters with labels indicating
    their median positions. If a grouping column is provided, separate subplots for
    each group will be created; otherwise, a single plot for the entire dataset is generated.

    Parameters:
    adata : anndata.AnnData
        An AnnData object containing the single-cell dataset to be plotted.
    basis : str
        The key in `adata.obsm` that corresponds to the basis for the embedding (e.g., 'X_umap', 'X_pca').
    color_attribute : str
        The name of the column in `adata.obs` that holds the categorical variable used to color the points.
    overlay_attribute : str
        The name of the column in `adata.obs` to use for labeling points with their median positions.
    group_col : str, optional
        The name of the column in `adata.obs` to use for creating subplots for each unique category.
        If None, a single plot for the entire dataset is generated (default: None).
    subplot_width : float, optional
        The base width of each subplot in inches (default: 5). The total width of the figure is automatically
        adjusted based on the number of groups.
    text_props_kw : dict, optional
        A dictionary of keyword arguments for configuring the appearance of the text labels (e.g., fontsize, color).
        See matplotlib's `text` documentation for valid options (default: None).
    color_map : dict, optional
    A dictionary mapping the `color_attribute` values to colors.

    Returns:
    None: The function directly renders the plot.

    Usage Example:
    ------------------------------
    text_properties = {
        'fontsize': 8,
        'color': 'black',
        'weight': 'bold'
    }
    plot_embeddings_with_labels(
        adata=adata,
        basis='X_umap',
        color_attribute='cell_type',
        overlay_attribute='leiden',
        group_col='batch',
        subplot_width=5,
        text_props_kw=text_properties
    )
    
    """
    
    # If group_col is provided and it's in adata.obs, plot for each group
    if group_col and group_col in adata.obs:
        groups = adata.obs[group_col].unique()
        fig, axes = plt.subplots(1, len(groups), figsize=(np.ceil(subplot_width * (len(groups) + 0.5)) , subplot_width), squeeze=False)
        axes = axes.flatten()
        for ax, group in zip(axes, groups):
            subset = adata[adata.obs[group_col] == group]
            _plot_one_embedding_with_label(subset, basis, color_attribute, overlay_attribute, ax, text_props_kw)
            ax.set_title(f'{group_col} = {group}')
    else:
        # Plot for the entire dataset
        fig, ax = plt.subplots(figsize=(subplot_width + 1, subplot_width))
        _plot_one_embedding_with_label(adata, basis, color_attribute, overlay_attribute, ax, text_props_kw, color_map)
        ax.set_title('All data')
    
    
    # If suptitle is provided, add it to the figure
    if suptitle:
        plt.gcf().suptitle(suptitle, fontsize=15)
        
    plt.tight_layout()
    plt.show()
    
    
    
def calculate_median_position(adata_subset, basis, overlay_attribute):
    """Calculate median positions of clusters in the embedding."""
    unique_labels = adata_subset.obs[overlay_attribute].unique()
    median_positions = {}
    for label in unique_labels:
        label_subset = adata_subset[adata_subset.obs[overlay_attribute] == label]
        median_positions[label] = np.median(label_subset.obsm[basis], axis=0)
    return median_positions

def add_text_labels(ax, median_positions, text_props_kw=None):
    """Add text labels at median positions on the axes object."""
    if text_props_kw is None:
        text_props_kw = {}
    for label, position in median_positions.items():
        ax.text(position[0], position[1], label, ha='center', va='center', **text_props_kw)
        
def _plot_one_embedding_with_label(adata, basis, color_attribute, overlay_attribute, 
                                   ax, text_props_kw=None, color_map=None):
    """Plot the embedding for a single subset or the entire dataset with labels."""
 
    sc.pl.embedding(adata, basis=basis, color=color_attribute, ax=ax, show=False)
    median_positions = calculate_median_position(adata, basis, overlay_attribute)
    add_text_labels(ax, median_positions, text_props_kw=text_props_kw)
    