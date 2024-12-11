

def plot_violin_by_subset(adata, subset_by, gene_name, groupby, groupby_order=None,  ncol=1, figsize=(5, 12), y_range=None, layer=None, use_raw=True):
    """
    Plot violin plots for a specified gene across different subsets of data.

    Parameters:
    adata : AnnData
        Annotated data matrix.
    subset_by : dict
        Dictionary where the key is the column name to subset by and the value is a list of groups to subset.
    gene_name : str
        Name of the gene to plot.
    groupby : str
        Column name to group by for the violin plot.
    groupby_order : list of str, optional
        Order of categories for the x-axis (e.g., ["Y", "M", "O"]).
    ncol : int, optional (default=1)
        Number of plots in each row of the plot grid.
    figsize : tuple of int, optional (default=(5, 12))
        Size of the figure.
    y_range : tuple of float, optional
        Y-axis limits for all plots (y_min, y_max). If not provided, the range will be calculated from the data.

    Returns: None

    Example Usage:
    ---------------
    plot_violin_by_subset(age_adata, {"pca-leiden": ["3", "2", "5", "6", "15"]}, 
                          "Itga6", "age", ["Y", "M", "O"], ncol=2, y_range=(0, 1))


    """
    # Extract column name and groups from the dictionary
    column, groups = list(subset_by.items())[0]

    # Calculate the number of rows required
    nrow = (len(groups) + ncol - 1) // ncol

    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for ax, group in zip(axes, groups):
        adata_subset = adata[adata.obs[column] == group, :]
        sc.pl.violin(adata_subset, gene_name, groupby=groupby, layer=layer, use_raw=use_raw, order=groupby_order, ax=ax, show=False)
        if y_range is not None:
            ax.set_ylim(y_range)
        ax.set_title(f'{column} {group}')

    # Hide any unused axes
    for i in range(len(groups), len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()
    plt.show()

