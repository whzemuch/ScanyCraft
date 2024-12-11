import numpy as np
import pandas as pd
import scanpy as sc
import PyComplexHeatmap  as pch
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class GSEAConfig:
    """
    Configuration class for GSEA data processing and plotting.
    """
    index_col: str = "Term"         # Column name for the DataFrame index
    columns_col: str = "pca-leiden" # Column name for the DataFrame columns
    nes_col: str = "NES"            # Column name for NES values
    pval_col: str = "FDR p-value"   # Column name for P-value values
    obs_col_name: str = "term"      # Name for observation metadata in AnnData
    var_col_name: str = "pca-leiden" # Name for variable metadata in AnnData


class GSEAData:
    def __init__(self, gsea_dict, config: GSEAConfig, column_rename=None):
        """
        Initialize with GSEA data dictionary, configuration, and optional column renaming mapping.

        Parameters:
        - gsea_dict: dict, dictionary containing GSEA results.
        - config: GSEAConfig, configuration for column names and metadata.
        - column_rename: dict, optional, mapping to rename columns.
        """
        self.gsea_dict = gsea_dict
        self.config = config
        self.column_rename = column_rename
        self.nes_df = None
        self.pval_df = None
        self.adata = None
        self.annotations = None

    def prepare_data(self):
        """
        Prepare NES and P-value DataFrames from the GSEA dictionary.
        """
        self.nes_df = (pd.concat(self.gsea_dict.values())
                       .reset_index()
                       .pivot_table(index=self.config.index_col, 
                                    columns=self.config.columns_col, 
                                    values=self.config.nes_col))

        self.pval_df = (pd.concat(self.gsea_dict.values())
                        .reset_index()
                        .pivot_table(index=self.config.index_col, 
                                     columns=self.config.columns_col, 
                                     values=self.config.pval_col))

        # Apply column renaming if column_rename is provided
        if self.column_rename:
            self.nes_df.rename(columns=self.column_rename, inplace=True)
            self.pval_df.rename(columns=self.column_rename, inplace=True)

        # Sort columns for consistent order
        self.nes_df.sort_index(axis=1, inplace=True)
        self.pval_df.sort_index(axis=1, inplace=True)

        return self

    def create_ann_data(self):
        """
        Create an AnnData object with NES values and add P-value as a layer.
        """
        self.adata = sc.AnnData(
            X=self.nes_df.to_numpy(),
            obs=pd.DataFrame({self.config.obs_col_name: self.nes_df.index}, index=self.nes_df.index),
            var=pd.DataFrame({self.config.var_col_name: self.nes_df.columns}, index=self.nes_df.columns),
            dtype='float64'
        )
        self.adata.layers['pval'] = self.pval_df.to_numpy()
        return self

    def create_significance_annotations(self):
        """
        Create significance annotations based on P-values.
        """
        conditions = [self.pval_df < 0.01, self.pval_df < 0.05]
        choices = ['∗∗', '∗']
        self.annotations = pd.DataFrame(
            np.select(conditions, choices, default=''),
            index=self.pval_df.index, columns=self.pval_df.columns
        )



        return self




# Heatmap Plotting Class
class GSEAHeatmapPlotter:
    def __init__(self, adata, annotations):
        """
        Initialize with AnnData object and significance annotations.
        """
        self.adata = adata
        self.annotations = annotations
        self.top_annotation = None

    def create_top_annotation(self, annotation_name="cluster_id", custom_labels=None):
        """
        Create top annotation for the heatmap.

        Parameters:
        - annotation_name: str, the column name for the annotation in the top annotation DataFrame.
        - custom_labels: dict, optional, a mapping of column names to custom labels.

        Returns:
        - self: The updated object with `top_annotation`.
        """
        # Create the annotation DataFrame with dynamic column names
        annotation_df = pd.DataFrame(index=self.adata.var.index)

        if custom_labels:
            # Map custom labels to the annotation column
            annotation_df[annotation_name] = self.adata.var.index.map(custom_labels)
        else:
            # Use the default column names if no custom labels are provided
            annotation_df[annotation_name] = self.adata.var.index

        # Generate the top annotation object
        self.top_annotation = pch.HeatmapAnnotation(
            df=annotation_df,
            plot=False
        )
        return self

    def plot_heatmap(self, title, cmap='RdYlBu_r', vmin=-2, vmax=2, figsize=(4, 4), tight_layout=True, **kwargs):
        """
        Plot the heatmap with annotations and top annotation.

        Parameters:
        - title: str, the title of the heatmap.
        - cmap: str, the colormap to use.
        - vmin: float, minimum value for colormap scaling.
        - vmax: float, maximum value for colormap scaling.
        - figsize: tuple, size of the plot.
        - tight_layout: bool, whether to use tight layout for better spacing.
        - **kwargs: Additional keyword arguments for `pch.ClusterMapPlotter`.

        Example of Additional Arguments:
        - z_score: int, scale data (default: 0).
        - show_rownames: bool, display row names.
        - show_colnames: bool, display column names.
        - row_cluster: bool, cluster rows.
        - col_cluster: bool, cluster columns.
        - xticklabels_kws: dict, keyword arguments for x-axis tick labels.

        Returns:
        - None: Displays the heatmap.
        """
        plt.figure(figsize=figsize)
        pch.ClusterMapPlotter(
            data=self.adata.to_df(),
            top_annotation=self.top_annotation,
            annot=self.annotations,
            fmt=None,
            label="NES",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            **kwargs  # Pass all additional arguments to the ClusterMapPlotter
        )
        plt.suptitle(title)

        # Adjust layout if specified
        if tight_layout:
            plt.tight_layout()

        plt.show()




def analyze_and_plot_gsea_results(df_gsea_dict, index_col="Term", columns_col="pca-leiden", 
                                  nes_col="NES", pval_col="FDR p-value", 
                                  annotation_name="Aging Clusters", 
                                  title="GSEA Analysis of Aging Dataset", 
                                  cmap="coolwarm", vmin=-3, vmax=3, figsize=(6, 6)):
    """
    Prepare GSEA data and plot a heatmap with significance annotations.

    Parameters:
    - df_gsea_dict: dict
        Dictionary containing GSEA results for multiple clusters or conditions.
    - index_col: str, optional (default: "Term")
        Column name to use as the index for NES and P-value DataFrames.
    - columns_col: str, optional (default: "pca-leiden")
        Column name to use as the variable in NES and P-value DataFrames.
    - nes_col: str, optional (default: "NES")
        Column name containing NES (Normalized Enrichment Scores).
    - pval_col: str, optional (default: "FDR p-value")
        Column name containing P-values.
    - annotation_name: str, optional (default: "Aging Clusters")
        Name of the top annotation to display on the heatmap.
    - title: str, optional (default: "GSEA Analysis of Aging Dataset")
        Title of the heatmap plot.
    - cmap: str, optional (default: "coolwarm")
        Colormap to use for the heatmap.
    - vmin: float, optional (default: -3)
        Minimum value for colormap scaling.
    - vmax: float, optional (default: 3)
        Maximum value for colormap scaling.
    - figsize: tuple, optional (default: (6, 6))
        Size of the heatmap figure.

    Returns:
    - None
        The function generates and displays the heatmap.

    Example Usage:
    ---------------
    df_gsea_dict = {
        "cluster_1": pd.DataFrame(...),
        "cluster_2": pd.DataFrame(...),
    }
    analyze_and_plot_gsea_results(df_gsea_dict)
    """
    # Step 1: Prepare GSEA Data
    config = GSEAConfig(index_col=index_col, columns_col=columns_col, 
                        nes_col=nes_col, pval_col=pval_col)
    gsea_data = GSEAData(df_gsea_dict, config=config)
    gsea_data.prepare_data() \
             .create_ann_data() \
             .create_significance_annotations()

    # Step 2: Create Heatmap Plotter
    heatmap_plotter = GSEAHeatmapPlotter(gsea_data.adata, gsea_data.annotations)

    # Add top annotations
    heatmap_plotter.create_top_annotation(annotation_name=annotation_name)

    # Plot the heatmap
    heatmap_plotter.plot_heatmap(
        title=title,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        figsize=figsize
    )


