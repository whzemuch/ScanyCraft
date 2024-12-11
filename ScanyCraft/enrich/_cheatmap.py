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
    index_col: str = "Term"          # Column name for the DataFrame index
    columns_col: str = "pca-leiden"  # Column name for the DataFrame columns
    nes_col: str = "NES"             # Column name for NES values
    pval_col: str = "FDR p-value"    # Column name for P-value values
    obs_col_name: str = "term"       # Name for observation metadata in AnnData
    var_col_name: str = "pca-leiden" # Name for variable metadata in AnnData
    thresholds_labels: list = None   # List of tuples for P-value thresholds and labels



class GSEAData:
    def __init__(self, gsea_dict, config: GSEAConfig, column_rename=None):
        """
        Initialize with GSEA data dictionary, configuration, and optional column renaming mapping.

        Parameters:
        - gsea_dict: dict
            Dictionary containing GSEA results.
        - config: GSEAConfig
            Configuration for column names and metadata.
        - column_rename: dict, optional, 
            {'2': 'c3', '4': 'c2', '11': 'c15', '8': 'c5', '7': 'c6'}
            Mapping for renaming columns.
        """
        self.gsea_dict = gsea_dict
        self.config = config
        self.column_rename = column_rename
        self.adata = None

    def prepare_data(self):
        """
        Prepare NES, P-value, and significance annotation DataFrames from the GSEA dictionary.
        """
        nes_df = (pd.concat(self.gsea_dict.values())
                  .reset_index()
                  .pivot_table(index=self.config.index_col, 
                               columns=self.config.columns_col, 
                               values=self.config.nes_col))

        pval_df = (pd.concat(self.gsea_dict.values())
                   .reset_index()
                   .pivot_table(index=self.config.index_col, 
                                columns=self.config.columns_col, 
                                values=self.config.pval_col))

        # Apply column renaming if column_rename is provided
        if self.column_rename:
            nes_df.rename(columns=self.column_rename, inplace=True)
            pval_df.rename(columns=self.column_rename, inplace=True)

        # Sort columns for consistent order
        nes_df.sort_index(axis=1, inplace=True)
        pval_df.sort_index(axis=1, inplace=True)

        # Create significance annotations using internal method
        sig_anno = self.create_significance_annotations(pval_df)

        return nes_df, pval_df, sig_anno

    def create_significance_annotations(self, pval_df):
        """
        Create significance annotations based on P-values with customizable thresholds and labels.

        Parameters:
        - pval_df: pd.DataFrame
            DataFrame containing P-values.

        Returns:
        - pd.DataFrame
            DataFrame containing significance annotations.
        """
        thresholds_labels = self.config.thresholds_labels or [(0.01, '**'), (0.05, '*')]

        # Separate thresholds and labels
        thresholds, labels = zip(*thresholds_labels)

        # Create conditions dynamically based on thresholds
        conditions = [pval_df < threshold for threshold in thresholds]

        # Generate a DataFrame with annotations based on the conditions
        annotations = pd.DataFrame(
            np.select(conditions, labels, default=''),
            index=pval_df.index,
            columns=pval_df.columns
        )

        return annotations

    def create_ann_data(self):
        """
        Create an AnnData object and populate it with NES, P-value, and significance annotation data.
        """
        nes_df, pval_df, sig_anno = self.prepare_data()

        self.adata = sc.AnnData(
            X=nes_df.to_numpy(),
            obs=pd.DataFrame({self.config.obs_col_name: nes_df.index}, index=nes_df.index),
            var=pd.DataFrame({self.config.var_col_name: nes_df.columns}, index=nes_df.columns),
            dtype='float64'
        )
        self.adata.layers['pval'] = pval_df.to_numpy()
        self.adata.layers['sig_anno'] = sig_anno.to_numpy()

        return self

def analyze_and_plot_gsea_results(
    df_gsea_dict, 
    config=None,
    annotation_name="Clusters",
    title="GSEA Analysis",
    figsize=(4,4),
    **kwargs
):
    """
    Analyze GSEA results and plot a heatmap with annotations.

    Parameters:
    - df_gsea_dict: dict
        Dictionary containing GSEA results for multiple clusters or conditions.
    - config: GSEAConfig, optional
        Configuration object with column names and thresholds for significance annotations.
    - annotation_name: str, optional
        Name of the top annotation to display on the heatmap.
    - title: str, optional
        Title of the heatmap plot.
    - figsize: tuple, optional (default: (6, 6))
        Size of the heatmap figure.
    - **kwargs: Additional keyword arguments for ClusterMapPlotter.

    Example Usage:
    --------------
    >>> df_gsea_dict = {
    >>>     "cluster_1": pd.DataFrame({
    >>>         "Term": ["Pathway1", "Pathway2", "Pathway3"],
    >>>         "pca-leiden": [1, 1, 2],
    >>>         "NES": [2.5, -1.8, 0.5],
    >>>         "FDR p-value": [0.0005, 0.02, 0.1]
    >>>     }),
    >>>     "cluster_2": pd.DataFrame({
    >>>         "Term": ["Pathway1", "Pathway4", "Pathway5"],
    >>>         "pca-leiden": [1, 3, 3],
    >>>         "NES": [-2.3, 1.7, -0.4],
    >>>         "FDR p-value": [0.0008, 0.015, 0.2]
    >>>     })
    >>> }
    >>> config = GSEAConfig(
    >>>     thresholds_labels=[(0.001, '***'), (0.01, '**'), (0.05, '*')]
    >>> )
    >>> analyze_and_plot_gsea_results(
    >>>     df_gsea_dict,
    >>>     config=config,
    >>>     title="My GSEA Heatmap"
    >>> )

    Returns:
    - None
        The function generates and displays the heatmap.
    """
    if config is None:
        raise ValueError("A GSEAConfig object must be provided.")

    # Step 1: Prepare GSEA Data
    gsea_data = GSEAData(df_gsea_dict, config=config)
    gsea_data.create_ann_data()

    # Step 2: Create Top Annotation
    annotation_df = pd.DataFrame(
        {annotation_name: gsea_data.adata.var.index}, 
        index=gsea_data.adata.var.index
    )
    top_annotation = pch.HeatmapAnnotation(df=annotation_df, plot=False)

    # Step 3: Plot Heatmap
    plt.figure(figsize=figsize)
    pch.ClusterMapPlotter(
        data=gsea_data.adata.to_df(),
        top_annotation=top_annotation,
        annot=gsea_data.adata.to_df(layer='sig_anno'),
        fmt=None,
        **kwargs
    )
    plt.suptitle(title)
    plt.show()



