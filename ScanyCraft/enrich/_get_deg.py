import numpy as np
import scanpy as sc
from tqdm.auto import tqdm


def get_pairwise_deg(adata, comparisons, groupby_col='group', method='wilcoxon'):
    """
    Perform pairwise differential expression analysis for specified comparisons.
    
    Parameters:
    - adata: AnnData object
    - comparisons: list of tuples, each containing two groups to compare
    - groupby_col: str, name of the column in adata.obs to use for grouping
    - method: str, method to use for differential expression analysis
    
    Returns:
    - Dictionary containing results of all specified comparisons

    Example:
    --------
    # Assuming 'adata' is your AnnData object and 'groups' is your list of groups
    groups = ['O', 'M', 'Y']
    comparisons = list(combinations(groups, 2))

    # Filter the AnnData object for a specific condition
    adata_filtered = adata[adata.obs["pca-leiden"] == "1"]

    # Perform pairwise differential expression analysis
    result_dict = get_pairwise_deg(adata_filtered, comparisons, groupby_col='age')

    # Access results for a specific comparison
    o_vs_m_df = result_dict['O_vs_M']



    comparisons = [("O", "Y")]
    epi_clusters = ["2", "3", "5", "6", "15"]
    adata_clusters_deg = dict.fromkeys(epi_clusters) 
    {'2': None, '3': None, '5': None, '6': None, '15': None}

    for cluster in adata_clusters_deg:
        print(f"The current cluster is {cluster}")
        adata_sub = age_adata[age_adata.obs["pca-leiden"] == cluster]
    
        adata_clusters_deg[cluster] = get_pairwise_deg(adata_sub,
                                                      comparisons, 
                                                      groupby_col='age')
        

  
    """
    
    # Check and set log1p base if not present
    if 'log1p' not in adata.uns or 'base' not in adata.uns['log1p']:
        print("Warning: log1p base not found in adata.uns. Setting default base e.")
        adata.uns['log1p'] = {'base': np.e}
    
    results_dict = {}
    
    # Create a progress bar
    pbar = tqdm(total=len(comparisons), desc="Processing comparisons")
    
    for group1, group2 in comparisons:
        key = f'{group1}_vs_{group2}'
        
        # Update progress bar description with current comparison
        pbar.set_description(f"Processing {key}")
        
        sc.tl.rank_genes_groups(adata, groupby=groupby_col, groups=[group1], 
                                reference=group2, method=method, key_added=key, pts=True, use_raw=True)
        
        df = sc.get.rank_genes_groups_df(adata, group=group1, key=key)
        
        df['comparison'] = key
        results_dict[key] = df
        
        # Update progress bar
        pbar.update(1)
    
    pbar.close()
    
    return results_dict