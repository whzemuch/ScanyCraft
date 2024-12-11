def find_common_items(*sets):
    """
    Find common elements across multiple sets.

    This function takes an arbitrary number of sets as input and returns
    a new set containing elements that are common to all input sets.

    Parameters:
    *sets : Set objects
        An arbitrary number of sets to find common elements from.

    Returns:
    set
        A set containing elements common to all input sets.
        Returns an empty set if no sets are provided.

    Examples:
    >>> set1 = {1, 2, 3, 4}
    >>> set2 = {2, 3, 4, 5}
    >>> set3 = {3, 4, 5, 6}
    >>> find_common_genes(set1, set2, set3)
    {3, 4}
    >>> find_common_genes(set1, set2, set3, {1, 3, 4, 7})
    {3, 4}
    >>> find_common_genes()
    set()
    """
    if not sets:
        return set()
    return set.intersection(*sets)


def create_msigdb_df(marker_genes, collection_name='custom_pathways', geneset_name='GENESET'):
    """
    Create a custom DataFrame in MSigDB format from a list of marker genes.

    Parameters:
    marker_genes (list): A list of gene symbols to be included in the custom geneset.
    collection_name (str, optional): The name of the collection. Defaults to 'custom_pathways'.
    geneset_name (str, optional): The name of the geneset. Defaults to 'SENNET_GENES'.

    Returns:
    pandas.DataFrame: A DataFrame with columns 'genesymbol', 'collection', and 'geneset'.

    Example:
    >>> genes = ['GENE1', 'GENE2', 'GENE3']
    >>> df = create_msigdb_custom_df(genes)
    >>> print(df)
      genesymbol       collection      geneset
    0      GENE1  custom_pathways  SENNET_GENES
    1      GENE2  custom_pathways  SENNET_GENES
    2      GENE3  custom_pathways  SENNET_GENES
    """
    return (pd.DataFrame(marker_genes, columns=["genesymbol"])
            .assign(collection=collection_name, 
                    geneset=geneset_name))




