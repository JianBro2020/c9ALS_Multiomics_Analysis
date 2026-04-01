import scanpy as sc

def normalize_counts(
    adata: sc.AnnData,
    input_layer: str = "raw_counts",
    output_layer: str = "normalized_counts",
    target_sum: float = 1e4,
    overwrite: bool = False
) -> sc.AnnData:
    """
    Perform library-size normalization so that total counts per cell become comparable.
    -----------
    Parameters
    -----------
        adata: AnnData object. Must contain adata.layers["raw_counts"].
        input_layer: Layer containing the raw count matrix to normalize.
        output_layer: Layer name in which normalized counts will be stored.
        target_sum: Target total counts per cell after normalization.
        overwrite: Whether to overwrite 'output_layer' if it already exists.
    --------
    Returns
    --------
        adata: AnnData object, updated with adata.layers[output_layer]
    """
    if input_layer not in adata.layers:
        raise KeyError(
            f"'{input_layer}' not found in adata.layers. "
            f"Please ensure the input count layer exists."
        )
    if output_layer in adata.layers and not overwrite:
        raise ValueError(
            f"adata.layer['{output_layer}'] already exists. "
            f"Set overwrite = True if you want to replace it."
        )

    adata_tmp = adata.copy()
    adata_tmp.X = adata_tmp.layers[input_layer].copy()

    sc.pp.normalize_total(
        adata_tmp,
        target_sum = target_sum,
        inplace = True
    )

    adata.layers[output_layer] = adata_tmp.X.copy()

    print(
        f"Normalization complete.\n"
        f"Input layer used: {input_layer}\n"
        f"Output stored in: adata.layers['{output_layer}']\n"
        f"Target sum per cell: {target_sum}"
    )

    return adata

def log_transform(
    adata: sc.AnnData,
    input_layer: str = "normalized_counts",
    output_layer: str = "log1p_normalized_counts",
    overwrite: bool = False,
) -> sc.AnnData:
    """
    Apply log1p transformation to a normalized count matrix.
    -----------
    Parameters
    -----------
        adata: AnnData object. Must contain adata.layers[input_layer].
        input_layer: Layer containing normalized counts.
        output_layer: Layer name in which log-transformed values will be stored. 
        overwrite: Whether to overwrite 'output_layer' if it already exists.
    --------
    Returns
    --------
        adata: AnnData object, updated with adata.layers[output_layer]
    """

    if input_layer not in adata.layers:
        raise KeyError(
            f"'{input_layer}' not found in adata.layers. "
            f"Please run normalize_counts() first or specify a valid input_layer."
        )    
    
    if output_layer in adata.layers and not overwrite:
        raise ValueError(
            f"adata.layers['{output_layer}'] already exists. "
            f"Set overwrite = True if you want to replace it."
        )

    adata_tmp = adata.copy()
    adata_tmp.X = adata_tmp.layers[input_layer].copy()
    
    sc.pp.log1p(adata_tmp)

    adata.layers[output_layer] = adata_tmp.X.copy()

    print(
        f"Log transformation complete.\n"
        f"Input layer used: {input_layer}\n"
        f"Output stored in: adata.layers['{output_layer}']"
    )

    return adata

def select_hvgs(
    adata: sc.AnnData,
    input_layer: str = "log1p_normalized_counts",
    n_top_genes: int = 2000,
    flavor: str = "seurat",
) -> sc.AnnData:
    """
    Identify highly variable genes (HVGs).
    -----------
    Parameters
    -----------
        adata: AnnData object. Must contain adata.layers[input_layer].
        input_layer: Layer containing expression values used for HVG selection.
        n_top_genes: Number of top highly variable genes to retain.
        flavor: Method for HVG selection.
    --------
    Returns
    --------
        adata: AnnData object, updated with HVG annotations in adata.var["highly_variable"]
    """
    if input_layer not in adata.layers:
        raise KeyError(
            f"'{input_layer}' not found in adata.layers. "
            f"Please run log_transform() first or specify a valid input_layer.")
    
    adata_tmp = adata.copy()
    adata_tmp.X = adata_tmp.layers[input_layer].copy()

    sc.pp.highly_variable_genes(
        adata_tmp,
        n_top_genes = n_top_genes,
        flavor = flavor,
        inplace = True
    )

    adata.var["highly_variable"] = adata_tmp.var["highly_variable"].values
    n_hvgs = int(adata.var["highly_variable"].sum())

    print(
        f"HVG selection complete.\n"
        f"Input layer used: {input_layer}\n"
        f"Number of HVGs selected: {n_hvgs}\n"
    )

    return adata

def scale_expression(
    adata: sc.AnnData,
    input_layer: str = "log1p_normalized_counts",
    use_hvgs_only: bool = True,
    max_value: float = 10.0
) -> sc.AnnData: 
    """
    Scale expression values so that each gene has mean 0 and variance 1. 
    This function creates a new AnnData object for downstream dimensionality reduction. By default, only high variable genes are retained.
    -----------
    Parameters
    -----------
        adata: AnnData object. Must contain adata.layers[input_layer].
        input_layer: Layer containing the expression matrix to scale.
        use_hvgs_only: Whether to subset to 'adata.var["highly_variable"]' before scaling.
        max_value: Clip scaled values to this maximum absolute value.
    --------
    Returns
    --------
        adata_scaled: New AnnData object where X contains the scaled matrix. 
                      If 'use_hvg_only = True', only HVGs are retained.
    """

    if input_layer not in adata.layers:
        raise KeyError(
            f"'{input_layer}' not found in adata.layers. "
            f"Please run log_transform() first or specify a valid input_layer."
        )

    adata_scaled = adata.copy()
    adata_scaled.X = adata_scaled.layers[input_layer].copy()

    if use_hvgs_only:
        if "highly_variable" not in adata_scaled.var.columns:
            raise KeyError(
                "'highly_variable' not found in adata.var. "
                "Please run select_highly_variable_genes() first."
            )
        adata_scaled = adata_scaled[:, adata_scaled.var["highly_variable"]].copy()
    
    sc.pp.scale(
        adata_scaled,
        max_value = max_value,
        zero_center = True
    )

    print(
        f"Expression scaling complete.\n"
        f"Input layer used: {input_layer}\n"
        f"Using HVGs only: {use_hvgs_only}\n"
        f"Number of genes retained: {adata_scaled.n_vars}\n"
        f"Scaled matrix stored in adata_scaled.X"
    )

    return adata_scaled
