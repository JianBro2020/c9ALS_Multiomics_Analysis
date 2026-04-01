import scanpy as sc
import pandas as pd
import soupx

def estimate_ambient_rna(
        adata: sc.AnnData,
        adata_raw: sc.AnnData,
        cluster_key: str = "soupx_clusters",
        target_sum: float = 1e4,
        leiden_resolution: float = 0.4,
        n_pcs: int = 30,
        n_neighbors: int = 15
) -> tuple[sc.AnnData, object]:
    """
    Estimate ambient RNA contamination using SoupX.
    -----------
    Parameters
    -----------
        adata: High-quality cells AnnData object after removing low-quality cells.
        adata_raw: AnnData object imported from ".../Solo.out/GeneFull/raw", containing all droplets (cells + empty droplets).
        cluster_key: Name of the obs column in which helper SoupX clusters will be stored. 
        target_sum: Target library size for temporary normalisation, to be used in helper clustering.
        leiden_resolution: Resolution for helper Leiden clustering.
        n_pcs: Number of Principal components for helper clustering.
        n_neighbors: Number of neighbors for helper clustering.
    --------
    Returns
    --------
        adata: AnnData object, updated with helper cluster labels in 'adata.obs[cluster_key]'.
        sc_chan: SoupX Channel object with estimated contamination parameters.
    """

    # -----------------------------------------
    # 1. Prepare the raw count matrix
    # -----------------------------------------
    # Align genes between raw droplets and filtered cells.
    adata_raw = adata_raw[:, adata.var_names].copy()
    # Extra raw integer counts for the retained cells.
    adata_filtered_raw = adata_raw[adata.obs_names, :].copy()

    # -----------------------------------------
    # 2. Create the helper clustering for SoupX
    # -----------------------------------------
    adata_tmp = adata.copy()
    if "raw_counts" in adata_tmp.layers:
        adata_tmp.X = adata_tmp.layers["raw_counts"].copy()
    
    sc.pp.normalize_total(adata_tmp, target_sum = target_sum)
    sc.pp.log1p(adata_tmp)
    sc.pp.pca(adata_tmp, n_comps = max(n_pcs, 30))
    sc.pp.neighbors(adata_tmp, n_neighbors = n_neighbors, n_pcs = n_pcs)
    sc.tl.leiden(adata_tmp, resolution = leiden_resolution, key_added = cluster_key)

    adata.obs[cluster_key] = adata_tmp.obs[cluster_key].copy()
    del adata_tmp

    # -----------------------------------------
    # 3. Build SoupX Channel
    # -----------------------------------------
    sc_chan = soupx.SoupChannel(
        tod = adata_raw.X.T.tocsr(),
        toc = adata_filtered_raw.X.T.tocsr(),
        metaData = pd.DataFrame(index=adata.obs_names),
    )
    sc_chan.setClusters(adata.obs[cluster_key].astype(str).values)

    # -----------------------------------------
    # 4. Estimate ambient contamination
    # -----------------------------------------
    sc_chan = soupx.autoEstCont(sc_chan, verbose=True)

    print("Ambient RNA estimation complete.")

    return adata, sc_chan

def remove_ambient_rna(
    adata: sc.AnnData,
    sc_chan,
    layer_name: str = "soupx_counts",
    overwrite: bool = False
) -> sc.AnnData:
    """
    -----------
    Parameters
    -----------
        adata: AnnData object.
        sc_chan: SoupChannel object.
        layer_name: Name of the layer in which corrected counts will be stored. 
        overwrite: Whether to overwrite an existing layer with the same name.
    --------
    Returns
    --------
        adata: AnnData object updated with adata.layers[layer_name].
    """

    if layer_name in adata.layers and not overwrite:
        raise ValueError(
            f"adata.layers['{layer_name}'] already exists. "
            f"Set overwrite=True if you want to replace it. "
        )
    corrected_matrix = soupx.adjustCounts(sc_chan)
    adata.layers[layer_name] = corrected_matrix.T.copy()

    print(f"Ambient RNA correction applied. Corrected counts stored in adata.layers['{layer_name}'].")

    original_sum = adata.layers["raw_counts"][0].sum() if "raw_counts" in adata.layers else adata.X[0].sum()
    corrected_sum = adata.layers[layer_name][0].sum()
    print(f"Orignal counts: {original_sum:.0f}")
    print(f"Corrected counts: {corrected_sum:.0f}")
    print(f"Counts removed: {original_sum - corrected_sum:.0f}")
    
    return adata