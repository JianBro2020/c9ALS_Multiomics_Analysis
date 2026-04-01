import scanpy as sc
import scrublet as scr

def detect_doublets(
    adata: sc.AnnData,
    expected_doublet_rate: float | None = None,
    min_counts: int = 2, 
    min_cells: int = 3,
    min_gene_variability_pctl: int = 85,
    n_prin_comps: int = 30,
    use_approx_neighbors: bool = False,
    distance_metric: str = "euclidean"
) -> sc.AnnData:
    """
    -----------
    Parameters
    -----------
        adata: AnnData object. Must contain adata.layers["raw_counts"].
        expected_doublet_rate: Expected doublet rate for this sample. If None, estimate as 0.8% per 1000 recovered cells.
        min_counts: Minimum counts for a gene to be considered by Scrublet.
        min_cells: Minimum number of cells in which a gene must be detected.
        min_gene_variability_pctl: Percentile threshold for selecting variable genes when simulating doublets.
        n_prin_comps: Number of principal components used by Scrublet.
        use_approx_neighbors: Whether to use approximate nearest neighbors in Scrublet. 
        distance_metric: Distance metric used by Scrublet.
    --------
    Returns
    --------
        adata: AnnData object, updated with:
            - adata.obs["doublet_score_scrublet"]
            - adata.obs["predicted_doublet_scrublet"]
    """
    if "raw_counts" not in adata.layers:
        raise KeyError(
            "'raw counts' not found in adata.layers. "
            "Please run freeze_raw_counts() first."
        )
    
    if expected_doublet_rate is None:
        expected_doublet_rate = 0.008 * (adata.n_obs / 1000)
    
    scrub = scr.Scrublet(
        adata.layers["raw_counts"],
        expected_doublet_rate = expected_doublet_rate
    )

    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts = min_counts,
        min_cells = min_cells,
        min_gene_variability_pctl = min_gene_variability_pctl,
        n_prin_comps = n_prin_comps,
        use_approx_neighbors = use_approx_neighbors,
        distance_metric = distance_metric
    )

    adata.obs["doublet_score_scrublet"] = doublet_scores
    adata.obs["predicted_doublet_scrublet"] = predicted_doublets

    n_cells = adata.n_obs
    n_doublets = int(adata.obs["predicted_doublet_scrublet"].sum())
    predicted_rate = n_doublets / n_cells
    print(
        f"Doublet detection complete.\n"
        f"Cells analysed: {n_cells}\n"
        f"Expected doublet rate used: {expected_doublet_rate:.3%}\n"
        f"Predicted doublets: {n_doublets} ({predicted_rate:.3%})"
    )

    return adata, scrub

def remove_doublets(adata: sc.AnnData) -> sc.AnnData:
    """
    Keep only cells not predicted as doublets.
    -----------
    Parameters
    -----------
        adata: AnnData object alreay annotated by detect_doublets().
    --------
    Returns
    -------- 
        adata_singlets: New AnnData object containing only predicted singlets.
    """

    if "predicted_doublet_scrublet" not in adata.obs.columns:
        raise KeyError(
            "'predicted_doublet_scrublet' not found in adata.obs. " 
            "Please run detect_doublets() first."
        )

    n_before = adata.n_obs
    adata_singlets = adata[~adata.obs["predicted_doublet_scrublet"]].copy()
    n_after = adata_singlets.n_obs
    n_removed = n_before - n_after

    print(
        f"Predicted doublets removed. "
        f"\nCells before: {n_before}, after: {n_after}, removed: {n_removed}."
    )

    return adata_singlets