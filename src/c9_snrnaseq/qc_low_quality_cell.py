from pathlib import Path
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

def compute_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """
    Add QC gene flags and compute per-cell QC metrics.
    -----------
    Parameters
    -----------
        adata: AnnData object
    --------
    Returns
    --------
        adata, updated with:
            - adata.var["mt"]
            - adata.var["ribo"]
            - adata.var["hb"]
            - QC metrics in adata.obs via scanpy.pp.calculate_qc_metrics()
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDEG]")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[20],
        log1p=True,
        inplace=True,
    )

    print("QC metrics computed")

    return adata

def detect_low_quality_cells(
    adata: sc.AnnData,
    nmads_counts: int = 5,
    nmads_genes: int = 5,
    nmads_top20: int = 5,
    nmads_mt: int = 3,
    mt_percent_threshold: float = 5.0
) -> sc.AnnData:
    """
    Detect low-quality cells using MAD-based outlier logic and mitochondrial thresholds.
    -----------
    Parameters
    -----------
        adata: with QC annotations from compute_qc_metrics().
        nmads_counts: MAD threshold for log1p_total_counts.
        nmads_genes: MAD threshold for log1p_n_genes_by_counts.
        nmads_top20: MAD threshold for pct_counts_in_top_20_genes.
        nmads_mt: MAD threshold for pct_counts_mt.
        mt_percent_threshold: Hard threshold for pct_counts
    --------
    Returns
    --------
        adata, updated with 
            - adata.obs["outliers"]
            - adata.obs["mt_outliers"]
            - adata.obs["high_quality_cells"]
    """

    def is_outlier(metric: str, nmads: int) -> pd.Series:
        values = adata.obs[metric]
        median = np.median(values)
        mad = median_abs_deviation(values)

        lower_bound = median - (nmads * mad)
        upper_bound = median + (nmads * mad)

        return (values < lower_bound) | (values > upper_bound)

    adata.obs["outliers"] = (
        is_outlier("log1p_total_counts", nmads_counts) |
        is_outlier("log1p_n_genes_by_counts", nmads_genes) |
        is_outlier("pct_counts_in_top_20_genes", nmads_top20)
    )

    mt_outlier_mad = is_outlier("pct_counts_mt", nmads_mt)
    mt_outlier_percent = adata.obs["pct_counts_mt"] > mt_percent_threshold

    adata.obs["mt_outliers"] = mt_outlier_mad | mt_outlier_percent

    adata.obs["high_quality_cells"] = (~adata.obs["outliers"]) & (~adata.obs["mt_outliers"])

    n_outliers = (~adata.obs.high_quality_cells).sum()
    print(f"Detection complete. Flagged {n_outliers} cells as low-quality.")

    return adata

def remove_low_quality_cells(adata: sc.AnnData) -> sc.AnnData:
    """
    Keep only cells where "high_quality_cells" == True.
    -----------
    Parameters
    -----------
        adata: AnnData object that has been annotated by detect_low_quality_cells().
    --------
    Returns
    --------
        adata_filtered: AnnData object containing only high-quality cells.
    """

    if "high_quality_cells" not in adata.obs.columns:
        raise KeyError(
            "'high_quality_cells' not found in adata.obs."
            "Make sure of running detect_low_quality_cells() first."
        )
    n_before = adata.n_obs
    adata_filtered = adata[adata.obs["high_quality_cells"]].copy()
    n_after = adata_filtered.n_obs
    n_removed = n_before - n_after
    print(
        f"Low-quality cells removed."
        f"Cells before: {n_before}, after: {n_after}, removed: {n_removed}."
    )
    return adata_filtered

def freeze_raw_counts(adata: sc.AnnData) -> sc.AnnData:
    """
    Freeze the current count matrix as raw integer counts in a layer.
    This is typically called after low-quality-cell removal and before ambient RNA correction or any normalization and transformation.
    -----------
    Parameters
    -----------
        adata: AnnData object whose current 'X' is still the raw count data.
    --------
    Returns
    --------
        adata: AnnData object, updated with adata.layers["raw_counts"]
    """
    adata.layers["raw_counts"] = adata.X.copy()
    print("Raw counts frozen in adata.layers['raw_counts'].")
    
    return adata