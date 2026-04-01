from pathlib import Path
import scanpy as sc
import pandas as pd

from c9_snrnaseq.io_utils import (
    PROJECT_ROOT,
    load_and_annotate_sheet,
    save_checkpoint,
    _stage_header,
    _log,
)
from c9_snrnaseq.qc_low_quality_cell import (
    compute_qc_metrics,
    detect_low_quality_cells,
    remove_low_quality_cells,
    freeze_raw_counts,
)
from c9_snrnaseq.ambient_rna import (
    estimate_ambient_rna,
    remove_ambient_rna,
)
from c9_snrnaseq.doublets_removal import (
    detect_doublets,
    remove_doublets,
)

from c9_snrnaseq.expression_preprocessing import (
    normalize_counts,
    log_transform,
    select_hvgs,
    scale_expression
)

from c9_snrnaseq.dimensionality_reduction import (
    run_pca,
    build_neighbor_graph,
    run_umap
)

from c9_snrnaseq.annotation import (
    run_leiden_clustering,
    find_cluster_markers
)


def process_one_sample(
    sample_row: pd.Series,
    processed_dir: str | Path = "data/processed/per_sample",
    save_intermediate: bool = True,
    verbose: bool = True
) -> dict:
    """
    Process one sample through sample-level low-quality-cell filtering, ambient RNA correction, and doublet detection / removal.
    Workflow:
        1. Load filtered count matrix and annotate metadata.
        2. Compute QC metrics.
        3. Detect low-quality cells.
        4. Save checkpoint before removing the low-quality cells.
        5. Remove low-quality cells.
        6. Freeze raw counts.
        7. Save checkpoint after removing the low-quality cells.
        8. Load matching raw droplet matrix.
        9. Estimate and remove ambient RNA.
        10. Save checkpoint after ambient RNA correction.
        11. Detect doublets.
        12. Save checkpoint before doublet removal.
        13. Remove doublets.
        14. Save checkpoint after doublet removal.
    -----------
    Parameters
    -----------
        sample_row: One row from the sample sheet (pandas Series). Must contain:
            - sample_id
            - condition 
            - filtered_matrix_dir
            - raw_matrix_dir
            Optional fileds such as dataset will be added to 'adata.obs' by load_and_annotate_sheet() if present.
        processed_dir: Base directory where per-sample checkpoint '.h5ad' files will be saved. 
        save_intermediate: Whether to save checkpoint files after key stages.
        verbose: Whether to print concise stage-level progress messages.
    --------
    Returns
    --------
        result: a dictionary object containing
            - "sample_id": sample_identifier 
            - "adata_postQC": final adata object after the three stages of cleaning.
            - "sc_chan": SoupX channel object
            - "scrub": Scrublet object
            - "checkpoint_paths": dict of saved checkpoint paths
    """
    total_stages = 5

    sample_id = str(sample_row["sample_id"])
    condition = str(sample_row['condition'])
    processed_dir = Path(processed_dir)
    if not processed_dir.is_absolute():
        processed_dir = PROJECT_ROOT / processed_dir
    sample_outdir = processed_dir / sample_id
    sample_outdir.mkdir(parents=True, exist_ok=True)

    checkpoint_paths: dict[str, Path] = {}

    _log(f"\nProcessing sample: {sample_id} ({condition})", verbose)

    # Load filtered count matrix.
    _stage_header(1, total_stages, "Load sample", verbose)
    adata = load_and_annotate_sheet(sample_row)
    _log(f"Loaded sample {sample_id}: {adata.n_obs} cells x {adata.n_vars} genes", verbose)
    # Compute QC metrics
    _stage_header(2, total_stages, "Low-quality-cell filtering", verbose)
    adata = compute_qc_metrics(adata)
    # Detect + Remove low-quality cells
    adata = detect_low_quality_cells(adata)
    if save_intermediate:
        checkpoint_paths["before_low_quality_cell_filtering"] = save_checkpoint(
            adata = adata,
            output_dir = sample_outdir,
            stage_name = "before_low_quality_cell_filtering" 
        )
    adata = remove_low_quality_cells(adata)
    _log(f"Low-quality-cell filtering complete: ", verbose)

    # Freeze raw counts
    _stage_header(3, total_stages, "Freeze raw counts and save checkpoint", verbose)
    adata = freeze_raw_counts(adata)
    if save_intermediate:
        checkpoint_paths["after_low_quality_cell_filtering"] = save_checkpoint(
            adata = adata,
            output_dir = sample_outdir,
            stage_name = "after_low_quality_cell_filtering",
        )
    _log("Raw counts frozen and post-filter checkpoint saved.", verbose)

    # Load raw droplet matrix
    _stage_header(4, total_stages, "Ambient RNA correction", verbose)
    raw_dir = Path(sample_row["raw_matrix_dir"])
    if not raw_dir.is_absolute():
        raw_dir = PROJECT_ROOT / raw_dir
    required_files = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
    missing_files = [name for name in required_files if not (raw_dir / name).exists()]
    if missing_files:
        raise FileNotFoundError(
            f"Raw matrix directory is incomplete: {raw_dir}. Missing files: {missing_files}"
        )
    
    adata_raw = sc.read_10x_mtx(
        raw_dir,
        var_names = "gene_symbols",
        cache = False
    )
    adata_raw.var_names_make_unique()

    # Estimate + Remove ambient RNA
    adata, sc_chan = estimate_ambient_rna(
        adata = adata,
        adata_raw = adata_raw
    )
    adata = remove_ambient_rna(
        adata = adata,
        sc_chan = sc_chan
    )
    if save_intermediate:
        checkpoint_paths["after_ambient_rna_removal"] = save_checkpoint(
            adata = adata,
            output_dir = sample_outdir,
            stage_name = "after_ambient_rna_removal"
        )
    _log("Ambient RNA correction complete.", verbose)

    # Detect + Remove doublets
    _stage_header(5, total_stages, "Doublet detection and removal", verbose)
    adata, scrub = detect_doublets(adata)
    if save_intermediate:
        checkpoint_paths["before_doublet_removal"] = save_checkpoint(
            adata = adata,
            output_dir = sample_outdir,
            stage_name = "before_doublet_removal",
        )
    adata = remove_doublets(adata)
    if save_intermediate:
        checkpoint_paths["after_doublet_removal"] = save_checkpoint(
            adata = adata,
            output_dir = sample_outdir,
            stage_name = "after_doublet_removal"
        )
    _log(f"Doublet removal complete.", verbose)
    
    _log("\nAvailable checkpoint objects for before/after inspection:", verbose)
    for key, path in checkpoint_paths.items():
        _log(f"  {key} -> {path}", verbose)

    _log(f"\nFinished processing sample: {sample_id}\n", verbose)
    
    return {
        "sample_id": sample_id,
        "adata_postQC": adata,
        "sc_chan": sc_chan,
        "scrub": scrub, 
        "checkpoint_paths": checkpoint_paths
    }


def process_all_samples(
        sample_sheet: str | Path = "config/samples.csv",
        processed_dir: str | Path = "data/processed/per_sample",
        merged_output_path: str | Path = "data/processed/merged/merged_after_clustering.h5ad",
        save_intermediate: bool = True,
        verbose: bool = True
) -> dict:
    """
    Process all samples listed in the sample sheet and save merged output.
    -----------
    Parameters
    -----------
        sample_sheet: Path to the sample sheet CSV file. Must contain columns:
            - sample_id
            - condition 
            - filtered_matrix_dir
            - raw_matrix_dir
        processed_dir: Base directory where per-sample checkpoint '.h5ad' files will be saved. 
        merged_output_path: Path where the merged adata object after processing all samples will be saved.
        save_intermediate: Whether to save checkpoint files after key stages for each sample.
        verbose: Whether to print concise stage-level progress messages.
    --------
    Returns
    --------
        result: a dictionary object containing
            - "per_sample_results": list of per-sample result dictionaries from process_one_sample()
            - "adata_merged": merged adata object after preprocessing, clustering, and marker detection
    """
    
    sample_sheet = Path(sample_sheet)
    if not sample_sheet.is_absolute():
        sample_sheet = PROJECT_ROOT / sample_sheet
    sample_df = pd.read_csv(sample_sheet)
    per_sample_results = []
    cleaned_adatas = []
    cleaned_paths = []

    for _, sample_row in sample_df.iterrows():
        result = process_one_sample(
            sample_row = sample_row,
            processed_dir = processed_dir,
            save_intermediate = save_intermediate,
            verbose = verbose
        )
        per_sample_results.append(result)
        cleaned_adatas.append(result["adata_postQC"])
        
        checkpoint_path = result["checkpoint_paths"].get("after_doublet_removal")
        if checkpoint_path is not None:
            cleaned_paths.append(checkpoint_path)
    
    adata_merged = sc.concat(
        cleaned_adatas,
        join="outer",
        label="concat_sample_id",
        keys=[adata.obs["sample"].iloc[0] for adata in cleaned_adatas],
        merge="same",
    )
    adata_merged.obs_names_make_unique()

    adata_merged = normalize_counts(adata_merged)
    adata_merged = log_transform(adata_merged)
    adata_merged = select_hvgs(adata_merged)
    adata_merged = scale_expression(adata_merged)
    adata_merged = run_pca(adata_merged)
    adata_merged = build_neighbor_graph(adata_merged)
    adata_merged = run_umap(adata_merged)
    adata_merged = run_leiden_clustering(adata_merged)
    adata_merged = find_cluster_markers(adata_merged)

    merged_output_path = Path(merged_output_path)
    if not merged_output_path.is_absolute():
        merged_output_path = PROJECT_ROOT / merged_output_path
    merged_output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_merged.write_h5ad(merged_output_path)
    _log(f"Merged dataset saved to: {merged_output_path}", verbose)

    return {
        "sample_table": sample_df,
        "per_sample_results": per_sample_results,
        "cleaned_paths": cleaned_paths,
        "adata_merged": adata_merged,
        "merged_output_path": merged_output_path,
    }
