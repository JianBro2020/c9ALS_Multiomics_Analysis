from pathlib import Path
import scanpy as sc
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]

def load_and_annotate(
    sample_dir: str|Path,
    sample_id: str,
    condition: str,
) -> sc.AnnData:
    """
    Load a sample's 10x-style count matrix and add basic metadata.
    -----------
    Parameters
    -----------
        sample_dir: Path to the directory that contains matrix.mtx.gz.
        sample_id: Unique sample identifier, e.g., 'GSM5292146' in this study for a c9ALS dataset.
        condition: Biological condition, e.g., 'c9ALS' or 'Control'.
    --------
    Returns
    --------
        adata: AnnData object for that sample.
    """
    sample_dir = Path(sample_dir)
    
    if not sample_dir.exists():
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

    adata = sc.read_10x_mtx(
        sample_dir,
        var_names="gene_symbols",
        cache=False 
    )
    adata.var_names_make_unique()

    adata.obs["sample"] = sample_id
    adata.obs["condition"] = condition

    print("Sample successfully loaded.")

    return adata

def load_and_annotate_sheet(sample_row: pd.Series) -> sc.AnnData:
    """
    Load a sample's 10x-style count matrix and add basic metadata.
    -----------
    Parameters
    -----------
        sample_row: One row from the sample sheet (pandas Series). Must contain:
            - sample_id
            - condition
            - filtered_matrix_dir
        add optional metadata if it is present:
            - dataset
    --------
    Returns
    --------
        adata: AnnData object for that sample.
    """
    sample_dir = Path(sample_row["filtered_matrix_dir"])
    if not sample_dir.is_absolute():
        sample_dir = PROJECT_ROOT / sample_dir
    if not sample_dir.exists():
        raise FileNotFoundError(f"Sample directory not found: {sample_dir}")

    adata = sc.read_10x_mtx(
        sample_dir,
        var_names="gene_symbols",
        cache=False 
    )
    adata.var_names_make_unique()

    adata.obs["sample"] = str(sample_row["sample_id"])
    adata.obs["condition"] = str(sample_row["condition"])

    optional_meta = ["dataset", "tissue", "source_name"]
    for meta in optional_meta:
        if meta in sample_row.index and pd.notna(sample_row[meta]):
            adata.obs[meta] = str(sample_row[meta])

    print(f"Sample successfully loaded: {sample_row['sample_id']}")

    return adata


def save_checkpoint(
    adata: sc.AnnData,
    output_dir: str | Path,
    stage_name: str
) -> Path:
    """
    Save an AnnData checkpoint as an .h5ad file.
    -----------
    Parameters
    -----------
        adata: AnnData object to save.
        output_dir: Directory where the checkpoint will be written to as a .h5ad file.
        stage_name: A description of the current stage, e.g., 'low_quality_cells_filtered', 'ambient_rna_removed'.
    --------
    Returns
    --------
        output_path: Path to the written .h5ad file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    condition = adata.obs["condition"].iloc[0] if "condition" in adata.obs.columns else "unknown_condition"
    sample_id = adata.obs["sample"].iloc[0] if "sample" in adata.obs.columns else "unknown_sample"
    output_path = output_dir / f"{condition}_{sample_id}_{stage_name}.h5ad"

    adata.write_h5ad(output_path)
    print(f"Checkpoint saved: {output_path}")

    return output_path

def _stage_header(
        stage_no: int,
        total_stages: int,
        title: str,
        verbose: bool = True
) -> None:
    if verbose:
        print(f"\n========== [{stage_no}/{total_stages}] {title} ==========\n")

def _log(
        msg: str,
        verbose: bool = True
) -> None:
    if verbose:
        print(msg)
