"""
Essential per-sample preprocessing for one STARsolo sample.

This script is intended for SLURM-array production runs. It performs the
sample-level work previously developed in notebooks 01-07, but saves only the
final object needed for later merging:

    data/processed/per_sample/{sample_id}/{condition}_{sample_id}_after_doublet_removal.h5ad
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc

SRC_DIR = Path(__file__).resolve().parents[1]
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from c9_snrnaseq.ambient_rna import estimate_ambient_rna, remove_ambient_rna
from c9_snrnaseq.doublets_removal import detect_doublets, remove_doublets
from c9_snrnaseq.io_utils import PROJECT_ROOT, load_and_annotate_sheet
from c9_snrnaseq.qc_low_quality_cell import (
    compute_qc_metrics,
    detect_low_quality_cells,
    freeze_raw_counts,
    remove_low_quality_cells,
)


def project_path(path: str | Path) -> Path:
    path = Path(path)
    return path if path.is_absolute() else PROJECT_ROOT / path


def select_sample(sample_sheet: Path, sample_id: str) -> pd.Series:
    samples = pd.read_csv(sample_sheet)
    samples["sample_id"] = samples["sample_id"].astype(str)

    match = samples.loc[samples["sample_id"].eq(str(sample_id))]
    if match.empty:
        raise ValueError(f"Sample not found in {sample_sheet}: {sample_id}")
    if len(match) > 1:
        raise ValueError(f"Sample appears more than once in {sample_sheet}: {sample_id}")

    return match.iloc[0]


def check_matrix_dir(matrix_dir: Path) -> None:
    missing = [
        name
        for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")
        if not (matrix_dir / name).exists()
    ]
    if missing:
        raise FileNotFoundError(f"Incomplete matrix directory: {matrix_dir}. Missing: {missing}")


def process_sample(sample_row: pd.Series, output_base: Path) -> Path:
    sample_id = str(sample_row["sample_id"])
    condition = str(sample_row["condition"])

    filtered_dir = project_path(sample_row["filtered_matrix_dir"])
    raw_dir = project_path(sample_row["raw_matrix_dir"])
    check_matrix_dir(filtered_dir)
    check_matrix_dir(raw_dir)

    print(f"\nProcessing {sample_id} ({condition})")
    print(f"Filtered matrix: {filtered_dir}")
    print(f"Raw matrix:      {raw_dir}")

    adata = load_and_annotate_sheet(sample_row)
    print(f"Loaded filtered cells: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    adata = compute_qc_metrics(adata)
    adata = detect_low_quality_cells(adata)
    adata = remove_low_quality_cells(adata)
    adata = freeze_raw_counts(adata)

    adata_raw = sc.read_10x_mtx(raw_dir, var_names="gene_symbols", cache=False)
    adata_raw.var_names_make_unique()
    print(f"Loaded raw droplets: {adata_raw.n_obs:,} droplets x {adata_raw.n_vars:,} genes")

    adata, _ = estimate_ambient_rna(adata=adata, adata_raw=adata_raw)
    adata = remove_ambient_rna(adata=adata, sc_chan=_)

    adata, _ = detect_doublets(adata)
    adata = remove_doublets(adata)

    sample_outdir = output_base / sample_id
    sample_outdir.mkdir(parents=True, exist_ok=True)
    output_path = sample_outdir / f"{condition}_{sample_id}_after_doublet_removal.h5ad"
    adata.write_h5ad(output_path)

    print(f"Saved final per-sample object: {output_path}")
    print(f"Final cells: {adata.n_obs:,}")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run essential per-sample preprocessing for one STARsolo sample."
    )
    parser.add_argument("--sample-id", required=True, help="Sample ID, for example GSM5292178.")
    parser.add_argument(
        "--sample-sheet",
        default="data/metadata/samples.csv",
        help="CSV with sample_id, condition, filtered_matrix_dir, and raw_matrix_dir.",
    )
    parser.add_argument(
        "--output-dir",
        default="data/processed/per_sample",
        help="Base directory for final per-sample h5ad outputs.",
    )
    args = parser.parse_args()

    sample_sheet = project_path(args.sample_sheet)
    output_dir = project_path(args.output_dir)

    sample_row = select_sample(sample_sheet, args.sample_id)
    process_sample(sample_row=sample_row, output_base=output_dir)


if __name__ == "__main__":
    main()
