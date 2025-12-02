#!/usr/bin/env python3
"""Simple module score calculation for VST expression matrices.

The script takes a VST expression table (genes x samples) and a set of
plain-text gene lists (one gene per line). For each gene list it calculates a
module score per sample using the mean VST value of the genes present in the
matrix. An optional global centering step subtracts the mean VST value across
all genes for each sample.

Example:
    python scripts/module_score.py \
        --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
        --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
        --gene-column external_gene_name \
        --output scores.tsv \
        --heatmap scores_heatmap.png
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def read_gene_list(path: Path) -> List[str]:
    """
    Read a plain-text gene list (one gene per line).

    Args:
        path: Path to the gene list file.

    Returns:
        List of gene identifiers (stripped, non-empty).
    """
    # Read file, strip whitespace, ignore empty lines.
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]



def compute_module_scores(
    expr: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    center: bool = True,
) -> pd.DataFrame:
    """
    Compute module scores for each gene set.

    For each gene set, calculate the mean VST expression per sample, optionally
    subtracting the global mean expression for centering.

    Args:
        expr: Expression matrix indexed by gene, columns are samples.
        gene_sets: Dict mapping gene set name to list of gene IDs.
        center: If True, subtract global mean per sample.

    Returns:
        DataFrame: gene sets as rows, samples as columns.
    """
    # Calculate per-sample mean for centering (if requested).
    background_mean = expr.mean(axis=0) if center else 0

    scores = {}  # Store module scores for each gene set.

    for name, genes in gene_sets.items():
        # Filter to genes present in the expression matrix.
        present = [g for g in genes if g in expr.index]
        if not present:
            raise ValueError(f"No genes from set '{name}' were found in the matrix.")

        # Mean expression per sample for this gene set.
        gene_mean = expr.loc[present].mean(axis=0)

        # Center by subtracting background mean if requested.
        scores[name] = gene_mean - background_mean

    # Return scores as DataFrame (gene sets x samples).
    return pd.DataFrame(scores).T



def plot_heatmap(scores: pd.DataFrame, output: Path) -> None:
    """
    Save a heatmap of module scores to file.

    Args:
        scores: DataFrame of module scores (gene sets x samples).
        output: Path to save the heatmap image.
    """
    # Create figure for heatmap.
    plt.figure(figsize=(8, 6))

    # Center colormap at zero for visualizing up/down regulation.
    sns.heatmap(scores, cmap="vlag", center=0, annot=True, fmt=".2f")

    # Adjust layout and save image.
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()



def parse_args() -> argparse.Namespace:
    """
    Set up and parse command-line arguments for the script.

    Returns:
        argparse.Namespace with parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Compute module scores from VST matrix and gene lists.")
    parser.add_argument("--vst", required=True, type=Path, help="Path to VST matrix (TSV with genes as rows).")
    parser.add_argument("--gene-lists", required=True, nargs="+", type=Path, help="Paths to gene list files.")
    parser.add_argument(
        "--gene-column",
        default="external_gene_name",
        help="Column name containing gene identifiers in the VST file (default: external_gene_name).",
    )
    parser.add_argument("--output", default=Path("module_scores.tsv"), type=Path, help="Output TSV for module scores.")
    parser.add_argument(
        "--no-center",
        action="store_true",
        help="Disable centering by subtracting the global mean expression per sample.",
    )
    parser.add_argument(
        "--heatmap",
        type=Path,
        help="Optional path to save a heatmap (PNG/PDF/SVG supported by matplotlib).",
    )
    return parser.parse_args()



def main() -> None:
    """
    Main entry point: parse arguments, load data, compute scores, and save outputs.
    """
    args = parse_args()

    # Load VST expression matrix (genes x samples).
    expr = pd.read_csv(args.vst, sep="\t")
    if args.gene_column not in expr.columns:
        raise ValueError(f"Column '{args.gene_column}' not found in VST file: {args.vst}")

    # Set gene column as index for easy lookup.
    expr = expr.set_index(args.gene_column)

    # Read gene lists, using file stem as set name.
    gene_sets = {path.stem: read_gene_list(path) for path in args.gene_lists}

    # Compute module scores (mean VST per gene set, per sample).
    scores = compute_module_scores(expr, gene_sets, center=not args.no_center)
    scores.to_csv(args.output, sep="\t")

    # Optionally plot heatmap if requested.
    if args.heatmap:
        plot_heatmap(scores, args.heatmap)


if __name__ == "__main__":
    main()
