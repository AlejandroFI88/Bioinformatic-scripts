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
    """Read a plain-text gene list (one gene per line)."""

    # Strip whitespace and ignore empty lines so downstream calculations only
    # receive valid identifiers.
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def compute_module_scores(
    expr: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    center: bool = True,
) -> pd.DataFrame:
    """Compute module scores for each gene set.

    Args:
        expr: Expression matrix indexed by gene with samples in columns.
        gene_sets: Mapping of gene set name to list of genes.
        center: If True, subtract the global mean expression per sample.

    Returns:
        DataFrame with gene sets as rows and samples as columns.
    """
    # Global per-sample mean that will be subtracted from each gene set if
    # centering is requested. Using 0 when `center` is False keeps the code
    # simple without an additional branch later.
    background_mean = expr.mean(axis=0) if center else 0

    # Collect the mean for each gene set in a dictionary keyed by gene set name.
    scores = {}

    for name, genes in gene_sets.items():
        # Only keep genes that are present in the expression matrix. This keeps
        # the script robust to mismatched identifiers between VST and lists.
        present = [g for g in genes if g in expr.index]
        if not present:
            raise ValueError(f"No genes from set '{name}' were found in the matrix.")

        # Compute the per-sample mean for the genes found in the list.
        gene_mean = expr.loc[present].mean(axis=0)

        # Subtract the background mean (or 0) to center values if requested.
        scores[name] = gene_mean - background_mean

    return pd.DataFrame(scores).T


def plot_heatmap(scores: pd.DataFrame, output: Path) -> None:
    """Save a clustered heatmap of module scores."""
    # Create a new figure with a small, publication-friendly footprint.
    plt.figure(figsize=(8, 6))

    # `center=0` places white at zero to quickly highlight positive/negative
    # deviations across samples.
    sns.heatmap(scores, cmap="vlag", center=0, annot=True, fmt=".2f")

    # Avoid clipped labels and then write the image to disk.
    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()


def parse_args() -> argparse.Namespace:
    """Configure the command-line interface and parse arguments."""

    # All parameters are optional except the VST matrix and gene lists. Defaults
    # are chosen to mimic the example provided in the script docstring.
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
    args = parse_args()

    # Load the VST matrix. It is expected to be a tab-separated file with genes
    # as rows and samples in columns.
    expr = pd.read_csv(args.vst, sep="\t")
    if args.gene_column not in expr.columns:
        raise ValueError(f"Column '{args.gene_column}' not found in VST file: {args.vst}")

    # Move the gene identifier column into the index to simplify lookups.
    expr = expr.set_index(args.gene_column)

    # Read all provided gene lists into a dictionary using the stem of each file
    # as the gene set name (e.g., `2c_UP.txt` -> `2c_UP`).
    gene_sets = {path.stem: read_gene_list(path) for path in args.gene_lists}

    # Compute module scores and write them to disk as a TSV.
    scores = compute_module_scores(expr, gene_sets, center=not args.no_center)
    scores.to_csv(args.output, sep="\t")

    if args.heatmap:
        # Optionally render a small heatmap for a quick visual summary.
        plot_heatmap(scores, args.heatmap)


if __name__ == "__main__":
    main()
