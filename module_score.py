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
    # Read file once, drop newline characters, remove blanks, and return clean IDs for downstream lookup.
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
    # Calculate per-sample mean for centering (if requested) so each score reflects deviation from overall expression.
    background_mean = expr.mean(axis=0) if center else 0

    # Store module scores for each gene set in a dict that we later convert to a DataFrame.
    scores = {}

    for name, genes in gene_sets.items():
        # Filter the list to genes present in the expression matrix so we do not select missing rows.
        present = [g for g in genes if g in expr.index]
        
        # Report gene overlap statistics
        n_present = len(present)
        n_total = len(genes)
        pct_found = (n_present / n_total * 100) if n_total > 0 else 0
        print(f"[{name}] {n_present}/{n_total} genes found ({pct_found:.1f}%)")
        
        if not present:
            raise ValueError(f"No genes from set '{name}' were found in the matrix.")
        
        # Warn if less than 50% of genes are found
        if pct_found < 50:
            print(f"  ⚠️  WARNING: Less than 50% of genes found for '{name}'")

        # Mean expression per sample for this gene set computed across all available genes.
        gene_mean = expr.loc[present].mean(axis=0)

        # Center by subtracting background mean if requested.
        scores[name] = gene_mean - background_mean

    # Return scores as DataFrame (gene sets x samples) for easy downstream export/plotting.
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

    # Center colormap at zero for visualizing up/down regulation and annotate exact values.
    sns.heatmap(scores, cmap="vlag", center=0, annot=True, fmt=".2f")

    # Adjust layout to avoid clipped labels and save the figure to disk.
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
    # Parse all CLI arguments once at start to capture configuration.
    args = parse_args()

    # Load VST expression matrix (genes x samples).
    expr = pd.read_csv(args.vst, sep="\t")
    # Ensure the requested gene identifier column exists before reindexing.
    if args.gene_column not in expr.columns:
        raise ValueError(f"Column '{args.gene_column}' not found in VST file: {args.vst}")

    # Set gene column as index for easy lookup.
    expr = expr.set_index(args.gene_column)
    
    # Report matrix dimensions
    print(f"\nExpression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")
    
    # Check for duplicate genes
    if expr.index.duplicated().any():
        n_duplicates = expr.index.duplicated().sum()
        print(f"⚠️  WARNING: {n_duplicates} duplicate gene identifiers found in matrix")
        print(f"   Keeping first occurrence of each duplicate.")
        expr = expr[~expr.index.duplicated(keep='first')]

    # Read gene lists, using file stem as set name.
    gene_sets = {path.stem: read_gene_list(path) for path in args.gene_lists}
    print(f"\nLoaded {len(gene_sets)} gene sets:")
    for name, genes in gene_sets.items():
        print(f"  - {name}: {len(genes)} genes")

    # Compute module scores (mean VST per gene set, per sample).
    print(f"\nComputing module scores (centering: {not args.no_center})...")
    scores = compute_module_scores(expr, gene_sets, center=not args.no_center)
    
    # Display summary statistics
    print(f"\nModule score summary:")
    print(f"  Mean ± SD across all scores: {scores.values.mean():.3f} ± {scores.values.std():.3f}")
    print(f"  Range: [{scores.values.min():.3f}, {scores.values.max():.3f}]")
    
    # Persist scores as TSV so downstream analyses can reuse them.
    scores.to_csv(args.output, sep="\t")
    print(f"\n✓ Scores saved to: {args.output}")

    # Optionally plot heatmap if requested.
    if args.heatmap:
        plot_heatmap(scores, args.heatmap)
        print(f"✓ Heatmap saved to: {args.heatmap}")


if __name__ == "__main__":
    main()
