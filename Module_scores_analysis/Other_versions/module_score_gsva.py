#!/usr/bin/env python3
"""GSVA-based module score calculation for VST expression matrices.

This script uses the GSVA (Gene Set Variation Analysis) method to calculate
module scores for gene sets. GSVA is a non-parametric, unsupervised method
that calculates sample-wise gene set enrichment scores as a function of genes
inside and outside the gene set.

Reference: Hänzelmann et al., BMC Bioinformatics (2013)

Example:
    python module_score_gsva.py \
        --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
        --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
        --gene-column external_gene_name \
        --output scores_gsva.tsv \
        --heatmap scores_gsva_heatmap.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def check_dependencies() -> None:
    """
    Check if required dependencies are installed.
    
    Raises:
        ImportError: If gseapy is not installed.
    """
    try:
        import gseapy
        print(f"✓ gseapy version: {gseapy.__version__}")
    except ImportError:
        print("ERROR: gseapy is not installed.")
        print("\nPlease install it using:")
        print("  pip install gseapy")
        print("\nOr with conda:")
        print("  conda install -c bioconda gseapy")
        sys.exit(1)


def read_gene_list(path: Path) -> List[str]:
    """
    Read a plain-text gene list (one gene per line).

    Args:
        path: Path to the gene list file.

    Returns:
        List of gene identifiers (stripped, non-empty).
    
    Raises:
        FileNotFoundError: If the file doesn't exist.
        ValueError: If the file is empty.
    """
    if not path.exists():
        raise FileNotFoundError(f"Gene list file not found: {path}")
    
    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    
    if not genes:
        raise ValueError(f"Gene list file is empty: {path}")
    
    return genes


def compute_gsva_scores(
    expr: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> pd.DataFrame:
    """
    Compute GSVA scores for each gene set.

    Uses the GSVA algorithm (Hänzelmann et al., 2013) to calculate sample-wise
    enrichment scores for each gene set. GSVA is a non-parametric approach that
    evaluates variation of gene set enrichment through the samples.

    Args:
        expr: Expression matrix indexed by gene, columns are samples.
        gene_sets: Dict mapping gene set name to list of gene IDs.

    Returns:
        DataFrame: gene sets as rows, samples as columns.
    """
    import gseapy as gp
    
    # Convert index to set for O(1) lookup instead of O(n)
    expr_genes = set(expr.index)
    
    # Report gene overlap statistics before running GSVA
    print("\nGene set overlap with expression matrix:")
    for name, genes in gene_sets.items():
        # Use set intersection for faster lookup
        present = [g for g in genes if g in expr_genes]
        n_present = len(present)
        n_total = len(genes)
        pct_found = (n_present / n_total * 100) if n_total > 0 else 0
        print(f"  [{name}] {n_present}/{n_total} genes found ({pct_found:.1f}%)")
        
        if not present:
            raise ValueError(f"No genes from set '{name}' were found in the matrix.")
        
        if pct_found < 50:
            print(f"    ⚠️  WARNING: Less than 50% of genes found for '{name}'")
    
    print("\nRunning GSVA algorithm...")
    print("  Method: GSVA")
    print("  Kernel: Gaussian (for continuous VST data)")
    
    # Run GSVA
    # Note: gseapy.gsva expects genes as index, samples as columns (same as our input)
    # kcdf='Gaussian' is appropriate for continuous expression data like VST
    # mx_diff=True: ES is difference between max positive and max negative deviations (standard GSVA)
    # abs_ranking=False: ES is magnitude difference (default)
    try:
        gsva_results = gp.gsva(
            data=expr,
            gene_sets=gene_sets,
            method='gsva',
            kcdf='Gaussian',
            mx_diff=True,  # Standard GSVA behavior
            abs_ranking=False,
            max_size=5000,
            min_size=1,
            permutation_num=0,  # GSVA doesn't use permutations
            outdir=None,  # Don't save intermediate files
            no_plot=True,  # Don't generate plots
            verbose=True
        )
    except Exception as e:
        raise RuntimeError(f"GSVA computation failed: {e}")
    
    # Extract the enrichment scores
    # gsva_results is a DataFrame with gene sets as index, samples as columns
    scores = gsva_results
    
    return scores


def plot_heatmap(scores: pd.DataFrame, output: Path) -> None:
    """
    Save a heatmap of module scores to file.

    Args:
        scores: DataFrame of module scores (gene sets x samples).
        output: Path to save the heatmap image.
    """
    # Dynamically adjust figure size based on data dimensions
    n_sets, n_samples = scores.shape
    fig_width = max(8, n_samples * 0.8)
    fig_height = max(6, n_sets * 0.6)
    
    plt.figure(figsize=(fig_width, fig_height))
    
    # Center colormap at zero for visualizing up/down regulation and annotate exact values.
    sns.heatmap(scores, cmap="vlag", center=0, annot=True, fmt=".2f", 
                cbar_kws={'label': 'GSVA Enrichment Score'})
    
    plt.title('GSVA Module Scores', fontsize=14, fontweight='bold')
    plt.xlabel('Samples', fontsize=12)
    plt.ylabel('Gene Sets', fontsize=12)
    plt.tight_layout()
    
    try:
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        plt.close()
        raise IOError(f"Failed to save heatmap to {output}: {e}")


def parse_args() -> argparse.Namespace:
    """
    Set up and parse command-line arguments for the script.

    Returns:
        argparse.Namespace with parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compute GSVA-based module scores from VST matrix and gene lists."
    )
    parser.add_argument(
        "--vst",
        required=True,
        type=Path,
        help="Path to VST matrix (TSV with genes as rows)."
    )
    parser.add_argument(
        "--gene-lists",
        required=True,
        nargs="+",
        type=Path,
        help="Paths to gene list files."
    )
    parser.add_argument(
        "--gene-column",
        default="external_gene_name",
        help="Column name containing gene identifiers in the VST file (default: external_gene_name).",
    )
    parser.add_argument(
        "--output",
        default=Path("module_scores_gsva.tsv"),
        type=Path,
        help="Output TSV for module scores."
    )
    parser.add_argument(
        "--heatmap",
        type=Path,
        help="Optional path to save a heatmap (PNG/PDF/SVG supported by matplotlib).",
    )
    return parser.parse_args()


def main() -> None:
    """
    Main entry point: parse arguments, load data, compute GSVA scores, and save outputs.
    """
    # Check dependencies first
    check_dependencies()
    
    # Parse all CLI arguments
    args = parse_args()

    # Load VST expression matrix (genes x samples)
    print(f"\nLoading expression matrix from: {args.vst}")
    try:
        expr = pd.read_csv(args.vst, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(f"VST matrix file not found: {args.vst}")
    except Exception as e:
        raise IOError(f"Failed to read VST matrix: {e}")
    
    # Ensure the requested gene identifier column exists
    if args.gene_column not in expr.columns:
        available_cols = ", ".join(expr.columns[:5].tolist())
        raise ValueError(
            f"Column '{args.gene_column}' not found in VST file.\n"
            f"Available columns (first 5): {available_cols}..."
        )

    # Set gene column as index
    expr = expr.set_index(args.gene_column)
    
    # Report matrix dimensions
    print(f"Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")
    
    # Check for duplicate genes
    if expr.index.duplicated().any():
        n_duplicates = expr.index.duplicated().sum()
        print(f"⚠️  WARNING: {n_duplicates} duplicate gene identifiers found in matrix")
        print(f"   Keeping first occurrence of each duplicate.")
        expr = expr[~expr.index.duplicated(keep='first')]

    # Read gene lists, using file stem as set name
    gene_sets = {path.stem: read_gene_list(path) for path in args.gene_lists}
    print(f"\nLoaded {len(gene_sets)} gene sets:")
    for name, genes in gene_sets.items():
        print(f"  - {name}: {len(genes)} genes")

    # Compute GSVA scores
    scores = compute_gsva_scores(expr, gene_sets)
    
    # Display summary statistics
    print(f"\nGSVA score summary:")
    print(f"  Mean ± SD across all scores: {scores.values.mean():.3f} ± {scores.values.std():.3f}")
    print(f"  Range: [{scores.values.min():.3f}, {scores.values.max():.3f}]")
    
    # Save scores
    try:
        scores.to_csv(args.output, sep="\t")
        print(f"\n✓ Scores saved to: {args.output}")
    except Exception as e:
        raise IOError(f"Failed to save scores to {args.output}: {e}")

    # Optionally plot heatmap
    if args.heatmap:
        plot_heatmap(scores, args.heatmap)
        print(f"✓ Heatmap saved to: {args.heatmap}")


if __name__ == "__main__":
    main()
