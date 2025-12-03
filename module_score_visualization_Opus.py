#!/usr/bin/env python3
"""
Module Score Visualization for Bulk RNA-seq Data
=================================================

This script calculates module scores from VST expression matrices and generates
publication-ready visualizations similar to Seurat's AddModuleScore plots
(e.g., Figure 3G/3H style from scRNA-seq papers).

Algorithm - Seurat AddModuleScore Adaptation for Bulk RNA-seq:
    1. Bin all genes by average expression (default: 24 bins)
    2. For each gene in the target set:
       - Find its expression bin
       - Sample control genes from the same bin (excluding target genes)
    3. Calculate: Score = mean(target genes) - mean(control genes)
    
    Key Adaptation for Bulk vs scRNA-seq:
    - Bulk RNA-seq (default): 1 control gene per target gene (1:1 ratio)
    - scRNA-seq (Seurat): 100 control genes per target gene (100:1 ratio)
    
    The 1:1 ratio is more appropriate for Bulk RNA-seq because:
    - Fewer samples (biological replicates) vs thousands of cells
    - Less need for extensive control averaging
    - Reduces risk of oversampling from limited gene pool

Output visualizations (similar to Figure 3G/3H from papers):
- Boxplots/Violin plots comparing module scores across conditions (with statistics)
- Heatmaps with hierarchical clustering
- Bar plots with error bars (mean ± SEM)
- Dot plots showing score magnitude and direction

Usage:
    # Basic usage with sample groups
    python module_score_visualization_Opus.py \
        --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
        --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
        --gene-column external_gene_name \
        --sample-groups groups.txt \
        --output-prefix results/module_scores \
        --plots boxplot heatmap

    # Auto-infer groups from sample names
    python module_score_visualization_Opus.py \
        --vst expression_vst.txt \
        --gene-lists pathway1.txt pathway2.txt \
        --infer-groups \
        --output-prefix results/figure3_style

References:
    - Seurat AddModuleScore: Tirosh et al., Science (2016)
    - Paper example: Lu et al., Cell (2025) Figure 3G/3H
    - GSVA alternative: Hänzelmann et al., BMC Bioinformatics (2013)
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)


# =============================================================================
# Gene List and Data Loading Functions
# =============================================================================

def read_gene_list(path: Path) -> List[str]:
    """Read a plain-text gene list (one gene per line)."""
    if not path.exists():
        raise FileNotFoundError(f"Gene list file not found: {path}")
    genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not genes:
        raise ValueError(f"Gene list file is empty: {path}")
    return genes


def read_sample_groups(path: Path) -> pd.DataFrame:
    """
    Read sample group assignments from a TSV file.
    
    Expected format (tab-separated):
        sample_id    group
        sample1      Control
        sample2      Treatment
        ...
    
    Args:
        path: Path to the groups file.
        
    Returns:
        DataFrame with 'sample_id' and 'group' columns.
    """
    if not path.exists():
        raise FileNotFoundError(f"Sample groups file not found: {path}")
    
    groups = pd.read_csv(path, sep='\t')
    
    # Check for required columns
    if 'sample_id' not in groups.columns or 'group' not in groups.columns:
        # Try to infer from first two columns
        if len(groups.columns) >= 2:
            groups.columns = ['sample_id', 'group'] + list(groups.columns[2:])
            print(f"  ⚠️  Assumed first column is 'sample_id' and second is 'group'")
        else:
            raise ValueError(
                "Sample groups file must have 'sample_id' and 'group' columns"
            )
    
    return groups[['sample_id', 'group']]


def infer_groups_from_sample_names(sample_names: List[str]) -> pd.DataFrame:
    """
    Attempt to infer groups from sample names using common patterns.
    
    Common patterns detected:
        - Prefix before underscore: Control_1, Treatment_2
        - Prefix before number: WT1, KO2
        - Replicate patterns: sample_rep1, sample_rep2
    
    Args:
        sample_names: List of sample names.
        
    Returns:
        DataFrame with 'sample_id' and 'group' columns.
    """
    groups = []
    
    for name in sample_names:
        # Try pattern: group_replicate (e.g., Control_1, shNsun2_rep1)
        match = re.match(r'^(.+?)_(?:rep)?(\d+)$', name, re.IGNORECASE)
        if match:
            groups.append({'sample_id': name, 'group': match.group(1)})
            continue
        
        # Try pattern: group + number at end (e.g., Control1, WT2)
        match = re.match(r'^(.+?)(\d+)$', name)
        if match:
            groups.append({'sample_id': name, 'group': match.group(1)})
            continue
        
        # Default: use the whole name as the group
        groups.append({'sample_id': name, 'group': name})
    
    df = pd.DataFrame(groups)
    
    # Report inferred groups
    unique_groups = df['group'].unique()
    print(f"  Inferred {len(unique_groups)} groups: {', '.join(unique_groups)}")
    
    return df


# =============================================================================
# Module Score Calculation (Seurat-like)
# =============================================================================

def _prepare_bins(expr: pd.DataFrame, bins: int) -> Tuple[pd.Series, Dict[int, List[str]]]:
    """Bin genes by average expression for Seurat-like control sampling."""
    gene_means = expr.mean(axis=1)
    ranked = gene_means.rank(method="first")
    n_bins = min(bins, ranked.nunique())
    gene_bins = pd.qcut(ranked, q=n_bins, labels=False, duplicates="drop")

    bin_to_genes: Dict[int, List[str]] = {}
    for gene, bin_id in gene_bins.items():
        bin_to_genes.setdefault(int(bin_id), []).append(gene)

    return gene_bins, bin_to_genes


def compute_module_scores(
    expr: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    method: str = "seurat-like",
    bins: int = 24,
    seed: int = 42,
    ctrl_genes_per_gene: int = 1,
) -> pd.DataFrame:
    """
    Compute module scores for each gene set.

    Methods:
    - ``seurat-like``: Seurat AddModuleScore-style with control genes sampled 
      from expression-matched bins. Score = mean(gene set) - mean(control genes).
      For Bulk RNA-seq, uses 1:1 ratio (1 control per gene) vs scRNA-seq's 100:1.
    - ``mean``: Simple mean expression for the gene set minus global mean.

    Args:
        expr: Expression matrix indexed by gene, columns are samples.
        gene_sets: Dict mapping gene set name to list of gene IDs.
        method: Either ``"seurat-like"`` or ``"mean"``.
        bins: Number of expression bins for control sampling (default: 24, matching Seurat).
        seed: Random seed for reproducible control gene sampling.
        ctrl_genes_per_gene: Number of control genes per gene in set. 
            For Bulk RNA-seq: 1 (default), For scRNA-seq: 100 (Seurat default).

    Returns:
        DataFrame: gene sets as rows, samples as columns.
    """
    rng = np.random.default_rng(seed)
    expr_genes = set(expr.index)

    # Prepare bins if using Seurat-like method
    gene_bins = None
    bin_to_genes: Dict[int, List[str]] = {}
    if method == "seurat-like":
        gene_bins, bin_to_genes = _prepare_bins(expr, bins)
        print(f"  Using {len(bin_to_genes)} expression bins for control genes")

    background_mean = expr.mean(axis=0) if method == "mean" else 0
    scores = {}

    for name, genes in gene_sets.items():
        present = [g for g in genes if g in expr_genes]
        n_present, n_total = len(present), len(genes)
        pct_found = (n_present / n_total * 100) if n_total > 0 else 0
        print(f"  [{name}] {n_present}/{n_total} genes found ({pct_found:.1f}%)")

        if not present:
            raise ValueError(f"No genes from set '{name}' found in matrix.")
        if pct_found < 50:
            print(f"    ⚠️  WARNING: Less than 50% of genes found")

        gene_mean = expr.loc[present].mean(axis=0)

        if method == "mean":
            scores[name] = gene_mean - background_mean
            continue

        # Seurat-like: sample control genes from expression-matched bins
        assert gene_bins is not None
        control_genes: List[str] = []
        for gene in present:
            bin_id = int(gene_bins[gene])
            # Select control genes from same expression bin, excluding target genes
            candidates = [g for g in bin_to_genes[bin_id] if g not in present]
            # If all genes in bin belong to the target set, fall back to the full bin
            if not candidates:
                candidates = bin_to_genes[bin_id]
            
            # Sample control genes: for Bulk RNA-seq, typically 1:1 ratio
            # For scRNA-seq compatibility, can increase ctrl_genes_per_gene to 100
            n_sample = min(ctrl_genes_per_gene, len(candidates))
            if n_sample == 1:
                # Optimized path for Bulk RNA-seq (default)
                control_genes.append(rng.choice(candidates))
            else:
                # scRNA-seq style: sample multiple controls per gene
                control_genes.extend(rng.choice(candidates, size=n_sample, replace=False))

        control_mean = expr.loc[control_genes].mean(axis=0)
        scores[name] = gene_mean - control_mean

    return pd.DataFrame(scores).T


# =============================================================================
# Visualization Functions
# =============================================================================

def plot_module_boxplot(
    scores: pd.DataFrame,
    groups: pd.DataFrame,
    output: Path,
    palette: Optional[Dict[str, str]] = None,
    show_points: bool = True,
    figsize: Tuple[int, int] = (10, 6),
) -> None:
    """
    Create boxplots of module scores grouped by condition.
    
    Similar to Figure 3G/3H style from scRNA-seq papers.
    
    Args:
        scores: Module scores (gene sets x samples).
        groups: Sample-to-group mapping.
        output: Output file path.
        palette: Color palette for groups.
        show_points: Whether to overlay individual data points.
        figsize: Figure dimensions.
    """
    # Convert scores to long format for plotting
    scores_long = scores.T.reset_index()
    scores_long.columns = ['sample_id'] + list(scores.index)
    scores_long = scores_long.melt(
        id_vars=['sample_id'],
        var_name='gene_set',
        value_name='module_score'
    )
    
    # Merge with group information
    scores_long = scores_long.merge(groups, on='sample_id', how='left')
    
    # Handle samples without group assignment
    if scores_long['group'].isna().any():
        missing = scores_long[scores_long['group'].isna()]['sample_id'].unique()
        print(f"  ⚠️  {len(missing)} samples without group assignment")
        scores_long['group'] = scores_long['group'].fillna('Unknown')
    
    # Create figure with subplots for each gene set
    gene_sets = scores.index.tolist()
    n_sets = len(gene_sets)
    n_cols = min(2, n_sets)
    n_rows = (n_sets + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0], figsize[1] * n_rows / 2))
    if n_sets == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Default palette if not provided
    if palette is None:
        unique_groups = scores_long['group'].unique()
        colors = sns.color_palette("Set2", len(unique_groups))
        palette = dict(zip(unique_groups, colors))
    
    for idx, gene_set in enumerate(gene_sets):
        ax = axes[idx]
        data = scores_long[scores_long['gene_set'] == gene_set]
        
        # Create boxplot
        sns.boxplot(
            data=data,
            x='group',
            y='module_score',
            palette=palette,
            ax=ax,
            width=0.6,
            linewidth=1.5,
        )
        
        # Overlay individual points if requested
        if show_points:
            sns.stripplot(
                data=data,
                x='group',
                y='module_score',
                color='black',
                alpha=0.6,
                size=5,
                ax=ax,
                jitter=True,
            )
        
        # Add statistical comparison (t-test between groups if exactly 2)
        unique_groups = data['group'].unique()
        if len(unique_groups) == 2:
            g1 = data[data['group'] == unique_groups[0]]['module_score']
            g2 = data[data['group'] == unique_groups[1]]['module_score']
            _, pval = stats.ttest_ind(g1, g2)
            
            # Add significance stars
            sig = _get_significance_symbol(pval)
            y_max = data['module_score'].max()
            y_range = data['module_score'].max() - data['module_score'].min()
            ax.text(0.5, y_max + 0.1 * y_range, sig, ha='center', fontsize=12)
            ax.set_title(f'{gene_set}\n(p={pval:.2e})', fontsize=11, fontweight='bold')
        else:
            ax.set_title(gene_set, fontsize=11, fontweight='bold')
        
        ax.set_xlabel('')
        ax.set_ylabel('Module Score', fontsize=10)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    
    # Hide unused subplots
    for idx in range(n_sets, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Module Scores by Condition', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def plot_module_violin(
    scores: pd.DataFrame,
    groups: pd.DataFrame,
    output: Path,
    palette: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (10, 6),
) -> None:
    """
    Create violin plots of module scores grouped by condition.
    
    Args:
        scores: Module scores (gene sets x samples).
        groups: Sample-to-group mapping.
        output: Output file path.
        palette: Color palette for groups.
        figsize: Figure dimensions.
    """
    # Convert scores to long format
    scores_long = scores.T.reset_index()
    scores_long.columns = ['sample_id'] + list(scores.index)
    scores_long = scores_long.melt(
        id_vars=['sample_id'],
        var_name='gene_set',
        value_name='module_score'
    )
    scores_long = scores_long.merge(groups, on='sample_id', how='left')
    scores_long['group'] = scores_long['group'].fillna('Unknown')
    
    gene_sets = scores.index.tolist()
    n_sets = len(gene_sets)
    n_cols = min(2, n_sets)
    n_rows = (n_sets + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0], figsize[1] * n_rows / 2))
    if n_sets == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    if palette is None:
        unique_groups = scores_long['group'].unique()
        colors = sns.color_palette("Set2", len(unique_groups))
        palette = dict(zip(unique_groups, colors))
    
    for idx, gene_set in enumerate(gene_sets):
        ax = axes[idx]
        data = scores_long[scores_long['gene_set'] == gene_set]
        
        sns.violinplot(
            data=data,
            x='group',
            y='module_score',
            palette=palette,
            ax=ax,
            inner='box',
            linewidth=1.5,
        )
        
        # Add individual points
        sns.stripplot(
            data=data,
            x='group',
            y='module_score',
            color='black',
            alpha=0.6,
            size=4,
            ax=ax,
            jitter=True,
        )
        
        ax.set_title(gene_set, fontsize=11, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Module Score', fontsize=10)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    
    for idx in range(n_sets, len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle('Module Scores by Condition', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def plot_module_heatmap(
    scores: pd.DataFrame,
    groups: Optional[pd.DataFrame],
    output: Path,
    cluster_samples: bool = True,
    cluster_gene_sets: bool = True,
) -> None:
    """
    Create a clustered heatmap of module scores with group annotations.
    
    Args:
        scores: Module scores (gene sets x samples).
        groups: Sample-to-group mapping (optional).
        output: Output file path.
        cluster_samples: Whether to cluster samples (columns).
        cluster_gene_sets: Whether to cluster gene sets (rows).
    """
    # Prepare group colors for annotation
    col_colors = None
    if groups is not None:
        sample_to_group = dict(zip(groups['sample_id'], groups['group']))
        unique_groups = groups['group'].unique()
        group_palette = dict(zip(unique_groups, sns.color_palette("Set2", len(unique_groups))))
        
        col_colors = pd.Series(
            [group_palette.get(sample_to_group.get(s, 'Unknown'), 'gray') for s in scores.columns],
            index=scores.columns,
            name='Group'
        )
    
    # Create clustermap
    g = sns.clustermap(
        scores,
        cmap='vlag',
        center=0,
        figsize=(max(10, len(scores.columns) * 0.5), max(6, len(scores) * 0.8)),
        col_cluster=cluster_samples,
        row_cluster=cluster_gene_sets,
        col_colors=col_colors,
        annot=True if scores.shape[1] <= 10 else False,
        fmt='.2f',
        linewidths=0.5,
        cbar_kws={'label': 'Module Score'},
        dendrogram_ratio=(0.1, 0.15),
    )
    
    g.ax_heatmap.set_xlabel('Samples', fontsize=12)
    g.ax_heatmap.set_ylabel('Gene Sets', fontsize=12)
    plt.suptitle('Module Score Heatmap', fontsize=14, fontweight='bold', y=1.02)
    
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def plot_combined_barplot(
    scores: pd.DataFrame,
    groups: pd.DataFrame,
    output: Path,
    palette: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (12, 6),
) -> None:
    """
    Create a grouped bar plot showing mean module scores per condition.
    
    This is useful when you have many samples and want to summarize by group.
    
    Args:
        scores: Module scores (gene sets x samples).
        groups: Sample-to-group mapping.
        output: Output file path.
        palette: Color palette for groups.
        figsize: Figure dimensions.
    """
    # Convert to long format and merge with groups
    scores_long = scores.T.reset_index()
    scores_long.columns = ['sample_id'] + list(scores.index)
    scores_long = scores_long.melt(
        id_vars=['sample_id'],
        var_name='gene_set',
        value_name='module_score'
    )
    scores_long = scores_long.merge(groups, on='sample_id', how='left')
    scores_long['group'] = scores_long['group'].fillna('Unknown')
    
    # Calculate mean and SEM per group
    summary = scores_long.groupby(['gene_set', 'group'])['module_score'].agg(['mean', 'std', 'count'])
    summary['sem'] = summary['std'] / np.sqrt(summary['count'])
    summary = summary.reset_index()
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=figsize)
    
    if palette is None:
        unique_groups = summary['group'].unique()
        colors = sns.color_palette("Set2", len(unique_groups))
        palette = dict(zip(unique_groups, colors))
    
    gene_sets = scores.index.tolist()
    unique_groups = summary['group'].unique()
    x = np.arange(len(gene_sets))
    width = 0.8 / len(unique_groups)
    
    for i, group in enumerate(unique_groups):
        group_data = summary[summary['group'] == group]
        # Ensure data is ordered by gene_set
        group_data = group_data.set_index('gene_set').loc[gene_sets].reset_index()
        
        bars = ax.bar(
            x + i * width - (len(unique_groups) - 1) * width / 2,
            group_data['mean'],
            width,
            label=group,
            color=palette.get(group, 'gray'),
            yerr=group_data['sem'],
            capsize=3,
        )
    
    ax.set_xlabel('Gene Set', fontsize=12)
    ax.set_ylabel('Mean Module Score ± SEM', fontsize=12)
    ax.set_title('Module Scores by Condition', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(gene_sets, rotation=45, ha='right')
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    ax.legend(title='Group', bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def plot_dotplot(
    scores: pd.DataFrame,
    groups: pd.DataFrame,
    output: Path,
    figsize: Tuple[int, int] = (10, 6),
) -> None:
    """
    Create a dot plot showing mean scores with size indicating magnitude.
    
    Similar to Seurat's DotPlot for visualizing module scores across conditions.
    
    Args:
        scores: Module scores (gene sets x samples).
        groups: Sample-to-group mapping.
        output: Output file path.
        figsize: Figure dimensions.
    """
    # Calculate mean scores per group
    scores_long = scores.T.reset_index()
    scores_long.columns = ['sample_id'] + list(scores.index)
    scores_long = scores_long.melt(
        id_vars=['sample_id'],
        var_name='gene_set',
        value_name='module_score'
    )
    scores_long = scores_long.merge(groups, on='sample_id', how='left')
    scores_long['group'] = scores_long['group'].fillna('Unknown')
    
    # Calculate mean per group
    mean_scores = scores_long.groupby(['gene_set', 'group'])['module_score'].mean().reset_index()
    
    # Create pivot for plotting
    pivot = mean_scores.pivot(index='gene_set', columns='group', values='module_score')
    
    fig, ax = plt.subplots(figsize=figsize)
    
    gene_sets = pivot.index.tolist()
    groups_list = pivot.columns.tolist()
    
    # Create dot plot
    for i, gene_set in enumerate(gene_sets):
        for j, group in enumerate(groups_list):
            score = pivot.loc[gene_set, group]
            # Size proportional to absolute score
            size = abs(score) * 200 + 50
            # Color based on direction
            color = 'red' if score > 0 else 'blue'
            alpha = min(abs(score) / pivot.abs().values.max(), 1.0)
            
            ax.scatter(j, i, s=size, c=color, alpha=max(0.3, alpha), edgecolors='black', linewidths=0.5)
    
    ax.set_xticks(range(len(groups_list)))
    ax.set_xticklabels(groups_list, rotation=45, ha='right')
    ax.set_yticks(range(len(gene_sets)))
    ax.set_yticklabels(gene_sets)
    
    ax.set_xlabel('Group', fontsize=12)
    ax.set_ylabel('Gene Set', fontsize=12)
    ax.set_title('Module Scores Dot Plot\n(Red: positive, Blue: negative, Size: magnitude)', 
                 fontsize=12, fontweight='bold')
    
    # Add grid
    ax.set_axisbelow(True)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()


def _get_significance_symbol(pval: float) -> str:
    """Convert p-value to significance symbol."""
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return 'ns'


# =============================================================================
# Command Line Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute module scores and generate publication-ready visualizations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # With explicit group file
    python module_score_visualization.py \\
        --vst expression_vst.txt \\
        --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \\
        --sample-groups groups.txt \\
        --output-prefix results/module_scores

    # Auto-infer groups from sample names
    python module_score_visualization.py \\
        --vst expression_vst.txt \\
        --gene-lists 2c_DOWN.txt 2c_UP.txt \\
        --infer-groups \\
        --output-prefix results/module_scores

    # Custom visualization types
    python module_score_visualization.py \\
        --vst expression_vst.txt \\
        --gene-lists 2c_DOWN.txt 2c_UP.txt \\
        --sample-groups groups.txt \\
        --plots boxplot violin heatmap \\
        --output-prefix results/module_scores
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--vst",
        required=True,
        type=Path,
        help="Path to VST expression matrix (TSV, genes as rows)."
    )
    parser.add_argument(
        "--gene-lists",
        required=True,
        nargs="+",
        type=Path,
        help="Paths to gene list files (one gene per line)."
    )
    
    # Gene column
    parser.add_argument(
        "--gene-column",
        default="external_gene_name",
        help="Column name containing gene identifiers (default: external_gene_name)."
    )
    
    # Sample groups
    group_args = parser.add_mutually_exclusive_group()
    group_args.add_argument(
        "--sample-groups",
        type=Path,
        help="TSV file with sample_id and group columns."
    )
    group_args.add_argument(
        "--infer-groups",
        action="store_true",
        help="Attempt to infer groups from sample names."
    )
    
    # Scoring method
    parser.add_argument(
        "--method",
        choices=["seurat-like", "mean"],
        default="seurat-like",
        help="Scoring method (default: seurat-like)."
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=24,
        help="Number of expression bins for control sampling (default: 24)."
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)."
    )
    
    # Output
    parser.add_argument(
        "--output-prefix",
        default="module_scores",
        type=str,
        help="Prefix for output files (default: module_scores)."
    )
    parser.add_argument(
        "--plots",
        nargs="+",
        choices=["boxplot", "violin", "heatmap", "barplot", "dotplot", "all"],
        default=["all"],
        help="Plot types to generate (default: all)."
    )
    parser.add_argument(
        "--format",
        choices=["png", "pdf", "svg"],
        default="png",
        help="Output format for plots (default: png)."
    )
    
    return parser.parse_args()


def main() -> None:
    """Main entry point."""
    args = parse_args()
    
    print("=" * 60)
    print("Module Score Visualization for Bulk RNA-seq")
    print("=" * 60)
    
    # Load expression matrix
    print(f"\n📂 Loading expression matrix: {args.vst}")
    expr = pd.read_csv(args.vst, sep="\t")
    
    if args.gene_column not in expr.columns:
        raise ValueError(f"Column '{args.gene_column}' not found in VST file.")
    
    expr = expr.set_index(args.gene_column)
    
    # Handle duplicates
    if expr.index.duplicated().any():
        n_dup = expr.index.duplicated().sum()
        print(f"  ⚠️  Removing {n_dup} duplicate genes (keeping first)")
        expr = expr[~expr.index.duplicated(keep='first')]
    
    # Keep only numeric columns (expression data)
    numeric_cols = expr.select_dtypes(include=[np.number]).columns.tolist()
    expr = expr[numeric_cols]
    
    print(f"  Matrix size: {expr.shape[0]} genes × {expr.shape[1]} samples")
    print(f"  Samples: {', '.join(expr.columns[:5])}{'...' if len(expr.columns) > 5 else ''}")
    
    # Load gene lists
    print(f"\n📂 Loading {len(args.gene_lists)} gene lists:")
    gene_sets = {}
    for path in args.gene_lists:
        genes = read_gene_list(path)
        gene_sets[path.stem] = genes
        print(f"  - {path.stem}: {len(genes)} genes")
    
    # Handle sample groups
    print("\n📊 Sample group assignment:")
    if args.sample_groups:
        groups = read_sample_groups(args.sample_groups)
        print(f"  Loaded from: {args.sample_groups}")
    elif args.infer_groups:
        groups = infer_groups_from_sample_names(expr.columns.tolist())
    else:
        # Create a default group (all samples in one group)
        groups = pd.DataFrame({
            'sample_id': expr.columns,
            'group': ['All'] * len(expr.columns)
        })
        print("  No groups specified - using single group 'All'")
    
    # Report group distribution
    group_counts = groups.groupby('group').size()
    print(f"  Groups: {dict(group_counts)}")
    
    # Compute module scores
    print(f"\n🔬 Computing module scores (method: {args.method})...")
    scores = compute_module_scores(
        expr,
        gene_sets,
        method=args.method,
        bins=args.bins,
        seed=args.seed,
    )
    
    # Summary statistics
    print(f"\n📈 Score summary:")
    print(f"  Mean ± SD: {scores.values.mean():.3f} ± {scores.values.std():.3f}")
    print(f"  Range: [{scores.values.min():.3f}, {scores.values.max():.3f}]")
    
    # Save scores
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    scores_file = output_prefix.with_suffix('.tsv')
    scores.to_csv(scores_file, sep='\t')
    print(f"\n✓ Scores saved to: {scores_file}")
    
    # Generate plots
    plot_types = args.plots if 'all' not in args.plots else ['boxplot', 'violin', 'heatmap', 'barplot', 'dotplot']
    fmt = args.format
    
    print(f"\n🎨 Generating visualizations...")
    
    if 'boxplot' in plot_types:
        output = Path(f"{output_prefix}_boxplot.{fmt}")
        plot_module_boxplot(scores, groups, output)
        print(f"  ✓ Boxplot: {output}")
    
    if 'violin' in plot_types:
        output = Path(f"{output_prefix}_violin.{fmt}")
        plot_module_violin(scores, groups, output)
        print(f"  ✓ Violin plot: {output}")
    
    if 'heatmap' in plot_types:
        output = Path(f"{output_prefix}_heatmap.{fmt}")
        plot_module_heatmap(scores, groups, output)
        print(f"  ✓ Heatmap: {output}")
    
    if 'barplot' in plot_types:
        output = Path(f"{output_prefix}_barplot.{fmt}")
        plot_combined_barplot(scores, groups, output)
        print(f"  ✓ Bar plot: {output}")
    
    if 'dotplot' in plot_types:
        output = Path(f"{output_prefix}_dotplot.{fmt}")
        plot_dotplot(scores, groups, output)
        print(f"  ✓ Dot plot: {output}")
    
    print("\n" + "=" * 60)
    print("✅ Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
