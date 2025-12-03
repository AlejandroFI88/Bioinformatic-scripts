# GSVA Module Score Script - README

## Quick Start

### Installation

```bash
# Install required package
pip install gseapy
```

### Usage

```bash
python3 module_score_gsva.py \
    --vst your_vst_matrix.txt \
    --gene-lists gene_set1.txt gene_set2.txt \
    --gene-column external_gene_name \
    --output scores_gsva.tsv \
    --heatmap scores_gsva_heatmap.png
```

## What is GSVA?

GSVA (Gene Set Variation Analysis) is a robust, publication-quality method for calculating gene set enrichment scores. It's more sophisticated than simple mean-based approaches.

**Reference:** Hänzelmann et al., BMC Bioinformatics (2013)

## Comparison with Simple Method

| Method | `module_score.py` | `module_score_gsva.py` |
|--------|------------------|----------------------|
| Algorithm | Mean(signature) - Mean(all genes) | GSVA (KS-like statistic) |
| For publication | ⚠️ Needs justification | ✅ Widely accepted |
| Speed | 🚀 Fast | 🐢 Slower |
| Robustness | ✅ Good for bulk | ✅✅ Very robust |

## When to Use GSVA

**Use GSVA (`module_score_gsva.py`) when:**
- ✅ Preparing manuscript for publication
- ✅ Reviewers request established methods
- ✅ Maximum statistical robustness needed

**Use simple method (`module_score.py`) when:**
- ✅ Exploratory analysis
- ✅ Quick results needed
- ✅ Easy interpretation is priority

## Citation

If using GSVA in your publication:

> Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.

## Troubleshooting

**Error: "gseapy is not installed"**
```bash
pip install gseapy
```

**Error: "No module named 'pandas'"**
```bash
pip install pandas seaborn matplotlib
```

## Files

- `module_score.py` - Original simple mean-based method
- `module_score_gsva.py` - New GSVA-based method (this script)
