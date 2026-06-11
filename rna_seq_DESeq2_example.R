# RNA-seq analysis in R with DESeq2
# Usage:
#   Rscript rna_seq_DESeq2_example.R --counts example_counts.tsv --coldata example_coldata.tsv --output results

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

# Parse simple CLI arguments
parse_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) {
    return(args[idx + 1])
  }
  return(default)
}

counts_file <- parse_arg("--counts", "example_counts.tsv")
coldata_file <- parse_arg("--coldata", "example_coldata.tsv")
output_dir <- parse_arg("--output", "rna_seq_results")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Reading counts from:", counts_file, "\n")
cat("Reading metadata from:", coldata_file, "\n")

counts <- read.delim(counts_file, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
coldata <- read.delim(coldata_file, sep = "\t", stringsAsFactors = FALSE)

# Basic checks
if (!all(c("sample", "condition") %in% names(coldata))) {
  stop("coldata must contain columns: sample, condition")
}

if (!all(coldata$sample %in% colnames(counts))) {
  stop("Some samples in coldata are not present in the counts matrix")
}

# Reorder counts to match metadata
counts <- counts[, coldata$sample, drop = FALSE]

# DESeq2 object
# IMPORTANT: condition should be a factor with the reference level first
coldata$condition <- factor(coldata$condition)
coldata$condition <- relevel(coldata$condition, ref = levels(coldata$condition)[1])

# If your experiment has more than 2 groups, replace 'condition' accordingly.
# For a two-group comparison, DESeq2 will compare the second level against the first.

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Keep genes with at least 10 counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DE analysis
cat("Running DESeq2...\n")
# Robust fallback for small toy datasets
if (ncol(dds) >= 4) {
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  dds <- nbinomWaldTest(dds)
} else {
  dds <- DESeq(dds)
}

# Differential expression results
res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df %>%
  arrange(pvalue)

write.table(res_df,
            file = file.path(output_dir, "differential_expression_results.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Top DE genes (optional)
top_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

write.table(top_sig,
            file = file.path(output_dir, "significant_genes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# PCA plot
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 13) +
  ggtitle("PCA of RNA-seq samples")

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 7, height = 6, dpi = 150)

# Volcano plot
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "yes", "no")

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  theme_bw(base_size = 13) +
  labs(title = "Volcano plot", x = "log2 fold change", y = "-log10(p-value)") +
  scale_color_manual(values = c("no" = "grey60", "yes" = "firebrick")) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "volcano_plot.png"), volcano, width = 7, height = 6, dpi = 150)

cat("\nDone. Files created in:", output_dir, "\n")
cat("- differential_expression_results.tsv\n")
cat("- significant_genes.tsv\n")
cat("- pca_plot.png\n")
cat("- volcano_plot.png\n")
