#!/usr/bin/env Rscript
# =============================================================================
# 08_sensitivity_analysis.R - how robust are the differential-expression results?
#
# 1) Leave-one-out: refit the model without SRR1782761, the control sample that
#    dominates PC1 (size factor 0.49), and check whether each DE gene of the
#    full model is still supported. Same sign is always required:
#      replicated     padj < 0.05 without the outlier
#      consistent     p < 0.05 without the outlier (padj >= 0.05)
#      not supported  otherwise
#    Dropping a control lowers the power, so "consistent" still counts as support.
# 2) Published signature: Haberman et al. (2014) reported DUOX2 up-regulation
#    and APOA1 down-regulation in the ileum of the full RISK cohort; check the
#    direction of these markers in this 3-vs-3 subset (log2FC = MLE, unshrunk).
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel)
})
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "scripts"
source(file.path(script_dir, "common.R"))

OUTLIER <- "SRR1782761"
OUT <- proj_path("results/sensitivity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

cnt <- load_counts()
coldata <- load_coldata(colnames(cnt))

fit <- function(samples) {
  dds <- DESeqDataSetFromMatrix(cnt[, samples], coldata[samples, ], design = ~ sex + condition)
  dds <- DESeq(prefilter(dds), quiet = TRUE)
  list(dds = dds, res = results(dds, contrast = c("condition", "Crohn", "Control"), alpha = 0.05))
}
full <- fit(colnames(cnt))
loo  <- fit(setdiff(colnames(cnt), OUTLIER))

## ---- 1) Leave-one-out support for the full-model DE genes ----
full_df <- annotate_results(full$res)
sig <- subset(full_df, !is.na(padj) & padj < 0.05)
loo_df <- as.data.frame(loo$res)[match(sig$gene_id, rownames(loo$res)), ]

same_sign <- !is.na(loo_df$log2FoldChange) & sign(loo_df$log2FoldChange) == sign(sig$log2FoldChange)
support <- ifelse(same_sign & !is.na(loo_df$padj) & loo_df$padj < 0.05, "replicated",
           ifelse(same_sign & !is.na(loo_df$pvalue) & loo_df$pvalue < 0.05, "consistent",
                  "not supported"))

# How far the outlier sits from the other two controls (normalised counts, log2 ratio).
norm <- counts(full$dds, normalized = TRUE)
other_controls <- setdiff(rownames(coldata)[coldata$condition == "Control"], OUTLIER)
outlier_shift <- log2((norm[sig$gene_id, OUTLIER] + 1) /
                      (rowMeans(norm[sig$gene_id, other_controls, drop = FALSE]) + 1))

loo_table <- data.frame(
  gene_id = sig$gene_id,
  symbol = sig$symbol,
  biotype = sig$biotype,
  direction = ifelse(sig$log2FoldChange > 0, "up", "down"),
  log2FC_full = round(sig$log2FoldChange, 3),
  padj_full = signif(sig$padj, 3),
  log2FC_without_outlier = round(loo_df$log2FoldChange, 3),
  pvalue_without_outlier = signif(loo_df$pvalue, 3),
  padj_without_outlier = signif(loo_df$padj, 3),
  outlier_log2_vs_other_controls = round(outlier_shift, 2),
  support = support
)
write.csv(loo_table, file.path(OUT, paste0("leave_one_out_", OUTLIER, ".csv")), row.names = FALSE)
support_counts <- table(direction = loo_table$direction,
                        support = factor(support, levels = c("replicated", "consistent", "not supported")))

## ---- 2) PCA without the outlier ----
vsd <- vst(loo$dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = c("condition", "sex"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
ggsave(file.path(OUT, paste0("PCA_without_", OUTLIER, ".png")),
  ggplot(pca, aes(PC1, PC2, color = condition, shape = sex)) +
    geom_point(size = 4) + geom_text_repel(aes(label = name), size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", pv[1], "%")) + ylab(paste0("PC2: ", pv[2], "%")) +
    theme_bw() + ggtitle(paste("PCA (VST) without", OUTLIER)), width = 7, height = 5, dpi = 150)

## ---- 3) Published ileal signature (Haberman et al. 2014) ----
markers <- data.frame(
  symbol   = c("DUOX2", "DUOXA2", "APOA1"),
  expected = c("up", "up", "down"),
  note     = c("Haberman et al. 2014", "maturation partner of DUOX2", "Haberman et al. 2014")
)
ids <- full_df$gene_id[match(markers$symbol, full_df$symbol)]
pick <- function(res) as.data.frame(res)[match(ids, rownames(res)), ]
m_full <- pick(full$res)
m_loo <- pick(loo$res)
markers_table <- cbind(markers,
  log2FC_full = round(m_full$log2FoldChange, 2),
  pvalue_full = signif(m_full$pvalue, 2),
  padj_full = signif(m_full$padj, 2),
  log2FC_without_outlier = round(m_loo$log2FoldChange, 2),
  pvalue_without_outlier = signif(m_loo$pvalue, 2),
  direction_as_published = ifelse(markers$expected == "up", m_full$log2FoldChange > 0,
                                   m_full$log2FoldChange < 0)
)
write.csv(markers_table, file.path(OUT, "published_markers.csv"), row.names = FALSE)

## ---- 4) Summary ----
report <- c(
  "=== Sensitivity analysis: Crohn vs Control ===",
  sprintf("Full model (6 samples): %d DE genes (padj < 0.05)", nrow(sig)),
  sprintf("Without %s (5 samples): %d DE genes (padj < 0.05)", OUTLIER,
          sum(loo$res$padj < 0.05, na.rm = TRUE)),
  "",
  "Support for the full-model DE genes without the outlier:",
  capture.output(print(support_counts)),
  "",
  "Published markers (log2FC = MLE, unshrunk):",
  capture.output(print(markers_table, row.names = FALSE))
)
writeLines(report, file.path(OUT, "sensitivity_summary.txt"))
cat(report, sep = "\n")
