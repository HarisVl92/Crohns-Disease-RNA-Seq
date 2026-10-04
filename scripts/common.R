# =============================================================================
# common.R - shared helpers for the R steps (06_deseq2.R, 08_sensitivity_analysis.R).
#
# The calling script defines `script_dir` (its own folder) before sourcing this
# file. PROJ, the repository root, can be overridden from the environment.
# =============================================================================
PROJ <- Sys.getenv("PROJ", unset = normalizePath(file.path(script_dir, "..")))
proj_path <- function(...) file.path(PROJ, ...)

# featureCounts matrix (genes x samples); the repository ships a gzipped copy.
load_counts <- function() {
  f <- proj_path("results/counts/counts.txt")
  if (!file.exists(f)) f <- paste0(f, ".gz")
  fc <- read.delim(f, comment.char = "#", check.names = FALSE)
  cnt <- as.matrix(fc[, 7:ncol(fc)])
  rownames(cnt) <- fc$Geneid
  colnames(cnt) <- sub("\\.sorted\\.bam$", "", basename(colnames(cnt)))
  cnt
}

# Sample information from the manifest, in the column order of the count matrix.
load_coldata <- function(samples) {
  man <- read.delim(proj_path("data/samples_manifest.tsv"), comment.char = "#")
  rownames(man) <- man$run_accession
  man <- man[samples, ]
  stopifnot(identical(rownames(man), samples))
  data.frame(
    row.names = rownames(man),
    condition = factor(man$condition, levels = c("Control", "Crohn")),
    sex       = factor(man$sex, levels = c("Male", "Female"))
  )
}

# Keep genes with >= 10 reads in at least as many samples as the smallest group.
prefilter <- function(dds) {
  min_samples <- min(table(dds$condition))
  dds[rowSums(DESeq2::counts(dds) >= 10) >= min_samples, ]
}

# apeglm shrinkage when the package is installed, else DESeq2's 'normal' estimator.
shrinkage_type <- function() {
  if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
}

# Results as a data frame with gene symbols and biotypes, ordered by adjusted p-value.
annotate_results <- function(res) {
  ann <- read.delim(proj_path("data/reference/gene_annotation.tsv"))
  rownames(ann) <- ann$gene_id
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df$symbol  <- ann[df$gene_id, "gene_name"]
  df$biotype <- ann[df$gene_id, "biotype"]
  df[order(df$padj), ]
}
