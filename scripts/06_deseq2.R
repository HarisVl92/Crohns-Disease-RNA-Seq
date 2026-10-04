#!/usr/bin/env Rscript
# =============================================================================
# 06_deseq2.R — Differential expression: Crohn vs Control (host ileal RNA-seq)
# Design: ~ sex + condition   (sex ως covariate: 1 Female στους Crohn -> ελέγχεται)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap)
  library(RColorBrewer); library(ggrepel)
})
PROJ <- "/home/haris/Desktop/Bioinformatics_Project"
OUT  <- file.path(PROJ, "results/DESeq2")
dir.create(file.path(OUT,"figures"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(OUT,"tables"),  showWarnings=FALSE, recursive=TRUE)
set.seed(1)

## ---- 1) Φόρτωση πίνακα featureCounts ----
fc <- read.delim(file.path(PROJ,"results/counts/counts.txt"), comment.char="#", check.names=FALSE)
rownames(fc) <- fc$Geneid
cnt <- as.matrix(fc[, 7:ncol(fc)])
colnames(cnt) <- sub("\\.sorted\\.bam$","", basename(colnames(cnt)))

## ---- 2) colData από το manifest ----
man <- read.delim(file.path(PROJ,"data/samples_manifest.tsv"), comment.char="#")
rownames(man) <- man$run_accession
man <- man[colnames(cnt), ]
stopifnot(identical(rownames(man), colnames(cnt)))
coldata <- data.frame(
  row.names = rownames(man),
  condition = factor(man$condition, levels=c("Control","Crohn")),
  sex       = factor(man$sex,       levels=c("Male","Female"))
)
cat("== colData ==\n"); print(coldata)

## ---- 3) DESeqDataSet + prefilter ----
dds <- DESeqDataSetFromMatrix(countData=cnt, colData=coldata, design = ~ sex + condition)
smallestGroup <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroup
cat(sprintf("Prefilter: %d -> %d genes\n", nrow(dds), sum(keep)))
dds <- dds[keep, ]

## ---- 4) DESeq ----
dds <- DESeq(dds)
saveRDS(dds, file.path(OUT,"dds.rds"))
sf <- sizeFactors(dds)
write.csv(data.frame(sample=names(sf), sizeFactor=sf), file.path(OUT,"tables/size_factors.csv"), row.names=FALSE)

## ---- 5) Results: Crohn vs Control ----
res <- results(dds, contrast=c("condition","Crohn","Control"), alpha=0.05)
resShrink <- lfcShrink(dds, coef="condition_Crohn_vs_Control", type="normal", res=res)  # apeglm απών -> normal

ann <- read.delim(file.path(PROJ,"data/reference/gene_annotation.tsv"), header=TRUE)
rownames(ann) <- ann$gene_id
res_df <- as.data.frame(resShrink); res_df$gene_id <- rownames(res_df)
res_df$symbol  <- ann[res_df$gene_id,"gene_name"]
res_df$biotype <- ann[res_df$gene_id,"biotype"]
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(OUT,"tables/DE_results_all.csv"), row.names=FALSE)
sig <- subset(res_df, !is.na(padj) & padj < 0.05)
write.csv(sig, file.path(OUT,"tables/DE_results_padj0.05.csv"), row.names=FALSE)
up <- subset(sig, log2FoldChange>0); dn <- subset(sig, log2FoldChange<0)
writeLines(unique(na.omit(up$symbol)),  file.path(OUT,"tables/genes_UP_symbols.txt"))
writeLines(unique(na.omit(dn$symbol)),  file.path(OUT,"tables/genes_DOWN_symbols.txt"))
writeLines(unique(na.omit(sig$symbol)), file.path(OUT,"tables/genes_SIG_symbols.txt"))

## ---- 6) Transforms + normalized counts ----
vsd  <- vst(dds, blind=FALSE)
norm <- counts(dds, normalized=TRUE)
write.csv(norm, file.path(OUT,"tables/normalized_counts.csv"))

## ---- 7) Διαγνωστικά plots ----
png(file.path(OUT,"figures/01_dispersion.png"),1000,800,res=130); plotDispEsts(dds); dev.off()

pcaData <- plotPCA(vsd, intgroup=c("condition","sex"), returnData=TRUE)
pv <- round(100*attr(pcaData,"percentVar"))
ggsave(file.path(OUT,"figures/02_PCA.png"),
  ggplot(pcaData, aes(PC1,PC2,color=condition,shape=sex)) +
    geom_point(size=4) + geom_text_repel(aes(label=name),size=3,show.legend=FALSE) +
    xlab(paste0("PC1: ",pv[1],"%")) + ylab(paste0("PC2: ",pv[2],"%")) +
    theme_bw() + ggtitle("PCA (VST) — Crohn vs Control"), width=7,height=5,dpi=150)

sampleDist <- dist(t(assay(vsd))); sdm <- as.matrix(sampleDist)
png(file.path(OUT,"figures/03_sample_distance.png"),1000,900,res=130)
pheatmap(sdm, clustering_distance_rows=sampleDist, clustering_distance_cols=sampleDist,
         col=colorRampPalette(rev(brewer.pal(9,"Blues")))(255),
         annotation_col=as.data.frame(coldata)); dev.off()

png(file.path(OUT,"figures/04_MA.png"),1000,800,res=130)
plotMA(resShrink, ylim=c(-5,5), main="MA (shrunk LFC) — Crohn vs Control"); dev.off()

vol <- res_df
vol$sig <- ifelse(!is.na(vol$padj)&vol$padj<0.05&abs(vol$log2FoldChange)>1,
                  ifelse(vol$log2FoldChange>0,"Up","Down"),"NS")
lab <- head(vol[order(vol$padj),],20)
ggsave(file.path(OUT,"figures/05_volcano.png"),
  ggplot(vol, aes(log2FoldChange,-log10(pvalue),color=sig)) +
    geom_point(alpha=0.6,size=1) +
    scale_color_manual(values=c(Up="#c0392b",Down="#2471a3",NS="grey70")) +
    geom_text_repel(data=lab, aes(label=symbol), size=3, max.overlaps=20, show.legend=FALSE) +
    geom_vline(xintercept=c(-1,1),linetype="dashed",color="grey50") +
    geom_hline(yintercept=-log10(0.05),linetype="dashed",color="grey50") +
    theme_bw() + ggtitle("Volcano — Crohn vs Control"), width=8,height=6,dpi=150)

topg <- head(rownames(sig), 40)
if(length(topg)>=2){
  mat <- assay(vsd)[topg,]; mat <- mat-rowMeans(mat)
  rn <- ann[topg,"gene_name"]; rn[is.na(rn)|rn==""] <- topg[is.na(rn)|rn==""]
  rownames(mat) <- rn
  png(file.path(OUT,"figures/06_top40_heatmap.png"),1000,1300,res=130)
  pheatmap(mat, annotation_col=as.data.frame(coldata), scale="none",
           main="Top 40 DE genes (VST, centered)"); dev.off()
}

## ---- 8) Έλεγχος φύλου (recorded vs expression) ----
sexgenes <- c(XIST="ENSG00000229807", RPS4Y1="ENSG00000129824",
              DDX3Y="ENSG00000067048", UTY="ENSG00000183878")
present <- sexgenes[sexgenes %in% rownames(norm)]
if(length(present)>0){
  sm <- norm[present,,drop=FALSE]; rownames(sm) <- names(present)
  write.csv(t(round(sm,1)), file.path(OUT,"tables/sex_check_normcounts.csv"))
  cat("== Sex-linked normalized counts ==\n"); print(round(t(sm),1))
}

## ---- 9) Summary ----
sink(file.path(OUT,"DESeq2_summary.txt"))
cat("=== DESeq2 — Crohn vs Control ===\nDesign: ~ sex + condition\n\n")
print(coldata); cat("\nGenes after prefilter:", nrow(dds), "\n")
cat("Size factors:\n"); print(round(sf,3))
cat("\n"); summary(res)
cat("\nSignificant padj<0.05:", nrow(sig), "| Up:", nrow(up), " Down:", nrow(dn), "\n\nTop 25:\n")
print(head(res_df[,c("symbol","baseMean","log2FoldChange","padj","biotype")],25))
sink()
writeLines(capture.output(sessionInfo()), file.path(OUT,"sessionInfo.txt"))
cat("DESeq2 DONE\n")
