#!/usr/bin/env Rscript
.libPaths("/storage/vannostrand/software/rnaseq/venv/lib/R/library")

args <- commandArgs(trailingOnly=TRUE)
usage <- "
deseq.R - A script for using DESeq2 to detect differentially expressed genes.

USAGE:
    Rscript deseq.R COUNT_MATRIX DESIGN_MATRIX [separator]

    Positional Parameters:
      COUNT_MATRIX   -  Path to a text file contains gene level count data
      DESIGN_MATRIX  -  Path to a text file contains the design of the experiment

    Optional Parameters:
      separator      -  A separator for count and design matrix, default: tab

"
if (length(args)==0) {
  cat(usage)
  quit(save = "no", status = 1, runLast = FALSE)
} else if (length(args)==1) {
  if (args[1]=="-h" | args[1]=='--help') {
    cat(usage)
    quit(save = "no", status = 1, runLast = FALSE)
  } else {
    cat("At least two argument must be supplied, see the usage below", usage, sep='\n')
    quit(save = "no", status = 1, runLast = FALSE)
  }
} else if (length(args)==2) {
  count_matrix <- args[1]
  design_matrix <- args[2]
  separator <- "\t"
} else if (length(args)==3) {
  count_matrix_file <- args[1]
  design_matrix <- args[2]
  separator <- args[3]
} else {
  cat("Too many arguments provided, see the usage below.", usage, sep='\n')
  quit(save = "no", status = 1, runLast = FALSE)
}

suppressPackageStartupMessages({
  library("DESeq2")
  library("pheatmap")
  library("purrr")
  library("gridExtra")
  library("vsn")
  library("ashr")
  library("RColorBrewer")
  library("ggplot2")
})

plot_basic_qc <- function (count_data, output='BasicQC.png') {
  png(file=output, width = 1000, height = 800, res=200)
  barplot(colSums(count_data)*1e-6, names=colnames(count_data), col='lightblue',
          xlab="Samples", ylab="Library size (millions)", main="Basic QC")
  dev.off()
}

plot_ma <- function(dds, result, output='MA_plot.png') {
  png(file=output, width = 1800, height = 600, res=200)
  par(mfrow=c(1, 4), mar=c(4, 4, 2, 1))
  xlim <- c(1, 1e5); ylim <- c(-3, 3)
  result_LFC <- lfcShrink(dds, coef = 2, type = "apeglm")
  result_Normal <- lfcShrink(dds, coef = 2, type = "normal")
  result_Ash <- lfcShrink(dds, coef = 2, type = "ashr")
  plotMA(result, xlim=xlim, ylim=ylim, main='No-Shrinkage')
  plotMA(result_LFC, xlim=xlim, ylim=ylim, main='APEGLM')
  plotMA(result_Normal, xlim=xlim, ylim=ylim, main='Normal')
  plotMA(result_Ash, xlim=xlim, ylim=ylim, main='Ash')
  dev.off()
}

plot_count <- function(dds, result) {
  png(file='minimum_counts.png', width = 1000, height = 800, res=200)
  plotCounts(dds, gene = which.min(result$padj), intgroup = "group")
  dev.off()
}

write_result_to_csv <- function(ordered_result, significant_result) {
  write.csv(as.data.frame(ordered_result), file="result.csv")
  write.csv(as.data.frame(significant_result), file="significant_result.csv")
}

data_transform <- function(ntd, vsd, rld) {
  png("data_transformation.png", width = 1800, height = 600, res=200)
  mean_sd_ntd <- meanSdPlot(assay(ntd))
  mean_sd_vsd <- meanSdPlot(assay(vsd))
  mean_sd_rld <- meanSdPlot(assay(rld))
  plot_list <- ls(pattern = "mean_sd_") %>% map(~eval(as.name(.))$gg)
  grid.arrange(arrangeGrob(grobs = plot_list, ncol=3,
                           top="Data Transformation (Normal|VST|RLOG)"))
  dev.off()
}

count_matrix_heatmap <- function(dds, ntd, vsd, rld) {
  png("count_matrix_heatmap.png", width = 1800, height = 600, res=200)
  select <- order(rowMeans(counts(dds, normalized=TRUE)),
                  decreasing=TRUE)[1:30]
  df <- as.data.frame(colData(dds)["group"])
  heatmap_ntd <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                          cluster_cols = F, annotation_col = df)
  heatmap_vsd <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                          cluster_cols = F, annotation_col = df)
  heatmap_rld <- pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
                          cluster_cols = F, annotation_col = df)
  plot_list <- ls(pattern = "heatmap_") %>% map(~eval(as.name(.))[[4]])
  grid.arrange(arrangeGrob(grobs = plot_list, ncol=3,
                           top="Data Transformation (Normal|VST|RLOG)"))
  dev.off()
}

plot_distance <- function(data) {
  sample_distance <- dist(t(assay(data)))
  distance_matrix <- as.matrix(sample_distance)
  rownames(distance_matrix) <- data$group
  colnames(distance_matrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  img <- pheatmap(distance_matrix,
                  clustering_distance_rows=sample_distance,
                  clustering_distance_cols=sample_distance,
                  col=colors)
  return(img)
}
plot_sample_to_sample_distance <- function(ntd, vsd, rld) {
  png(file="sample_to_sample_distance.png", width = 1800, height = 1000, res=200)
  heatmap_ntd <- plot_distance(ntd)
  heatmap_vsd <- plot_distance(vsd)
  heatmap_rld <- plot_distance(rld)
  plot_list <- ls(pattern = "heatmap_") %>% map(~eval(as.name(.))[[4]])
  grid.arrange(arrangeGrob(grobs = plot_list, ncol=3, top="Data Transformation (Normal|VST|RLOG)"))
  dev.off()
}

plot_principal_component <- function(ntd, vsd, rld) {
  png(file="principal_component.png", width = 1800, height = 600, res=200)
  p1 <- plotPCA(ntd, intgroup="group", returnData=F) + ggtitle("NTD") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  p2 <- plotPCA(vsd, intgroup="group", returnData=F) + ggtitle("VST") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
  p3 <- plotPCA(rld, intgroup="group", returnData=F) + ggtitle("RLD") + theme(aspect.ratio = 1,plot.title = element_text(hjust = 0.5))
  grid.arrange(p1, p2, p3, ncol=3)
  dev.off()
}

plot_despersion <- function(dds) {
  png(file="dispersion_estimates.png", width = 1200, height = 1200, res=200)
  plotDispEsts(dds)
  dev.off()
}

plot_independent_filtering <- function(result) {
  png(file="independent_filtering.png", width = 800, height = 800, res=200)
  plot(metadata(result)$filterNumRej,
       type="b", ylab="Number of rejections",
       xlab="Quantiles of filter")
  lines(metadata(result)$lo.fit, col="red")
  abline(v=metadata(result)$filterTheta)
  dev.off()
}

gene_differential_express <- function(count_matrix, design_matrix, minimum_count=0, sep='\t',
                                      fdr_alpha=0.1, abs_fold_change=1) {
  count_data <- read.csv(count_matrix, header=T, row.names=1, sep=sep)
  design_matrix <- read.csv(design_matrix, header=T, row.names=1, sep=sep)
  setwd(dirname(count_matrix))
  col_data <- DataFrame(group=factor(design_matrix$group))
  dds <- DESeqDataSetFromMatrix(count_data, col_data, formula(~group))
  dds <- DESeq(dds)
  keep <- rowSums(counts(dds)) >= minimum_count
  dds <- dds[keep, ]
  result <- results(dds)
  ordered_result <- result[order(result$padj), ]
  significant_result <- ordered_result[!is.na(ordered_result$padj) &
                                         ordered_result$padj < fdr_alpha &
                                         abs(ordered_result$log2FoldChange) >= abs_fold_change, ]
  write_result_to_csv(ordered_result, significant_result)
  plot_basic_qc(count_data)

  ntd <- normTransform(dds)
  vsd <- vst(dds, blind = F)
  rld <- rlog(dds, blind = F)

  plot_ma(dds, result)
  plot_principal_component(ntd, vsd, rld)
  data_transform(ntd, vsd, rld)
  count_matrix_heatmap(dds, ntd, vsd, rld)
  plot_sample_to_sample_distance(ntd, vsd, rld)
  plot_count(dds, result)
  plot_despersion(dds)
  plot_independent_filtering(result)
}


gene_differential_express(count_matrix, design_matrix)
quit(save = "no", status = 0, runLast = FALSE)