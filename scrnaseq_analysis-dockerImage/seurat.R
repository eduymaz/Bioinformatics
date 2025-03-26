library("optparse")
library(DropletUtils)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(scico)
library(scales)

option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.txt", 
              help="output file name ", metavar="character"),
  make_option(c("-t", "--t2g"), type="character", default="out.txt", 
              help="output file name", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


theme_set(theme_bw())


# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}


# Read matrix into R
res_mat <- read_count_output(opt$dir, name = "cells_x_genes")

dim(res_mat)

tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

bc_rank <- barcodeRanks(res_mat, lower = 10)


#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}


options(repr.plot.width=9, repr.plot.height=6)

fullPathPDF <- paste(opt$output, "/1_kneeplot.pdf", sep="")
pdf(file=fullPathPDF)
knee_plot(bc_rank)
dev.off()


res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
res_mat <- res_mat[Matrix::rowSums(res_mat) > 0,]
dim(res_mat)


tr2g <- read_tsv(opt$t2g, col_names = c("transcript", "gene", "gene_symbol")) %>%
  select(-transcript) %>%
  distinct()

# Convert from Ensembl gene ID to gene symbol
rownames(res_mat) <- tr2g$gene_symbol[match(rownames(res_mat), tr2g$gene)]
unique_matris <- res_mat[!duplicated(rownames(res_mat)), ] #TODO

#START SEURAT
seu <- CreateSeuratObject(unique_matris, min.cells = 3, min.features = 200)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

# Visualize QC metrics as a violin plot
options(repr.plot.width=12, repr.plot.height=6)
fullPathPDF <- paste(opt$output, "/2_violinplot.pdf", sep="")
pdf(file=fullPathPDF)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

options(repr.plot.width=9, repr.plot.height=6)
fullPathPDF <- paste(opt$output, "/3_ggplot1.pdf", sep="")
pdf(file=fullPathPDF)

ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA)) +
  geom_hex(bins = 100) +
  scale_fill_scico(palette = "devon", direction = -1, end = 0.9) +
  scale_x_log10(breaks = breaks_log(12)) + 
  scale_y_log10(breaks = breaks_log(12)) + annotation_logticks() +
  labs(x = "Total UMI counts", y = "Number of genes detected") +
  theme(panel.grid.minor = element_blank())

dev.off()


fullPathPDF <- paste(opt$output, "/4_ggplot2.pdf", sep="")
pdf(file=fullPathPDF)
ggplot(seu@meta.data, aes(nCount_RNA, percent.mt)) +
  geom_pointdensity() +
  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  labs(x = "Total UMI counts", y = "Percentage mitochondrial")
dev.off()

seu <- subset(seu, subset = percent.mt < 3)
seu <- NormalizeData(seu) %>% ScaleData()

fullPathPDF <- paste(opt$output, "/5_FeaturePlot.pdf", sep="")
pdf(file=fullPathPDF)
seu <- FindVariableFeatures(seu, nfeatures = 3000)
top10 <- head(VariableFeatures(seu), 10)
plot1 <- VariableFeaturePlot(seu, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

dev.off()


seu <- RunPCA(seu, verbose = FALSE, npcs = 20) # uses HVG by default
fullPathPDF <- paste(opt$output, "/6_elbow.pdf", sep="")
pdf(file=fullPathPDF)
ElbowPlot(seu, ndims = 20)
dev.off()


fullPathPDF <- paste(opt$output, "/7_pca_first.pdf", sep="")
pdf(file=fullPathPDF)
PCAPlot(seu)
dev.off()

seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu)

#PCA
fullPathPDF <- paste(opt$output, "/8_pca.pdf", sep="")
pdf(file=fullPathPDF)
PCAPlot(seu)
dev.off()

#tSNE

seu <- RunTSNE(seu, dims = 1:10)
fullPathPDF <- paste(opt$output, "/9_tsne.pdf", sep="")
pdf(file=fullPathPDF)
TSNEPlot(seu)
dev.off()

#UMAP

seu <- RunUMAP(seu, dims = 1:10, verbose = FALSE)
fullPathPDF <- paste(opt$output, "/10_UMAP.pdf", sep="")
pdf(file=fullPathPDF)
UMAPPlot(seu)
dev.off()


fullPathPDF <- paste(opt$output, "/11_UMAP_celltype.pdf", sep="")
pdf(file=fullPathPDF)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()