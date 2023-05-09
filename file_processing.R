# This R script is intended for preprocessing the files downloaded from GSE64810

# Extract sample metadata and output as a csv file
library(rio)
library(GEOquery)
library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)


# gse <- rio::import("data/GSE64810_series_matrix.txt")
gse = getGEO(filename="data/GSE64810_series_matrix.txt")
metadata <- pData(gse)
rio::export(metadata,file = "data/metadata.csv", row.names = TRUE)


# Convert normalized counts to .csv format
norm <- rio::import("data/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt") %>%
  dplyr::rename(gene = V1)
rio::export(norm, "data/norm_counts.csv")

# Convert DE to .csv format
deg <- rio::import("data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt") %>%
  dplyr::rename(gene = V1)
rio::export(deg, "data/deg.csv")


# GSEA calculation
set.seed(1234)
deg <- rio::import("data/GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt") %>%
  dplyr::rename(gene = V1)
msigdb_gene_sets <- msigdbr(species = "Homo sapiens")
pathways <- msigdb_gene_sets %>%
  dplyr::select(gs_name, human_gene_symbol) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(human_gene_symbol)) %>%
  tibble::deframe()
ranked_genes <- deg %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::select(symbol, log2FoldChange) %>%
  tibble::deframe()
fgsea_res <- fgsea(pathways = pathways, 
                   stats = ranked_genes,
                   minSize = 15,
                   maxSize = 500) 
fgsea_res_sorted <- fgsea_res %>%
  dplyr::arrange(padj)
rio::export(fgsea_res_sorted, "data/fgsea_res_sorted.csv")

# head(ranked_genes)
# head(msigdb_gene_sets)
# head(pathways)
# head(fgsea_res)




