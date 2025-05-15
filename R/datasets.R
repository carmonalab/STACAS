#' Sample dataset to test STACAS installation
#'
#' A Seurat object containing single-cell transcriptomes
#' (scRNA-seq) for 50 cells and 20729 genes.
#' Single-cell UMI counts were normalized using a standard log-normalization:
#' counts for each cell were divided by the total counts for that cell and
#' multiplied by 10,000, then natural-log transformed using `log1p`. \cr\cr
#' This a subsample of 25 predicted B cells and 25 predicted NK cells from
#' the large scRNA-seq PBMC dataset published
#' by Hao et al. (\doi{10.1016/j.cell.2021.04.048}) and
#' available as UMI counts at
#' \url{https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat}
#' 
#' @format A sparse matrix of 50 cells and 20729 genes.
#' @source \doi{10.1016/j.cell.2021.04.048}
"sampleObj"

#' Genes blocklists for excluding HVGs
#'
#' A list of gene signatures, including cycling, heat-shock response,
#' mitochondrial and risobomal genes, interferon response; for mouse and
#' human. Derived from the SignatuR R package:
#' \url{https://github.com/carmonalab/SignatuR}
#' 
#' @format A list of gene signatures
#' @source \url{https://github.com/carmonalab/SignatuR}
"genes.blocklist"

#' Standardized gene list from ENSEMBL (human)
#'
#' A reference of stable gene names for Homo Sapiens
#' 
#' @format A dataframe of ENSEMBL and gene symbols
#' @source \url{https://www.ensembl.org/Homo_sapiens/Info/Index}
"EnsemblGeneTable.Hs"

#' Standardized gene list from ENSEMBL (mouse)
#'
#' A reference of stable gene names for Mus Musculus
#' 
#' @format A dataframe of ENSEMBL and gene symbols
#' @source \url{https://www.ensembl.org/Mus_musculus/Info/Index}
"EnsemblGeneTable.Mm"