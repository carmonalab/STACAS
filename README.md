# STACAS: Accurate integration of single-cell transcriptomics data with cell type unbalance

<p align="center">
  <img height="80" src="docs/RSticker_STACAS.png">
</p>

[STACAS](https://github.com/carmonalab/STACAS) is a method for scRNA-seq integration, specifically designed to align scRNA-seq datasets that are composed of only partially overlapping cell populations or sub-types.
It is based on the [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) integration framework, and adds the following innovations:

* anchors are down-weighted based on their distance in reciprocal PCA space
* guide trees for pairwaise integration are constructed based on the 'centrality' of datasets, as measured by the sum of their re-weighted anchor scores
* Prior knowledge, given as cell labels, can be used by the algorithm to remove inconsistent anchors, and perform semi-supervised integration 

## Package Installation

To install STACAS directly from the Git repository, run the following code from within RStudio:

```r
if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/STACAS")
```

## STACAS basic usage
### Standard integration (more [here](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#standard-integration))
```r
library(STACAS)

# get the test dataset "pbmcsca" from SeuratData package
if (!requireNamespace("remotes")) install.packages("remotes")
if (!requireNamespace("SeuratData")) install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
InstallData("pbmcsca")
data("pbmcsca")

# Integrate scRNA-seq datasets generated with different methods/technologies
pbmcsca.integrated <- NormalizeData(pbmcsca) |>
    SplitObject(split.by = "Method")|>
    Run.STACAS()

pbmcsca.integrated <- RunUMAP(pbmcsca.integrated, dims = 1:30) 

# Visualize
DimPlot(pbmcsca.integrated, group.by = c("Method","CellType")) 
```

### Semi-supervised integration (more [here](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html#semi-supervised-integration))

```r
pbmcsca.semisup <- NormalizeData(pbmcsca) |>
    SplitObject(split.by = "Method")|>
    Run.STACAS(cell.labels = "CellType")

pbmcsca.semisup <- RunUMAP(pbmcsca.semisup, dims = 1:30) 
```

## STACAS integration DEMOS

Find a tutorial for `STACAS` in a complete Seurat integration pipeline at: [STACAS demo](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html) (code and instructions [here](https://github.com/carmonalab/STACAS.demo))

See also how `STACAS` compares to Seurat for the integration of heterogeneos data sets: [STACAS vs Seurat](https://carmonalab.github.io/STACAS.demo/Tcell.demo.html)


### Citation

Massimo Andreatta, Santiago J Carmona "STACAS: Sub-Type Anchor Correction for Alignment in Seurat to integrate single-cell RNA-seq data", *Bioinformatics* **(2020)** - https://doi.org/10.1093/bioinformatics/btaa755

<p align="center">
  <img height="60" src="docs/RSticker_STACAS.png">
</p>
