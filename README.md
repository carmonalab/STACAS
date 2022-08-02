# STACAS: Sub-Type Anchoring Correction for Alignment in Seurat

<p align="center">
  <img height="80" src="docs/white.sq.png">
</p>

`STACAS` is a method for anchor identification, designed to be easily incorporated in Seurat pipelines for batch correction and integration of scRNA-seq datasets ([Andreatta & Carmona, Bioinformatics 2020](http://dx.doi.org/10.1093/bioinformatics/btaa755)).

`STACAS` is ideal to align scRNA-seq datasets that are composed of only partially overlapping cell populations or sub-types, where other methods tend to under-perform.

To see `STACAS` in action on larger scale integration tasks towards the construction of reference T cell maps in cancer and infection, please refer to these papers: [Andreatta et al. Nat Comm 2021](https://www.nature.com/articles/s41467-021-23324-4) and [Andreatta et al. eLife 2022](https://elifesciences.org/articles/76339)

Multi-study integrated atlases:

* tumor-infiltrating T cell atlas: https://spica.unil.ch/refs/TIL (Seurat object available at https://doi.org/10.6084/m9.figshare.12478571)

* viral infection CD8 T cell atlas: https://spica.unil.ch/refs/viral-CD8-T (Seurat object available at https://doi.org/10.6084/m9.figshare.12489518)

* viral infection CD4 T cell atlas: https://spica.unil.ch/refs/viral-CD4-T (Seurat object available at https://doi.org/10.6084/m9.figshare.16592693)


Find the installation instructions for the package below, and a vignette detailing its functions at [Demo (html)](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html) and [Demo (repository)](https://github.com/carmonalab/STACAS.demo)


### Package Installation

To install STACAS directly from the Git repository, run the following code from within RStudio:

```
if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)

remotes::install_github("carmonalab/STACAS")
```

### STACAS basic usage
```
library(STACAS)

# get the test dataset "pbmcsca" from SeuratData package
if (!requireNamespace("remotes")) install.packages("remotes")
if (!requireNamespace("SeuratData")) install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
InstallData("pbmcsca")
data("pbmcsca")

# Integrate scRNA-seq datasets generated with different methods/technologies
pbmcsca <- pbmcsca |> NormalizeData() |> SplitObject(split.by = "Method") |> Run.STACAS() |> RunUMAP(dims = 1:30) 

# Visualize
DimPlot(pbmcsca, group.by = c("Method","CellType")) 
```

### STACAS integration DEMOS

Find a tutorial for `STACAS` in a complete Seurat integration pipeline at: [STACAS demo](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html)

See also how `STACAS` compares to Seurat for the integration of heterogeneos data sets: [STACAS vs Seurat](https://carmonalab.github.io/STACAS.demo/Tcell.demo.html)


### Citation

Massimo Andreatta, Santiago J Carmona, STACAS: Sub-Type Anchor Correction for Alignment in Seurat to integrate single-cell RNA-seq data, Bioinformatics, 2020, btaa755, https://doi.org/10.1093/bioinformatics/btaa755

<p align="center">
  <img height="60" src="docs/white.sq.png">
</p>
