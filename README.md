# STACAS: Sub-Type Anchoring Correction for Alignment in Seurat

<p align="center">
  <img height="80" src="docs/white.sq.png">
</p>

`STACAS` is a method for anchor identification, designed to be easily incorporated in Seurat pipelines for batch correction and integration of scRNA-seq datasets ([Andreatta & Carmona, Bioinformatics 2020](http://dx.doi.org/10.1093/bioinformatics/btaa755)).

`STACAS` is ideal to align scRNA-seq datasets that are composed of only partially overlapping cell populations or sub-types, where other methods tend to under-perform.

To see `STACAS` in action on larger scale integration tasks towards the construction of reference T cell maps in cancer and infection, please refer to this paper: [Andreatta et al BioRxiv 2020](https://doi.org/10.1101/2020.06.23.166546).

Multi-study integrated atlases:

* tumor-infiltrating T cell atlas: http://TILatlas.unil.ch (Seurat object available at https://doi.org/10.6084/m9.figshare.12478571)

* viral infection CD8 T cell atlas: http://virusTcellAtlas.unil.ch/ (Seurat object available at https://doi.org/10.6084/m9.figshare.12489518)


Find the installation instructions for the package below, and a vignette detailing its functions at [Tutorial (html)](https://carmonalab.github.io/STACAS/tutorial.html) and [Tutorial (repository)](https://gitlab.unil.ch/carmona/STACAS.demo)

If you prefer to avoid installing R packages, you can run `STACAS` in Docker.
A ready-to-use Docker image with usage instructions is available on [DockerHub](https://hub.docker.com/repository/docker/mandrea1/stacas_demo)

### Package Installation

To install STACAS directly from the Git repository, run the following code from within RStudio:

```
if (!requireNamespace("remotes")) install.packages("remotes")
library(remotes)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Seurat", quietly = TRUE)) {
   BiocManager::install('multtest')
   install.packages("Seurat")
}

remotes::install_github("carmonalab/STACAS")
```

### Test the package

Load sample data and test the package:
```
library(STACAS)

data(STACAS.sampledata)

STACAS.anchors <- Run.STACAS(STACAS.sampledata)
```

### STACAS integration TUTORIAL

Find a tutorial for `STACAS` in a complete Seurat integration pipeline at: [STACAS Tutorial](https://carmonalab.github.io/STACAS/tutorial.html)

To run the code of the tutorial on your machine, download the demo repository: [STACAS Tutorial repo](https://gitlab.unil.ch/carmona/STACAS.demo) or obtain a [Docker image](https://hub.docker.com/repository/docker/mandrea1/stacas_demo) with all dependencies pre-installed.

### Documentation

See a description of the functions implemented in STACAS at: [STACAS functions](docs/functions.md)

### Citation

Massimo Andreatta, Santiago J Carmona, STACAS: Sub-Type Anchor Correction for Alignment in Seurat to integrate single-cell RNA-seq data, Bioinformatics, 2020, btaa755, https://doi.org/10.1093/bioinformatics/btaa755

<p align="center">
  <img height="60" src="docs/white.sq.png">
</p>
