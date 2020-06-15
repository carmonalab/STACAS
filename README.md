# STACAS: Sub-Type Anchoring Correction for Alignment in Seurat

<p align="center">
  <img height="80" src="docs/white.sq.png">
</p>

`STACAS` is a method for anchor identification, designed to be easily incorporated in Seurat pipelines for batch correction and integration of scRNA-seq datasets.

Find the installation instructions for the package below, and a vignette detailing its functions at [Tutorial (html)](https://carmonalab.github.io/STACAS/tutorial.html) and [Tutorial (repository)](https://gitlab.unil.ch/carmona/STACAS.demo)

A Docker image for STACAS is available on [DockerHub](https://hub.docker.com/repository/docker/mandrea1/stacas_demo)

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

<p align="center">
  <img height="60" src="docs/white.sq.png">
</p>
