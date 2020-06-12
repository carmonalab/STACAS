# STACAS: Sub-Type Anchoring Correction for Alignment in Seurat
`STACAS` is a method for anchor identification, meant to be incorporated in Seurat pipelines for batch correction and integration of scRNA-seq datasets.

Find the installation instructions for the package below, and a vignette detailing its functions at [STACAS.demo](https://gitlab.unil.ch/carmona/STACAS.demo)

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

remotes::install_git("https://gitlab.unil.ch/carmona/STACAS.git")
```

### Test the package

Load sample data and test the package:
```
library(STACAS)

data(STACAS.sampledata)

STACAS.anchors <- Run.STACAS(STACAS.sampledata)
```

You can find a tutorial for `STACAS` in a complete Seurat integration pipeline at [STACAS.demo](https://gitlab.unil.ch/carmona/STACAS.demo)
