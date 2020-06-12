# STACAS - Functions

* `Run.STACAS`   Wrapper function to run the complete anchor-finding algorithm in STACAS. For more control over individual steps, see the functions below.

* `FindAnchors.STACAS`  Computes anchors between datasets for integration in Seurat

* `PlotAnchors.STACAS`  Generates distribution plots which can be useful to inspect the anchor distance distributions between datasets
     and select a threshold for anchor filtering

* `FilterAnchors.STACAS`  Filters integration anchors based on pairwise rPCA distance, calculated using `FindAnchors.STACAS`

* `SampleTree.STACAS`   Determine automatically a hierarchical tree for Seurat integration using anchors calculated in STACAS

Find more information, syntax and examples using the R help function e.g. `?FindAnchors.STACAS`

