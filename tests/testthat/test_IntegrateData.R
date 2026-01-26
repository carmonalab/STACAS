test_that("IntegrateData.STACAS comprehensive integration tests", {
  # Shared base setup with 2 datasets
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 2,
    n_cells_per_dataset = 120,
    n_genes = 120,
    seed = 123
  )
  suppressMessages({
    anchors <- FindAnchors.STACAS(
      object.list = obj_list,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      verbose = FALSE
    )
  })
  
  # Test 1: Returns integrated Seurat object with default assay
  suppressMessages({
    integrated1 <- IntegrateData.STACAS(
      anchorset = anchors,
      dims = 5,
      k.weight = 10,
      verbose = FALSE
    )
  })
  
  expect_s4_class(integrated1, "Seurat")
  expect_true("integrated" %in% Seurat::Assays(integrated1))
  integrated_data <- Seurat::GetAssayData(integrated1, assay = "integrated", layer = "data")
  expect_gt(ncol(integrated_data), 0)
  expect_gt(nrow(integrated_data), 0)
  
  # Test 2: With custom assay name
  suppressMessages({
    integrated2 <- IntegrateData.STACAS(
      anchorset = anchors,
      new.assay.name = "batch_corrected",
      dims = 5,
      k.weight = 10,
      verbose = FALSE
    )
  })
  
  expect_true("batch_corrected" %in% Seurat::Assays(integrated2))
  
  # Test 3: Different k.weight values
  suppressMessages({
    integrated3a <- IntegrateData.STACAS(
      anchorset = anchors,
      dims = 5,
      k.weight = 5,
      verbose = FALSE
    )
  })
  expect_s4_class(integrated3a, "Seurat")
  
  suppressMessages({
    integrated3b <- IntegrateData.STACAS(
      anchorset = anchors,
      dims = 5,
      k.weight = "max",
      verbose = FALSE
    )
  })
  expect_s4_class(integrated3b, "Seurat")
  
  # Test 4: With custom sample tree
    custom_tree <- SampleTree.STACAS(
      anchorset = anchors,
      plot = FALSE
    )
  
  suppressMessages({
    integrated4 <- IntegrateData.STACAS(
      anchorset = anchors,
      dims = 5,
      k.weight = 10,
      sample.tree = custom_tree,
      verbose = FALSE
    )
  })
  expect_s4_class(integrated4, "Seurat")
  
  # Test 5: Preserves cell identities and metadata
  suppressMessages({
    integrated5 <- IntegrateData.STACAS(
      anchorset = anchors,
      dims = 5,
      k.weight = 10,
      verbose = FALSE
    )
  })
  
  expect_true("batch" %in% colnames(integrated5@meta.data))
  expect_true("cell_type" %in% colnames(integrated5@meta.data))
  
  # Test 6: Semi-supervised integration
  suppressMessages({
    anchors_semi <- FindAnchors.STACAS(
      object.list = obj_list,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      cell.labels = "cell_type",
      label.confidence = 0.8,
      verbose = FALSE
    )
  })
  
  suppressMessages({
    integrated6a <- IntegrateData.STACAS(
      anchorset = anchors_semi,
      dims = 5,
      k.weight = 10,
      semisupervised = TRUE,
      verbose = FALSE
    )
  })
  expect_s4_class(integrated6a, "Seurat")
  
  suppressMessages({
    integrated6b <- IntegrateData.STACAS(
      anchorset = anchors_semi,
      dims = 5,
      k.weight = 10,
      semisupervised = FALSE,
      verbose = FALSE
    )
  })
  expect_s4_class(integrated6b, "Seurat")
})
