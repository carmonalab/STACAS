test_that("SampleTree.STACAS comprehensive tree generation and parameter tests", {
  # Shared base setup with 3 datasets for tree generation
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 3,
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
      min.sample.size = 50,
      verbose = FALSE
    )
  })
  
  # Test 1: Returns valid tree structure (matrix for multiple datasets)
  tree1 <- SampleTree.STACAS(
    anchorset = anchors,
    plot = FALSE
  )
  expect_true(is.matrix(tree1))
  
  # Test 2: Different hclust.method values produce valid trees
  methods <- c("single", "complete", "average", "ward.D2")
  for (method in methods) {
    tree_method <- SampleTree.STACAS(
      anchorset = anchors,
      hclust.method = method,
      plot = FALSE
    )
    expect_true(is.matrix(tree_method))
  }
  
  # Test 3: Different usecol parameters (score vs dist.mean)
  tree_score <- SampleTree.STACAS(
    anchorset = anchors,
    usecol = "score",
    plot = FALSE
  )
  expect_true(is.matrix(tree_score))
  
  tree_dist <- SampleTree.STACAS(
    anchorset = anchors,
    usecol = "dist.mean",
    plot = FALSE
  )
  expect_true(is.matrix(tree_dist))
  
  # Test 4: Different method parameters (weight.sum vs counts)
  tree_weight <- SampleTree.STACAS(
    anchorset = anchors,
    method = "weight.sum",
    plot = FALSE
  )
  expect_true(is.matrix(tree_weight))
  
  tree_counts <- SampleTree.STACAS(
    anchorset = anchors,
    method = "counts",
    plot = FALSE
  )
  expect_true(is.matrix(tree_counts))
  
  # Test 5: Semi-supervised parameter with labeled data
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
  
  tree_semi_true <- SampleTree.STACAS(
    anchorset = anchors_semi,
    semisupervised = TRUE,
    plot = FALSE
  )
  expect_true(is.matrix(tree_semi_true))
  
  tree_semi_false <- SampleTree.STACAS(
    anchorset = anchors_semi,
    semisupervised = FALSE,
    plot = FALSE
  )
  expect_true(is.matrix(tree_semi_false))
  
  # Test 6: plot parameter (TRUE should still return tree matrix)
  tree_plot <- SampleTree.STACAS(
    anchorset = anchors,
    plot = TRUE
  )
  expect_true(is.matrix(tree_plot))
  
  # Test 7: With custom object names
  obj_names <- c("Sample1", "Sample2", "Sample3")
  tree_names <- SampleTree.STACAS(
    anchorset = anchors,
    obj.names = obj_names,
    plot = FALSE
  )
  expect_true(is.matrix(tree_names))
  
  # Test 8: Combination of multiple parameters
  tree_combo <- SampleTree.STACAS(
    anchorset = anchors,
    obj.names = c("DatasetA", "DatasetB", "DatasetC"),
    hclust.method = "ward.D2",
    usecol = "score",
    method = "weight.sum",
    semisupervised = TRUE,
    plot = FALSE
  )
  expect_true(is.matrix(tree_combo))
})
