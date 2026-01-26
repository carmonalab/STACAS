test_that("PlotAnchors.STACAS comprehensive visualization and parameter tests", {
  # Shared base setup with 2 datasets
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 2,
    n_cells_per_dataset = 120,
    n_genes = 120,
    seed = 123
  )
  
  suppressMessages({
    anchors_2 <- FindAnchors.STACAS(
      object.list = obj_list,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      verbose = FALSE
    )
  })
  
  # Test 1: Returns list of ggplot objects
  plots1 <- PlotAnchors.STACAS(
    ref.anchors = anchors_2,
    anchor.coverage = 0.5
  )
  
  expect_type(plots1, "list")
  expect_gt(length(plots1), 0)
  expect_true(all(sapply(plots1, function(x) inherits(x, "ggplot"))))
  
  # Test 2: With custom object names
  obj_names <- c("Sample_A", "Sample_B")
  plots2 <- PlotAnchors.STACAS(
    ref.anchors = anchors_2,
    obj.names = obj_names,
    anchor.coverage = 0.5
  )
  
  expect_type(plots2, "list")
  expect_gt(length(plots2), 0)
  
  # Test 3: Validates anchor.coverage parameter (errors)
  expect_error(
    PlotAnchors.STACAS(
      ref.anchors = anchors_2,
      anchor.coverage = 1.5
    )
  )
  
  expect_error(
    PlotAnchors.STACAS(
      ref.anchors = anchors_2,
      anchor.coverage = -0.1
    )
  )
  
  # Test 4: Different anchor.coverage values
  for (coverage in c(0.25, 0.5, 0.75)) {
    plots_cov <- PlotAnchors.STACAS(
      ref.anchors = anchors_2,
      anchor.coverage = coverage
    )
    expect_type(plots_cov, "list")
  }
  
})
