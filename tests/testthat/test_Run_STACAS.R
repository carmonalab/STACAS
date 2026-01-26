test_that("Run.STACAS wrapper function integration tests", {
  # Setup shared datasets
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 2,
    n_cells_per_dataset = 120,
    n_genes = 120,
    seed = 123
  )
  
  # Test 1: End-to-end pipeline completes successfully
  suppressMessages({
    integrated1 <- Run.STACAS(
      object.list = obj_list,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      k.weight = 10,
      min.sample.size = 50,
      verbose = FALSE
    )
  })
  
  expect_s4_class(integrated1, "Seurat")
  expect_true("integrated" %in% Seurat::Assays(integrated1))
  expect_equal(ncol(integrated1), sum(sapply(obj_list, ncol)))
  expect_true("pca" %in% Seurat::Reductions(integrated1))

})
