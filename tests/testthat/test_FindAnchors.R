test_that("FindAnchors.STACAS comprehensive anchor finding tests", {
  # Test 1: Returns valid AnchorSet with basic parameters
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 3,
    n_cells_per_dataset = 120,
    n_genes = 120,
    seed = 123
  )
  
  suppressMessages({
    anchors1 <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      min.sample.size = 10,
      verbose = FALSE
    )
  })
  
  expect_s4_class(anchors1, "AnchorSet")
  anchor_df <- slot(anchors1, "anchors")
  expected_cols <- c("cell1", "cell2", "dataset1", "dataset2", "score")
  expect_true(all(expected_cols %in% colnames(anchor_df)))
  expect_gt(nrow(anchor_df), 0)
  
  # Test 2: Filters by min.sample.size (removes small datasets)
  obj_list_size <- obj_list
  obj_list_size[[3]] <- obj_list_size[[3]][, 1:5]
  
  suppressMessages({
    anchors2 <- FindAnchors.STACAS(
      object.list = obj_list_size,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      min.sample.size = 50,
      verbose = FALSE
    )
  })
  
  anchor_df2 <- slot(anchors2, "anchors")
  datasets_in_anchors <- unique(c(anchor_df2$dataset1, anchor_df2$dataset2))
  expect_lte(length(datasets_in_anchors), 2)
  
  # Test 3: Parameter validation (alpha, label.confidence, anchor.coverage)
  expect_error(
    FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      alpha = 1.5,
      verbose = FALSE
    )
  )
  
  expect_error(
    FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      label.confidence = 1.5,
      verbose = FALSE
    )
  )
  
  expect_error(
    FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      anchor.coverage = 1.5,
      verbose = FALSE
    )
  )
  
  # Test 4: Handles dims parameter (single and vector)
  suppressMessages({
    anchors4a <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      verbose = FALSE
    )
  })
  expect_s4_class(anchors4a, "AnchorSet")
  
  suppressMessages({
    anchors4b <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = c(1, 2, 3, 4, 5),
      k.anchor = 3,
      k.score = 5,
      verbose = FALSE
    )
  })
  expect_s4_class(anchors4b, "AnchorSet")
  
  # Test 5: Semi-supervised with complete cell labels
  suppressMessages({
    anchors5 <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      cell.labels = "cell_type",
      label.confidence = 0.8,
      min.sample.size = 10,
      verbose = FALSE
    )
  })
  
  anchor_df5 <- slot(anchors5, "anchors")
  expect_true("Consistent" %in% colnames(anchor_df5))
  
  # Test 6: Semi-supervised with partial NA cell labels
  # Introduce NA values in cell labels
  obj_list_na <- obj_list[1:2]
  for (i in seq_along(obj_list_na)) {
    n_cells <- ncol(obj_list_na[[i]])
    set.seed(22)
    na_indices <- sample(1:n_cells, size = round(n_cells * 0.3))
    obj_list_na[[i]]@meta.data[na_indices, "cell_type"] <- NA
  }
  
  suppressMessages({
    anchors6 <- FindAnchors.STACAS(
      object.list = obj_list_na,
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      cell.labels = "cell_type",
      label.confidence = 0.8,
      verbose = FALSE
    )
  })
  
  expect_s4_class(anchors6, "AnchorSet")
  expect_true("Consistent" %in% colnames(slot(anchors6, "anchors")))
  
  # Test 7: scale.data parameter (FALSE and TRUE)
  suppressMessages({
    anchors7a <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      scale.data = FALSE,
      verbose = FALSE
    )
  })
  expect_s4_class(anchors7a, "AnchorSet")
  expect_gt(nrow(slot(anchors7a, "anchors")), 0)
  
  suppressMessages({
    anchors7b <- FindAnchors.STACAS(
      object.list = obj_list[1:2],
      anchor.features = 50,
      dims = 5,
      k.anchor = 3,
      k.score = 5,
      scale.data = TRUE,
      verbose = FALSE
    )
  })
  expect_s4_class(anchors7b, "AnchorSet")
  expect_gt(nrow(slot(anchors7b, "anchors")), 0)
  
})
