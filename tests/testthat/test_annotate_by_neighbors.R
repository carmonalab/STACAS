test_that("annotate.by.neighbors propagates labels correctly with various parameters", {
  # Single shared setup for all test cases
  obj <- create_dummy_seurat(n_cells = 100, n_genes = 120, seed = 22)
  obj <- process_Seurat(obj, ndim = 15)
  
  # Test 1: Assigns labels using provided reference cells
  set.seed(22)
  ref.cells_1 <- sample(colnames(obj), 25)
  obj_test1 <- obj
  obj_test1$functional.cluster <- NA
  obj_test1@meta.data[ref.cells_1, "functional.cluster"] <- 
    obj_test1@meta.data[ref.cells_1, "cell_type"]
  
  result1 <- annotate.by.neighbors(
    obj_test1,
    ref.cells = ref.cells_1,
    labels.col = "functional.cluster",
    reduction = "pca",
    ndim = 10,
    k = 5,
    ncores = 1
  )
  
  expect_s4_class(result1, "Seurat")
  expect_false(any(is.na(result1@meta.data$functional.cluster)))
  expect_equal(
    as.character(result1@meta.data[ref.cells_1, "functional.cluster"]),
    as.character(obj_test1@meta.data[ref.cells_1, "cell_type"])
  )
  
  # Test 2: Preserves dimensions after propagation
  obj_test2 <- obj
  ref.cells_2 <- colnames(obj_test2)[1:40]
  obj_test2$cell_type_partial <- NA
  obj_test2@meta.data[ref.cells_2, "cell_type_partial"] <- 
    obj_test2@meta.data[ref.cells_2, "cell_type"]
  
  original_cells <- ncol(obj_test2)
  original_genes <- nrow(obj_test2)
  
  result2 <- annotate.by.neighbors(
    obj_test2,
    labels.col = "cell_type_partial",
    ndim = 8,
    k = 7,
    ncores = 1
  )
  
  expect_equal(ncol(result2), original_cells)
  expect_equal(nrow(result2), original_genes)
  expect_false(any(is.na(result2@meta.data$cell_type_partial)))
  
  # Test 3: Works across different parameter choices
  obj_test3 <- obj
  ref.cells_3 <- colnames(obj_test3)[seq(1, 100, by = 2)]
  obj_test3$functional.cluster <- NA
  obj_test3@meta.data[ref.cells_3, "functional.cluster"] <- 
    obj_test3@meta.data[ref.cells_3, "cell_type"]
  
  result3 <- annotate.by.neighbors(
    obj_test3,
    ref.cells = ref.cells_3,
    labels.col = "functional.cluster",
    ndim = 6,
    k = 3,
    ncores = 1
  )
  
  expect_s4_class(result3, "Seurat")
  expect_false(any(is.na(result3@meta.data$functional.cluster)))
})
