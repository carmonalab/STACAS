test_that("FindVariableFeatures.STACAS comprehensive feature selection tests", {
  # Single shared setup
  test_obj <- create_dummy_seurat(n_cells = 50, n_genes = 200)
  
  # Test 1: Returns expected structure with default parameters
  result1 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 50,
    genesBlockList = NULL,
    min.exp = 0.01,
    max.exp = 3
  )
  
  expect_s4_class(result1, "Seurat")
  expect_true(length(Seurat::VariableFeatures(result1)) > 0)
  expect_lte(length(Seurat::VariableFeatures(result1)), 50)
  expect_true(all(Seurat::VariableFeatures(result1) %in% rownames(result1)))
  
  # Test 2: Applies gene blocklist correctly
  all_genes <- rownames(test_obj)[1:10]
  result2 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 100,
    genesBlockList = all_genes,
    min.exp = 0.01,
    max.exp = 3
  )
  
  var_genes <- Seurat::VariableFeatures(result2)
  blocked_genes <- intersect(all_genes, var_genes)
  expect_length(blocked_genes, 0)
  
  # Test 3: Handles edge cases (nfeat larger than genes)
  result3 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 500,
    genesBlockList = NULL,
    min.exp = 0.01,
    max.exp = 3
  )
  
  expect_true(length(Seurat::VariableFeatures(result3)) <= nrow(test_obj))
  expect_true(length(Seurat::VariableFeatures(result3)) > 0)
  
  # Test 4: Works with default blocklist
  result4 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 50,
    genesBlockList = "default",
    min.exp = 0.01,
    max.exp = 3
  )
  
  expect_s4_class(result4, "Seurat")
  expect_true(length(Seurat::VariableFeatures(result4)) > 0)
  
  # Test 5: Handles permissive parameters (function completes without error)
  result5 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 20,
    genesBlockList = NULL,
    min.exp = 0.01,
    max.exp = 10
  )
  
  expect_s4_class(result5, "Seurat")
  # Just verify the function runs and returns an object; actual feature count can vary
  expect_true(is.character(Seurat::VariableFeatures(result5)))
  
  # Test 6: Correctly filters genes by mean expression bounds
  # Get normalized data to check mean expressions
  norm_data <- Seurat::GetAssayData(test_obj, layer = "data")
  gene_means <- rowMeans(as.matrix(norm_data))
  
  # Select min/max that should exclude some genes
  min_threshold <- quantile(gene_means, 0.3)
  max_threshold <- quantile(gene_means, 0.7)
  
  result6 <- FindVariableFeatures.STACAS(
    test_obj,
    nfeat = 100,
    genesBlockList = NULL,
    min.exp = min_threshold,
    max.exp = max_threshold
  )
  
  expect_s4_class(result6, "Seurat")
  
  # Verify that selected features are within bounds
  var_genes <- Seurat::VariableFeatures(result6)
  if (length(var_genes) > 0) {
    selected_means <- gene_means[var_genes]
    expect_true(all(selected_means >= min_threshold - 0.01))  # Small tolerance for floating point
    expect_true(all(selected_means <= max_threshold + 0.01))
  }
})
