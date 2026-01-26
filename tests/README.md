# STACAS Test Suite

## Overview

Comprehensive unit tests for all STACAS exported functions using testthat Edition 3. Tests use lightweight, sparse matrix test data for CRAN compliance.

## Quick Start

```r
devtools::test()                          # Run all tests
devtools::test(filter = "FindAnchors")    # Run specific test file
```

## Test Strategy

### Single test_that per file
Each test file contains ONE comprehensive `test_that()` block with multiple subtests. This ensures test data is generated only once per file, improving efficiency.

### Test Coverage

| Function | Subtests | Key Features Tested |
|----------|----------|---------------------|
| FindVariableFeatures.STACAS | 6 | Blocklist filtering, expression bounds, mean-based removal |
| FindAnchors.STACAS | 7 | AnchorSet structure, min.sample.size, semi-supervised, NA labels, scale.data |
| SampleTree.STACAS | 8 | Tree generation, usecol, method, semisupervised, plot, hclust.method |
| IntegrateData.STACAS | 6 | Integration output, custom assays, k.weight, semi-supervised |
| Run.STACAS | 8 | Wrapper behavior, PCA computation, assay propagation |
| PlotAnchors.STACAS | 6 | ggplot output, anchor.coverage parameter |
| annotate.by.neighbors | 4 | Label transfer, ref.cells, ndim, k, ncores |

**Total**: 45+ subtests across 7 core functions

## Test Data Helpers

### `create_dummy_seurat()`
Creates a single Seurat object with sparse matrices (~80% zeros).

```r
create_dummy_seurat(
  n_cells = 100,
  n_genes = 500,
  cell_types = c("TypeA", "TypeB"),
  batch_id = "batch1",
  dataset_id = "dataset1",
  seed = 123
)
```

**Features**: Poisson counts, cell-type markers, log-normalized data, realistic metadata

### `create_dummy_datasets_with_batch()`
Creates multiple Seurat objects with batch effects (applied to all datasets except first).

```r
create_dummy_datasets_with_batch(
  n_datasets = 3,
  n_cells_per_dataset = 80,
  n_genes = 150,
  batch_effect_strength = 0.3,
  seed = 123
)
```

**Batch effects**: Gene-specific multiplicative and additive noise simulating technical variation

### `process_Seurat()`
Safely prepares objects for PCA with automatic dimension capping.

```r
process_Seurat(obj, n_features = 100, ndim = 30)
```

**Safety**: Caps PCs to `min(ndim, ncol-1, length(VariableFeatures))`

## Test Patterns

### Standard Test Block
```r
test_that("Function comprehensive tests", {
  # Shared setup
  obj_list <- create_dummy_datasets_with_batch(
    n_datasets = 3, n_cells_per_dataset = 120, n_genes = 120, seed = 123
  )
  
  # Subtest 1: Basic functionality
  result1 <- my_function(obj_list, param = value)
  expect_s4_class(result1, "ExpectedClass")
  
  # Subtest 2: Parameter validation
  expect_error(my_function(obj_list, invalid_param = 999))
  
  # Subtest 3: Edge cases
  # ... more subtests
})
```

### Parameter Testing
Each test file validates:
- Default parameter behavior
- Parameter value ranges (e.g., alpha 0-1)
- Optional parameters (e.g., cell.labels, scale.data)
- Edge cases (NA labels, small datasets)

## CRAN Compliance

- ✅ Execution time: ~20-30 seconds (< 60s limit)
- ✅ No external dependencies or downloads
- ✅ All data generated in-memory with sparse matrices
- ✅ Reproducible (seed-based)

## Performance Characteristics

| Dataset Size | Generation Time | Typical Use |
|--------------|-----------------|-------------|
| 50 cells, 100 genes | ~20ms | Fast parameter tests |
| 120 cells, 120 genes | ~100ms | Standard function tests |
| 3×80 cells, 150 genes | ~300ms | Integration workflows |

## Maintenance

### Adding Tests
1. Add subtests to existing `test_that()` block in relevant file
2. Test new parameters with both valid and invalid values
3. Include edge cases (NA values, empty inputs, boundary conditions)
4. Suppress output with `suppressMessages()` for clean test runs

### Extending Test Data
Edit `helper-create_test_data.R` while maintaining:
- Sparse matrix format
- Small data sizes (CRAN compliance)
- Reproducibility (seeds)
- Realistic distributions

---

**Last Updated**: January 26, 2026  
**Test Count**: 45+ subtests  
**Coverage**: All exported functions
