# Helper functions to create dummy scRNAseq data for testing STACAS
# These functions create light-weight, reproducible test datasets

#' Create a minimal dummy Seurat object for testing
#'
#' @param n_cells Number of cells (default: 100)
#' @param n_genes Number of genes (default: 500)
#' @param cell_types Vector of cell types to sample from
#' @param seed Random seed for reproducibility
#' @param batch_id Batch identifier for metadata (default: "Batch1")
#' @param dataset_id Dataset identifier for metadata (default: NULL)
#'
#' @return A Seurat object with minimal data
#'
#' @keywords internal
create_dummy_seurat <- function(n_cells = 100,
                                n_genes = 500,
                                cell_types = c("CellType_A", "CellType_B", "Celltype_C"),
                                seed = 123,
                                batch_id = "Batch1",
                                dataset_id = NULL) {
  set.seed(seed)
  
  # Create sparse expression matrix using Matrix package
  # Use realistic poisson-derived counts (zero-inflated)
  counts_data <- integer(0)
  i_indices <- integer(0)
  j_indices <- integer(0)
  
  # Simulate counts with realistic sparsity (~80%)
  for (i in seq_len(n_genes)) {
    # Some genes more highly expressed
    lambda <- if (i %% 5 == 0) 2 else 0.5
    gene_counts <- rpois(n_cells, lambda = lambda)
    
    # Only store non-zero values for sparse matrix
    non_zero <- which(gene_counts > 0)
    if (length(non_zero) > 0) {
      counts_data <- c(counts_data, gene_counts[non_zero])
      i_indices <- c(i_indices, rep(i, length(non_zero)))
      j_indices <- c(j_indices, non_zero)
    }
  }
  
  # Create sparse matrix
  counts <- Matrix::sparseMatrix(
    i = i_indices,
    j = j_indices,
    x = counts_data,
    dims = c(n_genes, n_cells)
  )
  
  # Add some cell-type specific genes for biological signal
  cell_assignments <- sample(cell_types, n_cells, replace = TRUE)
  
  # Dynamically assign marker genes to each cell type
  n_cell_types <- length(cell_types)
  for (gene_idx in seq_len(n_genes)) {
    # Distribute marker genes across cell types
    # Each cell type gets genes at different intervals
    for (ct_idx in seq_len(n_cell_types)) {
      # Use different modulo values for different cell types
      marker_interval <- 10 + (ct_idx - 1) * 5
      if (gene_idx %% marker_interval == 0) {
        idx <- which(cell_assignments == cell_types[ct_idx])
        if (length(idx) > 0) {
          extra_counts <- rpois(length(idx), lambda = 1.5)
          counts[gene_idx, idx] <- counts[gene_idx, idx] + extra_counts
        }
      }
    }
  }
  
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- paste0("Cell", seq_len(n_cells))
  
  # Create metadata
  meta_data <- data.frame(
    cell_type = cell_assignments,
    batch = batch_id,
    row.names = colnames(counts)
  )
  
  if (!is.null(dataset_id)) {
    meta_data$dataset <- dataset_id
  }
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = meta_data
  )
  
  # Normalize data
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  
  return(seurat_obj)
}


#' Create dummy datasets with batch effects
#'
#' @param n_datasets Number of datasets to create
#' @param n_cells_per_dataset Number of cells per dataset
#' @param n_genes Number of genes (shared across datasets)
#' @param batch_effect_strength Strength of batch effect (0-1)
#' @param cell_types Cell types to use
#' @param seed Random seed
#'
#' @return A list of Seurat objects with batch effects
#'
#' @keywords internal
create_dummy_datasets_with_batch <- function(
    n_datasets = 2,
    n_cells_per_dataset = 100,
    n_genes = 500,
    batch_effect_strength = 0.3,
    cell_types = c("CellType_A", "CellType_B", "CellType_C"),
    seed = 123) {

  datasets <- list()
  
  for (ds_idx in seq_len(n_datasets)) {
    # fix number of cells
    myseed <- ds_idx + seed
    set.seed(myseed)
    n_cells <- rnorm(1,
                     mean = n_cells_per_dataset,
                     sd = n_cells_per_dataset/10)
    # Create base dataset using create_dummy_seurat
    base_obj <- create_dummy_seurat(
      n_cells = n_cells,
      n_genes = n_genes,
      cell_types = cell_types,
      seed = myseed,
      batch_id = paste0("Batch_", ds_idx)
    )
    
    # First dataset (ds_idx == 1) is used as-is without batch effects
    if (ds_idx != 1) {
      # Apply batch effects to subsequent datasets
      # Get counts as sparse matrix
      counts <- Seurat::GetAssayData(base_obj, layer = "counts")
      
      # Generate batch-specific scaling and offset
      batch_offset <- rnorm(n_genes, mean = 0, sd = batch_effect_strength)
      batch_scale <- exp(rnorm(n_genes, mean = 0, sd = batch_effect_strength * 0.5))
      
      # Apply batch effects while maintaining sparsity
      # Convert to matrix temporarily for operations
      counts_matrix <- as.matrix(counts)
      for (i in seq_len(n_genes)) {
        counts_matrix[i, ] <- pmax(0, round(counts_matrix[i, ] * batch_scale[i] + batch_offset[i]))
      }
      
      # Convert back to sparse matrix
      counts_sparse <- Matrix::Matrix(counts_matrix, sparse = TRUE)
      rownames(counts_sparse) <- rownames(counts)
      colnames(counts_sparse) <- colnames(counts)
      
      # Update cell names to include dataset index
      new_cell_names <- paste0("Cell_", ds_idx, "_", seq_len(n_cells))
      colnames(counts_sparse) <- new_cell_names
      rownames(base_obj@meta.data) <- new_cell_names
      
      # Create new Seurat object with batch-affected counts
      suppressWarnings({
        seurat_obj <- Seurat::CreateSeuratObject(
          counts = counts_sparse,
          meta.data = base_obj@meta.data
        )
      })
      
      # Normalize
      seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
    } else {
      seurat_obj <- base_obj
    }
    
    datasets[[paste0("Dataset_", ds_idx)]] <- seurat_obj
  }
  
  return(datasets)
}


#' Verify basic structure of a Seurat object
#'
#' @param seurat_obj A Seurat object
#'
#' @return Logical: TRUE if object has expected structure
#'
#' @keywords internal
verify_seurat_structure <- function(seurat_obj) {
  has_counts <- !is.null(SeuratObject::Layers(seurat_obj, search = "counts"))
  has_metadata <- nrow(seurat_obj@meta.data) > 0
  has_cells <- ncol(seurat_obj) > 0
  has_genes <- nrow(seurat_obj) > 0
  
  all(c(has_counts, has_metadata, has_cells, has_genes))
}


#' Minimal preprocessing pipeline for PCA-ready Seurat objects
#'
#' Runs NormalizeData (optional), FindVariableFeatures, ScaleData, and RunPCA
#' with sensible defaults for lightweight test datasets.
#'
#' @param seurat_obj A Seurat object
#' @param nfeatures Number of variable features to select (default: 1000)
#' @param ndim Number of principal components to compute (default: 20)
#' @param normalize Whether to run NormalizeData before feature selection
#'
#' @return A Seurat object with PCA embeddings
#'
#' @keywords internal
process_Seurat <- function(seurat_obj,
                           nfeatures = 1000,
                           ndim = 20,
                           normalize = FALSE) {
  if (normalize) {
    seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  }
  
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    nfeatures = nfeatures,
    verbose = FALSE
  )
  
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  
  # Cap number of PCs to available dimensions to avoid failures on tiny test data
  max_pcs <- min(
    ndim,
    max(1L, ncol(seurat_obj) - 1L),
    length(Seurat::VariableFeatures(seurat_obj))
  )
  
  seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    npcs = max_pcs,
    verbose = FALSE
  )
  
  return(seurat_obj)
}
