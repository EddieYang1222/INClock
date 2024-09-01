#' Estimate dispersion parameters
#'
#' @param counts An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells. It can be sparse.
#' @param manifold The computed manifold from compute_manifold. It should be a matrix
#' of the same dimension as counts.
#' @param cell.types A vector with the same length as the columns of counts. Each value should be
#' a string corresponding to the cell type of each cell in counts.
#' @param ages A vector with the same length as the columns of counts. Each value should be
#' a string corresponding to the age group of each cell in counts.
#' @param model Model to be used for estimating dispersion parameters. Can be either "cCV" (constant CV),
#' "cFF" (constant Fano factor), or "cVar" (constant variance).
#' @param preprocess Whether the expression count matrix should be preprocessed,
#' where all genes with zero counts will removed. Default is TRUE.
#' @param ncores Number of cores to use. Default is 1.
#' @param size.factor A vector of cell size normalization factors.
#' Default uses mean library size normalization.
#' @param cell.types.to.use A vector of unique cell types to be used in the estimation process.
#' Default is all cell types.
#' @param ages.to.use A vector of unique age groups to be used in the estimation process.
#' Default is all age groups.
#' @param cell.types.cutoff Minimum number of cells in a cell type x age group. Default is 10.
#'
#' @return A named list for all dispersion parameters. Each element is a table
#' of dispersion parameters for a cell type. The columns in each table contain
#' dispersion parameters for each age group.
#' @importFrom Matrix rowSums colSums
#' @export
#'
#' @examples
#' # Loading test data from Tabula Muris Senis
#' library(TabulaMurisSenisData)
#' tms_marrow <- TabulaMurisSenisDroplet(tissues = "Marrow")
#' tms_marrow_counts <- tms_marrow$Marrow@assays@data$counts
#' rownames(tms_marrow_counts) <- rownames(tms_marrow$Marrow)
#' colnames(tms_marrow_counts) <- colnames(tms_marrow$Marrow)
#' tms_marrow_counts <- as(tms_marrow_counts, "dgCMatrix")
#'
#' tms_marrow_ages <- which(tms_marrow$Marrow$age %in% c("3m", "30m"))
#' tms_marrow_cell_types <- which(tms_marrow$Marrow$cell_ontology_class %in% c(
#' 'hematopoietic precursor cell', 'megakaryocyte-erythroid progenitor cell',
#' 'precursor B cell'))
#' tms_marrow_counts <- tms_marrow_counts[, intersect(tms_marrow_ages, tms_marrow_cell_types)]
#'
#' tms_marrow_manifold_neighbor <- compute_manifold(tms_marrow_counts, method = "neighbor")
#'
#'
#' tms_marrow_dispersions <- estimate_dispersion(tms_marrow_counts, tms_marrow_manifold_neighbor,
#'   tms_marrow$Marrow$cell_ontology_class[intersect(tms_marrow_ages, tms_marrow_cell_types)],
#'   tms_marrow$Marrow$age[intersect(tms_marrow_ages, tms_marrow_cell_types)],
#'   model = "cCV", ncores = 4, cell.types.cutoff = 10)
#'
estimate_dispersion <- function(counts, manifold, cell.types, ages, model = "cCV", preprocess = TRUE, ncores = 1,
                                size.factor = NULL, cell.types.to.use = NULL, ages.to.use = NULL,
                                cell.types.cutoff = 10) {
  if (preprocess) {
    message("Preprocessing the count matrix")
    message("The initial matrix size is ", nrow(counts), " genes and ", ncol(counts), " cells.")
    counts <- counts[rowSums(counts) != 0, ]
    message(
      "After removing genes with zero counts, the new matrix size is ",
      nrow(counts), " genes and ", ncol(counts), " cells."
    )
  }

  if (is.null(size.factor)) {
    message("Calculating normalization factors")
    size.factor <- colSums(counts) / mean(colSums(counts))
  }

  # Keep only specific cell types and ages, if necessary
  if (!is.null(cell.types.to.use)) {
    message("Selecting specific cell types")
    cell.types.index <- which(cell.types %in% c(cell.types.to.use))
    counts <- counts[, cell.types.index]
    manifold <- manifold[, cell.types.index]
    cell.types <- cell.types[cell.types.index]
    ages <- ages[cell.types.index]
  }

  if (!is.null(ages.to.use)) {
    message("Selecting specific age groups")
    ages.index <- which(ages %in% c(ages.to.use))
    counts <- counts[, ages.index]
    manifold <- manifold[, ages.index]
    cell.types <- cell.types[ages.index]
    ages <- ages[ages.index]
  }

  cell.types <- factor(cell.types)
  ages <- factor(ages)
  cell.types.levels <- levels(cell.types)
  ages.levels <- levels(ages)
  cell.types.nums <- as.numeric(cell.types)
  ages.nums <- as.numeric(ages)

  # Estimate dispersion parameters
  dispersion.list <- list()
  for (i in 1:max(cell.types.nums)) {
    for (j in 1:max(ages.nums)) {
      message("Estimating dispersion parameters for ", ages.levels[j], " ", cell.types.levels[i])
      index <- rep(0, ncol(counts))
      for (k in 1:ncol(counts)) {
        index[k] <- cell.types.nums[k] == i && ages.nums[k] == j
      }
      if (sum(index) >= cell.types.cutoff) {
        # Subset counts and normalize to match original manifold fitting
        celltype.counts.age.norm <- sweep(counts[, index == 1], 2, size.factor[index == 1], "/")
        # Subset the manifold
        celltype.mu.age <- manifold[, index == 1]
        # Run SAVER-D on the cell type and age specific subset
        celltype.saver.age <- SAVER::saver(celltype.counts.age.norm, mu = celltype.mu.age, ncores = ncores)

        # assign(paste0("data.size",j), sum(index))
        if (model == "cCV") {
          assign(paste0("data.disp", j), celltype.saver.age$a)
        } else if (model == "cFF") {
          assign(paste0("data.disp", j), celltype.saver.age$b)
        } else if (model == "cVar") {
          assign(paste0("data.disp", j), 1 / celltype.saver.age$k)
        } else {
          stop("Invalid model specified. Choose either 'cCV', 'cFF', or 'cVar'.")
        }
      } else {
        message("Skipped due to low cell counts")
        assign(paste0("data.disp", j), vector())
      }
    }
    # Put together cell type specific dispersion tables
    celltype.disp <- data.frame("gene" = rownames(counts))
    for (j in 1:max(ages.nums)) {
      if (length(get(paste0("data.disp", j))) > 0) {
        celltype.disp[ages.levels[j]] <- get(paste0("data.disp", j))
      }
    }
    dispersion.list[[cell.types.levels[i]]] <- celltype.disp
  }
  return(dispersion.list)
}
