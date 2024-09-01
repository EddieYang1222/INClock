#' Compute the manifold of gene expressions
#'
#' @param counts An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells. It can be sparse.
#' @param method Method for estimating the manifold. Can be either "SAVER",
#' which uses SAVER to recover the true expression, or 'neighbor',
#' which uses nearest neighbors to approximate the true expression
#' @param preprocess Whether the expression count matrix should be preprocessed,
#' where all genes with zero counts will removed. Default is TRUE.
#' @param ncores Number of cores to use. Default is 1. (SAVER only)
#' @param nfeatures Number of features to use for FindVariableFeatures. Default is 3000. (Neighbor-based only)
#' @param dims Number of dimensions for nearest neighbors. Default is 20.(Neighbor-based only)
#' @param neighbors Number of neighbors to output. Default is 20. (Neighbor-based only)
#'
#' @return A matrix of estimated gene expressions.
#' @importFrom Matrix sparseMatrix rowSums rowMeans
#' @importFrom pbapply pbapply
#' @importFrom methods as
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures ScaleData RunPCA FindNeighbors GetAssayData VariableFeatures
#' @import TabulaMurisSenisData
#' @export
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
#' # tms_marrow_manifold_saver <- compute_manifold(tms_marrow_counts, method = "SAVER")
#' tms_marrow_manifold_neighbor <- compute_manifold(tms_marrow_counts, method = "neighbor")
compute_manifold <- function(counts, method = "SAVER", preprocess = TRUE,
                             ncores = 1, nfeatures = 3000, dims = 20, neighbors = 20) {
  message("Starting manifold computation...\n")

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Count matrix is missing row names or column names.")
  }

  if (preprocess) {
    message("The initial matrix size is ", nrow(counts), " genes and ", ncol(counts), " cells.")
    counts <- counts[rowSums(counts) != 0, ]
    message(
      "After removing genes with zero counts, the new matrix size is ",
      nrow(counts), " genes and ", ncol(counts), " cells.\n"
    )
  }

  if (method == "SAVER") {
    # Estimate manifold with SAVER
    message("Estimating the manifold using the SAVER method...\n")
    start.time <- Sys.time()
    manifold <- SAVER::saver(counts, ncores = ncores)
    end.time <- Sys.time()
    time.taken <- round(difftime(end.time, start.time, units = "secs"), 2)
    message("Finished computing manifold using SAVER. Time taken: ", time.taken, " seconds.\n")
    return(manifold)
  } else if (method == "neighbor") {
    # Find neighbors
    message("Estimating the manifold using the neighbor-based method...\n")
    manifold_obj <- CreateSeuratObject(counts = counts)
    manifold_obj <- NormalizeData(manifold_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    manifold_obj <- FindVariableFeatures(manifold_obj, selection.method = "vst", nfeatures = nfeatures)
    manifold_obj <- ScaleData(manifold_obj, features = rownames(manifold_obj))
    manifold_obj <- RunPCA(manifold_obj, features = VariableFeatures(object = manifold_obj))
    manifold_obj <- FindNeighbors(manifold_obj, dims = 1:dims, k.param = neighbors + 1, return.neighbor = TRUE)

    # Calculate mean count for each cell
    manifold_counts <- GetAssayData(object = manifold_obj, assay = "RNA", layer = "counts")
    print(paste0("Successfully extracted the nearest neighbors. Estimation procedure is starting."))
    manifold <- sparseMatrix(
      i = integer(0), j = integer(0),
      dims = manifold_counts@Dim, dimnames = manifold_counts@Dimnames
    )
    manifold <- as(manifold, "dMatrix")
    manifold_neighbor_idx <- manifold_obj@neighbors$RNA.nn@nn.idx

    start.time <- Sys.time()

    manifold <- pbapply(manifold_neighbor_idx[, -1], 1, function(idx) {
      rowMeans(manifold_counts[, idx])
    })

    dimnames(manifold) <- dimnames(manifold_counts)

    end.time <- Sys.time()
    time.taken <- round(difftime(end.time, start.time, units = "secs"), 2)
    message("Finished computing manifold using the neighbor-based method. Time taken: ", time.taken, " seconds.\n")
    return(manifold)
  } else {
    stop("Invalid method specified. Choose either 'SAVER' or 'neighbor'.")
  }
}
