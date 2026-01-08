# #' Preprocess the single-cell raw data using functions in the \code{Seurat} package
# #'
# #' #' This function provide a simplified-version of Seurat analysis pipeline for single-cell RNA-seq data. It contains the following steps in the pipeline:
# #' \itemize{
# #'    \item Normalize the count data present in a given assay.
# #'    \item Identify the variable features.
# #'    \item Scales and centers features in the dataset.
# #'    \item Run a PCA dimensionality reduction.
# #'    \item Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
# #'    \item Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
# #' }
# #'
# #' @title run seurat preprocessing function
# #'
# #' @description Single-cell preprocess
# #'
# #' @details you can use this function to caculate x+1,then return the value of x+1.
# #'
# #' @param Obejct Obejct is a Seurat object
# #' @param normalization.method Method for normalization.
# #'   \itemize{
# #'   \item LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
# #'   This is then natural-log transformed using log1p.
# #'   \item CLR: Applies a centered log ratio transformation.
# #'   \item RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
# #'   No log-transformation is applied. For counts per million (CPM) set \code{scale.factor = 1e6}.
# #' }
# #' @param scale.factor Sets the scale factor for cell-level normalization.
# #' @param selection.method How to choose top variable features. Choose one of :
# #'   \itemize{
# #'   \item vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
# #'   Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
# #'   Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
# #'   \item mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function)
# #'   for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates
# #'   z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong
# #'   relationship between variability and average expression.
# #'   \item dispersion (disp): selects the genes with the highest dispersion values
# #'   }
# #' @param dims_Neighbors Dimensions of reduction to use as input.
# #' @param dims_UMAP Which dimensions to use as input features for UMAP.
# #' @param verbose Print output.
# #'
# #' @return a Seurat object
# #' @import Seurat
# #' @export

# ### Single-cell preprocess
# run_seurat <- function(Obejct,
#                        normalization.method = "LogNormalize",
#                        scale.factor = 10000,
#                        selection.method = "vst",
#                        dims_Neighbors = 1:40,
#                        dims_UMAP = 1:10,
#                        verbose = TRUE){

#   Obejct <- Seurat::NormalizeData(object = Obejct, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
#   Obejct <- Seurat::FindVariableFeatures(object = Obejct, nfeatures = 3000,selection.method = selection.method, verbose = verbose)
#   Obejct <- Seurat::ScaleData(object = Obejct, verbose = verbose)
#   Obejct <- Seurat::RunPCA(object = Obejct, features = VariableFeatures(Obejct), verbose = verbose)
#   Obejct <- Seurat::FindNeighbors(object = Obejct, dims = dims_Neighbors, verbose = verbose)
#   Obejct <- Seurat::RunUMAP(object = Obejct, dims = dims_UMAP, verbose = verbose)
# }

#' @title The scAB_data Class
#'
#' @description
#' The scAB_data object is created from a single-cell RNA-seq data, bulk RNA-seq data, and phenotype data.
#' Its input consists of a Seurat object, a digital matrix, and a matrix or vector providing phenotypic information. Among them, the Seurat object should contain Shared Nearest Neighbor (SNN) Graph information.
#' The class provides functions for data preprocessing, integrative analysis, and visualization.
#'
#' The key slots used in the scAB_data object are described below.
#'
#' @slot X Correlation matrix between bulk data and single-cell data (individuals should be in rows and cells in columns)
#' @slot S Diagonal matrix for information of individuals phenotype
#' @slot A Shared Nearest Neighbor (SNN) Graph matrix
#' @slot D Degree Matrix of Shared Neighbor Graph
#' @slot L Laplacian Matrix
#' @slot phenotype Matrix of phenotype information
#' @slot method method
#'
#' @exportClass scAB_data
#' @useDynLib scAB
#'
#' @import methods
#'
setClass(
  "scAB_data",
  slots = list(
    X = "matrix",
    S = "matrix",
    L = "matrix",
    D = "matrix",
    A = "matrix",
    phenotype = "matrix",
    method = "character"
  )
)


#' @title scAB_data preprocess
#' @description
#' preprocess the single-cell data, bulk data, and phenotype data.
#'
#' @param Object Seurat object
#' @param bulk_dataset matrix of bulk data
#' @param phenotype Phenotype data, a matrix with
#'    two columns "time" and "state", or a vector
#' @param method method "survival" or "binary"
#' @param verbose Logical, whether to print messages.
#' @param ... For future updates.
#'
#' @return a scAB_data
#'
#' @family scAB
#'
#' @export
#' @seealso [guanrank2()]
#'
create_scAB.v5 <- function(
  Object,
  bulk_dataset,
  phenotype,
  method = c("survival", "binary"),
  verbose = SigBridgeRUtils::getFuncOption("verbose") %||% TRUE,
  ...
) {
  dots <- rlang::list2()
  assay <- dots$assay %||% "RNA"

  graph_name <- paste0(assay, "_snn")
  # cell neighbors
  A <- if (graph_name %in% names(Object@graphs)) {
    if (verbose) {
      cli::cli_alert_info(
        " Using {.val {graph_name}} graph for network."
      )
    }
    SeuratObject::Graphs(
      object = Object,
      slot = graph_name
    )
  } else {
    cli::cli_abort(c(
      "x" = "{.val {graph_name}} graph not found in the given Seurat object. 
      Please check `Object@graphs`."
    ))
  }

  if (verbose) {
    ts_cli$cli_alert_info(
      "Calculating L, D and A"
    )
  }

  cal_res <- ComputeNormalizedLaplacian(A)

  # similarity matrix
  sc_exprs <- Matrix::Matrix(SeuratObject::LayerData(Object, assay = assay))
  common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))

  bulk_mat <- Matrix::Matrix(as.matrix(bulk_dataset[common, ]))

  dataset0 <- Matrix::cbind2(bulk_mat, sc_exprs[common, ]) # Dataset before quantile normalization.

  if (verbose) {
    ts_cli$cli_alert_info(
      "Normalizing quantiles of data"
    )
  }

  dataset1 <- SigBridgeRUtils::normalize.quantiles(as.matrix(dataset0)) # Dataset after quantile normalization.

  rownames(dataset1) <- common
  colnames(dataset1) <- colnames(dataset0)

  ncol_bulk <- ncol(bulk_dataset)

  Expression_bulk <- dataset1[, seq_len(ncol_bulk)]
  Expression_cell <- dataset1[, (ncol_bulk + 1):ncol(dataset1)]

  rm(dataset0, dataset1, bulk_mat, sc_exprs, common)
  gc(verbose = FALSE)
  if (verbose) {
    ts_cli$cli_alert_info(
      "Calculating correlation"
    )
  }

  X <- if (rlang::is_installed("WGCNA")) {
    WGCNA::cor(Expression_bulk, Expression_cell)
  } else {
    stats::cor(Expression_bulk, Expression_cell)
  }
  X <- X / Matrix::norm(X, "F")

  if (verbose) {
    ts_cli$cli_alert_info(
      "Guanranking"
    )
  }

  # phenotype ranking
  if (method == "survival") {
    ss <- guanrank2(phenotype[, c("time", "status")])
    S <- Matrix::diag(1 - ss[rownames(phenotype), 3]) # 3 is the rank column
  } else {
    S <- Matrix::diag(1 - phenotype)
  }

  # return
  obj <- list(
    X = Matrix::Matrix(X, sparse = TRUE),
    S = Matrix::Matrix(S, sparse = TRUE),
    L = Matrix::Matrix(cal_res$L, sparse = TRUE),
    D = Matrix::Matrix(cal_res$Dhat, sparse = TRUE),
    A = Matrix::Matrix(cal_res$Ahat, sparse = TRUE),
    phenotype = phenotype,
    method = method
  )
  class(obj) <- "scAB_data"

  obj
}

#' @title Compute Normalized Graph Laplacian and Related Matrices
#' @description
#' This function computes the symmetric normalized graph Laplacian (\eqn{L = I - D^{-1/2} A D^{-1/2}}),
#' along with the normalized adjacency matrix (\eqn{A_{\text{hat}} = D^{-1/2} A D^{-1/2}}) and the
#' normalized degree matrix (\eqn{D_{\text{hat}} = D^{-1/2} D D^{-1/2}}), for a given sparse adjacency matrix.
#'
#' The implementation is optimized for large-scale single-cell SNN graphs (e.g., Seurat's RNA_snn).
#' It avoids explicit sparse matrix multiplications (e.g., \code{D12 \%*\% A \%*\% D12}), which are
#' computationally expensive and may cause intermediate fill-in or memory explosion.
#' Instead, it directly manipulates the non-zero entries of \code{A} via its internal \code{dgCMatrix}
#' structure, achieving \eqn{O(\text{nnz}(A))} time complexity.
#'
#' @param A A square sparse adjacency matrix of class \code{"dgCMatrix"} (recommended) or coercible to it.
#'   Typically represents a k-nearest neighbor (KNN) or shared nearest neighbor (SNN) graph.
#'   Self-loops (diagonal entries) are automatically zeroed if present.
#'   Non-zero entries are binarized to 1 (i.e., unweighted graph), consistent with common GNN/SNN usage.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{L}}{The symmetric normalized graph Laplacian: \eqn{L = D_{\text{hat}} - A_{\text{hat}}}.
#'       Class \code{"dgCMatrix"}. For isolated nodes (degree = 0), the diagonal remains 0 (not 1),
#'       ensuring mathematical correctness: \eqn{L_{ii} = 0} when \eqn{\deg(i) = 0}.}
#'     \item{\code{Dhat}}{The normalized degree matrix: \eqn{D_{\text{hat}} = \text{diag}(d_1, \dots, d_n)},
#'       where \eqn{d_i = 1} if \eqn{\deg(i) > 0}, else \eqn{0}. Class \code{"ddiMatrix"} (diagonal).}
#'     \item{\code{Ahat}}{The symmetric normalized adjacency matrix:
#'       \eqn{(A_{\text{hat}})_{ij} = A_{ij} / \sqrt{\deg(i) \cdot \deg(j)}} for \eqn{A_{ij} \neq 0},
#'       and 0 otherwise. Class \code{"dgCMatrix"}.}
#'   }
#'
#' @details
#' \strong{Key Optimizations and Numerical Safeguards:}
#' \enumerate{
#'   \item \strong{No Matrix Multiplications:} Instead of computing
#'     \code{D12 \%*\% A \%*\% D12} (which involves two sparse-sparse products),
#'     we directly scale each non-zero entry \eqn{A_{ij}} by \eqn{1/\sqrt{\deg(i) \cdot \deg(j)}}.
#'     This leverages the internal structure of \code{dgCMatrix}:
#'     \itemize{
#'       \item \code{A@i}: 0-based row indices of non-zeros.
#'       \item \code{A@p}: column pointers; \code{diff(A@p)} gives non-zero counts per column.
#'       \item \code{j_idx <- rep.int(seq_len(ncol(A)), diff(A@p))} efficiently constructs 1-based column indices
#'         without copying the full index matrix (critical for scalability).
#'     }
#'     Time complexity is strictly \eqn{O(\text{nnz}(A))}, with minimal memory overhead.
#'
#'   \item \strong{Correct Handling of Isolated Nodes:}
#'     For a node \eqn{i} with \eqn{\deg(i) = 0}:
#'     \itemize{
#'       \item Its row/column in \code{A} is all zeros (after zeroing diagonal).
#'       \item \eqn{(D_{\text{hat}})_{ii} = (D^{-1/2} D D^{-1/2})_{ii} = 0}, since \eqn{D_{ii} = 0}.
#'       \item Thus, \eqn{L_{ii} = 0 - 0 = 0}, not 1.
#'     }
#'     The implementation sets \code{dhat_diag[i] <- 0} exactly when \code{degrees[i] == 0}.
#'
#'   \item \strong{Numerical Stability:} A tiny epsilon (\eqn{\epsilon = 2.2204 \times 10^{-256}})
#'     is used to prevent division by zero when computing \eqn{1 / \sqrt{\deg(i)}}.
#'     Note: For any edge \eqn{(i,j)} with \eqn{A_{ij} \neq 0}, we have \eqn{\deg(i) \ge 1} and \eqn{\deg(j) \ge 1}
#'     (since degree is the sum of incident edges), so the scaling factor for existing edges is always finite.
#'     The epsilon only affects the (unused) scaling for zero-degree nodes, ensuring well-defined \code{deg_sqrt_inv}.
#' }
#'
#' \strong{Assumptions:}
#' \itemize{
#'   \item Input matrix is square and represents an undirected graph (symmetry is not enforced but expected).
#'   \item Non-zero entries in \code{A} are positive (binarization sets them to 1).
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate a small SNN graph (500 cells, 5% sparsity, 50 isolated)
#' library(Matrix)
#' set.seed(123)
#' n <- 500
#' A <- rsparsematrix(n, n, density = 0.05)
#' A <- abs(A)  # Ensure positive weights
#' # Force 50 isolated nodes
#' A[(n-49):n, ] <- 0; A[, (n-49):n] <- 0
#' A <- as(A, "dgCMatrix")
#'
#' res <- ComputeNormalizedLaplacian(A)
#' # Verify: L = Dhat - Ahat
#' stopifnot(all.equal(res$L, res$Dhat - res$Ahat, tolerance = 1e-12))
#' }
#'
#' @keywords internal
ComputeNormalizedLaplacian <- function(A) {
  on.exit(gc(verbose = FALSE))
  if (any(Matrix::diag(A) != 0)) {
    Matrix::diag(A) <- 0
  }

  A@x[A@x != 0] <- 1

  degrees <- Matrix::rowSums(A) # numeric vector, length = n

  eps <- 2.2204e-256
  deg_sqrt_inv <- 1 / sqrt(pmax(degrees, eps)) # numeric vector

  # A@i: 0-based row indices (integer vector)
  i_idx <- A@i + 1L # 0-based → 1-based
  # j_idx: 1-based column indices, constructed without full matrix expansion
  # rep.int(seq_len(ncol(A)), diff(A@p)) repeats column number j for nnz_in_col_j times
  j_idx <- rep.int(seq_len(ncol(A)), diff(A@p))

  # Scale each non-zero A[i,j] by 1 / sqrt(deg(i) * deg(j))
  # Critical insight: For any (i,j) with A[i,j] != 0, degrees[i] > 0 and degrees[j] > 0,
  # so deg_sqrt_inv[i] and [j] are finite; no 0 * Inf risk.
  scale_factor <- deg_sqrt_inv[i_idx] * deg_sqrt_inv[j_idx]
  A@x <- A@x * scale_factor # -> Ahat

  # Construct Dhat = D^{-1/2} D D^{-1/2}
  # Diagonal element i: (1/sqrt(deg_i)) * deg_i * (1/sqrt(deg_i)) = 1 if deg_i > 0, else 0
  dhat_diag <- ifelse(degrees > 0, 1, 0)
  Dhat <- Matrix::Diagonal(x = dhat_diag) # ddiMatrix

  L <- Dhat - A

  list(L = L, Dhat = Dhat, Ahat = A)
}
