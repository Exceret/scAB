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
#' @param phenotype Phenotype data, a matrix with two columns "time" and "state", or a vector
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
    verbose = SigBridgeRUtils::getFuncOption("verbose"),
    ...
) {
    # cell neighbors
    if ("RNA_snn" %in% names(Object@graphs)) {
        A <- Matrix::Matrix(SeuratObject::Graphs(
            object = Object,
            slot = "RNA_snn"
        ))
        if (verbose) {
            cli::cli_alert_info(
                " Using {.val RNA_snn} graph for network."
            )
        }
    } else if ("integrated_snn" %in% names(Object@graphs)) {
        A <- Matrix::Matrix(SeuratObject::Graphs(
            object = Object,
            slot = "integrated_snn"
        ))

        if (verbose) {
            cli::cli_alert_info(
                "Using {.val integrated_snn} graph for network."
            )
        }
    } else {
        cli::cli_abort(c(
            "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check `Object@graphs`."
        ))
    }
    Matrix::diag(A) <- 0
    A@x[which(A@x != 0)] <- 1
    degrees <- Matrix::rowSums(A)
    D <- Matrix::diag(degrees)
    eps <- 2.2204e-256 # In the original implementation of scAB, a custom-defined `eps` is used instead of `.Machine$double.eps`.
    D12 <- Matrix::diag(1 / sqrt(pmax(degrees, eps)))

    L <- D12 %*% (D - A) %*% D12 # Normalized Graph Laplacian
    Dhat <- D12 %*% (D) %*% D12
    Ahat <- D12 %*% (A) %*% D12

    # similarity matrix
    sc_exprs <- Matrix::Matrix(SeuratObject::LayerData(Object))
    common <- intersect(rownames(bulk_dataset), rownames(sc_exprs))
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ]) # Dataset before quantile normalization.
    dataset1 <- SigBridgeRUtils::normalize.quantiles(as.matrix(dataset0)) # Dataset after quantile normalization.
    dataset1 <- Matrix::Matrix(dataset1)
    rownames(dataset1) <- common
    colnames(dataset1) <- colnames(dataset0)

    ncol_bulk <- ncol(bulk_dataset)
    Expression_bulk <- as.matrix(dataset1[, seq_len(ncol_bulk)])
    Expression_cell <- as.matrix(dataset1[, (ncol_bulk + 1):ncol(dataset1)])

    X <- stats::cor(Expression_bulk, Expression_cell)
    X <- X / Matrix::norm(X, "F")

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
        L = Matrix::Matrix(L, sparse = TRUE),
        D = Matrix::Matrix(Dhat, sparse = TRUE),
        A = Matrix::Matrix(Ahat, sparse = TRUE),
        phenotype = phenotype,
        method = method
    )
    class(obj) <- "scAB_data"

    obj
}
