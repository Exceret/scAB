#' @title Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization.
#' @description
#' Non-negative Matrix Factorization with phenotype and cell-cell similarity regularization,
#' for identifing phenotype-associated cell states at different resolutions.
#'
#' @param Object  a scAB_data object
#' @param K  the rank of matrix factorization
#' @param maxiter the maximum number of iterations
#' @param alpha Coefficient of phenotype regularization
#' @param alpha_2 Coefficient of cell-cell similarity regularization,
#' @param convergence_threshold The threshold of convergence when computing the loss
#'
#' @return a list with the submatrix and loss value
#' @export
#' @family scAB
#'
#'
scAB.optimized <- function(
  Object,
  K,
  alpha = 0.005,
  alpha_2 = 0.005,
  maxiter = 2000L,
  convergence_threshold = 1e-5
) {
  seed <- ifelse(Object$method == "survival", 7L, 5L)
  if (Object$method != "") {
    set.seed(seed)
  }
  # Cpp func
  result <- scAB_inner(
    X = Matrix::Matrix(Object$X, sparse = TRUE),
    A = Matrix::Matrix(Object$A, sparse = TRUE),
    D = Matrix::Matrix(Object$D, sparse = TRUE),
    L = Matrix::Matrix(Object$L, sparse = TRUE),
    S = Matrix::Matrix(Object$S, sparse = TRUE),
    K = K,
    alpha = alpha,
    alpha_2 = alpha_2,
    maxiter = maxiter,
    convergence_threshold = convergence_threshold
  )

  result$method <- Object$method

  result
}
