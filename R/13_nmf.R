#' @title Non-negative Matrix Factorization from scAB
#'
#' @description
#' An R implementation of classical non-negative matrix factorization (NMF).
#'
#' @param Object a scAB_data object or a data matrix
#' @param K The rank of factorization, i.e., the number of components/latent
#'          features to extract. Must be a positive integer smaller than
#'          both dimensions of X.
#' @param maxiter The maximum number of iterations for the optimization algorithm.
#'                Defaults to 2000.
#' @param tol The convergence tolerance for the loss function. The algorithm
#'            stops when the absolute change in Euclidean distance between
#'            consecutive iterations is less than this value. Defaults to 1e-5.
#'
#' @return A list containing the factorization results:
#' \itemize{
#'   \item \code{W} - The basis matrix (features x K), representing the
#'                    learned features or components
#'   \item \code{H} - The coefficient matrix (K x samples), representing the
#'                    weights or activations of components for each sample
#'   \item \code{iter} - The number of iterations actually performed
#'   \item \code{loss} - The final value of the objective function (squared
#'                       Euclidean distance)
#' }
#'
#' @seealso
#' This function is from scAB package, and it is not recommended to use it because the computational efficiency of the R language is not very high.
#' For more advanced NMF implementations, see:
#' \code{\link[NMF]{nmf}} from the NMF package, which is written in C++ and is much faster than this function.
#'
#' @export
#' @family scAB
#'
NMF.optimized <- function(Object, K, maxiter = 2000L, tol = 1e-5) {
    X = if (inherits(Object, "scAB_data")) Object$X else Object
    # Cpp Func
    NMF_optimized(X = as.matrix(X), K = K, maxiter = maxiter, tol = tol)
}
