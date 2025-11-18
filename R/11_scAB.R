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
    X <- Object$X
    A <- Object$A
    L <- Object$L
    D <- Object$D
    S <- Object$S
    eps <- 2.2204e-256
    nr <- nrow(X)
    nc <- ncol(X)
    W <- Matrix::Matrix(stats::runif(nr * K), nrow = nr, ncol = K)
    H <- Matrix::Matrix(stats::runif(K * nc), nrow = K, ncol = nc)
    SS <- S %*% S

    loss_func <- function(X, W, H, S, L, alpha, alpha_2) {
        # loss <- norm(X - W %*% H, "F")^2 +
        #     alpha * (norm(S %*% W, "F")^2) +
        #     alpha_2 * sum(diag(H %*% L %*% t(H)))
        loss1 <- Matrix::norm(X - W %*% H, "F")^2
        loss2 <- alpha * (Matrix::norm(S %*% W, "F")^2)
        loss3 <- alpha_2 * sum(H * (H %*% L)) # diag(H %*% L %*% t(H)) is equivalent to rowSums(H * (H %*% L))

        return(loss1 + loss2 + loss3)
    }
    # tD <- t(D)
    # tA <- t(A)
    tX <- Matrix::t(X)

    for (iter in seq_len(maxiter)) {
        # W_old <- W # no use
        H <- H *
            (crossprod(W, X) + alpha_2 * tcrossprod(H, A)) /
            (crossprod(W, W) %*% H + alpha_2 * tcrossprod(H, D) + eps)
        Pena <- SS %*% W
        W <- W *
            tcrossprod(X, H) /
            (W %*% tcrossprod(H, H) + alpha * Pena + eps)

        # tW <- t(W) # no use

        if (iter != 1) {
            eucl_dist <- loss_func(X, W, H, S, L, alpha, alpha_2)
            d_eucl <- abs(eucl_dist - old_eucl)
            if (d_eucl < convergence_threshold) {
                break
            }
            old_eucl <- eucl_dist
        } else {
            old_eucl <- loss_func(X, W, H, S, L, alpha, alpha_2)
        }
    }
    list(
        W = W,
        H = H,
        iter = iter,
        loss = old_eucl,
        method = Object$method
    )
}
