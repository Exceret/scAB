#' @title Selection of Parameter K for Non-negative Matrix Factorization
#'
#' @description
#' Automatically determines the optimal rank (K) for non-negative matrix
#' factorization using an empirical indicator method. This function evaluates
#' multiple candidate ranks and selects the one that provides the best
#' trade-off between model complexity and reconstruction accuracy.
#'
#' @param Object A scAB_data object containing the data matrix to be factorized.
#' @param K_max The maximum rank value to consider in the search. Must be at
#'              least 2. Defaults to 20.
#' @param repeat_times The number of repeated NMF runs for each candidate rank
#'                     to account for random initialization variability.
#'                     Defaults to 10.
#' @param maxiter The maximum number of iterations for each NMF run.
#'                Defaults to 2000.
#' @param seed Random seed for reproducible results. Defaults to 0.
#' @param verbose Logical indicating whether to print progress messages and
#'                intermediate results. Defaults to FALSE.
#'
#' @return An integer value representing the selected optimal rank K.
#'
#' @note
#' This function is from scAB package,
#'
#' @family scAB
#'
#' @export
#' @seealso [NMF.optimized()]
#'
select_K.optimized <- function(
    Object,
    K_max = 20L,
    repeat_times = 10L,
    maxiter = 2000L,
    seed = SigBridgeRUtils::getFuncOption("seed"),
    verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
    # Cpp func
    select_K_optimized(
        X = Object$X,
        K_max = K_max,
        repeat_times = repeat_times,
        maxiter = maxiter,
        seed = seed,
        verbose = verbose
    )
}
