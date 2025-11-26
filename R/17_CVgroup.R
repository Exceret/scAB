#' @title Create subsets of cross-validation
#'
#' @param k  k-fold cross validation
#' @param datasize  the size of samples
#' @param seed random seed
#'
#' @return a list with subsets of cross-validation
#' @export
#' @family scAB
#' @examples
#' \dontrun{
#' CVgroup(10,10)
#' }
#'
#'
CVgroup2 <- function(
    k,
    datasize,
    seed = SigBridgeRUtils::getFuncOption("seed")
) {
    set.seed(seed)
    folds <- sample(rep_len(seq_len(k), datasize))
    split(seq_len(datasize), folds)
}
