#' @title identification cells above the threshold
#'
#' @param H cell matrix
#' @param tred threshold
#' @param do.dip Logical, whether to do.dip
#'
#' @return a list with cells above the threshold
#'
#' @export
#'
findModule <- function(H, tred = 2, do.dip = FALSE) {
  K <- nrow(H)
  meanH <- SigBridgeRUtils::rowMeans3(H)
  sdH <- SigBridgeRUtils::rowSds3(H)

  if (!do.dip) {
    threshold_matrix <- H - meanH > tred * sdH
    module <- apply(threshold_matrix, 1, which, simplify = FALSE)
    return(module)
  }

  module <- vector("list", K)
  for (i in seq_len(K)) {
    x <- H[i, ]
    if (diptest::dip.test(x)$p.value < 0.05) {
      modes <- multimode::locmodes(x, mod0 = 2)
      module[[i]] <- which(x > modes$locations[2])
    } else {
      module[[i]] <- which(x - meanH[i] > tred * sdH[i])
    }
  }

  module
}
