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
    K = dim(H)[1]
    I = length(H)
    module = list()
    meanH = rowMeans(H)
    sdH = apply(H, 1, stats::sd)
    for (i in 1:K) {
        x <- H[i, ]
        if (do.dip && diptest::dip.test(x)$p.value < 0.05) {
            modes <- multimode::locmodes(x, mod0 = 2)
            module = c(module, list(which(x > modes$locations[2])))
        } else {
            module = c(module, list(which(H[i, ] - meanH[i] > tred * sdH[i])))
        }
    }
    return(module)
}
