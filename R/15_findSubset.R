#' @title Subsets identification
#'
#' @param Object a Seurat object
#' @param scAB_Object a scAB_data object from scAB()
#' @param tred threshold
#'
#' @return a Seurat object
#' @export
#' @family scAB
#'
#'
findSubset.optimized <- function(Object, scAB_Object, tred = 2L) {
    do.dip <- ifelse(scAB_Object$method == "binary", 1L, 0L)
    H <- as.matrix(scAB_Object$H)
    module <- findModule(H, tred = tred, do.dip = do.dip)
    scAB_index <- unique(unlist(module))
    # Add scAB column in metadata
    n_cells <- ncol(Object)
    scAB_select <- rep("Other", n_cells)
    scAB_select[scAB_index] <- "Positive"
    Object$scAB <- scAB_select

    for (i in seq_along(module)) {
        M <- rep("Other", n_cells)
        M[as.numeric(module[[i]])] <- "Positive"
        Object <- Seurat::AddMetaData(
            Object,
            metadata = M,
            col.name = paste0("scAB_Subset", i)
        )
        Object <- Seurat::AddMetaData(
            Object,
            metadata = H[i, ],
            col.name = paste0("Subset", i, "_loading")
        )
    }

    Object
}
