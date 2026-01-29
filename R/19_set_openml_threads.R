#' @title Set OpenMP threads
#'
#' @param n Number of threads
#' @export
set_openmp_threads <- function(n) {
  if (n > future::availableCores()) {
    cli::cli_abort("Requested threads exceed CPU cores")
  }
  Sys.setenv("OMP_NUM_THREADS" = as.character(n))
}
