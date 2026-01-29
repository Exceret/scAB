# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
}

.onLoad <- function(libname, pkgname) {
  Sys.setenv("OMP_NUM_THREADS" = "1")
  invisible()
}

#' Add timestamp to cli functions
#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()
