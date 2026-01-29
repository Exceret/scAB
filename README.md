<!-- README.md is generated from README.Rmd. Please edit that file -->

# scAB detects multiresolution cell states with clinical significance by integrating single-cell genomics and bulk sequencing data

Forked from [scAB](https://github.com/jinworks/scAB), see
[SigBridgeR](https://github.com/WangLabCSU/SigBridgeR) for the usage of
scAB. Independently usage of scAB is also viable.

## News:

-   2025-11-21:
    -   Tweaks for performance when selecting optimal alpha
-   2025-11-18:
    -   tweaks for performance
    -   deleted tutorial data and preprocessing code
    -   writed some documentation
    -   added some functions

## Installation

scAB package can be easily installed from Github using `pak`:

    pak::pkg_install("Exceret/scAB")

Install with `ARMA_64BIT_WORD` enabled:

    repo_url <- "https://github.com/Exceret/scAB.git"
    pkg_path <- file.path(tempdir(), "scAB")

    system2("git", c("clone", repo_url, pkg_path))

    makevars_file <- if (.Platform$OS.type == "windows") {
      file.path(pkg_path, "src", "Makevars.win")
    } else {
      file.path(pkg_path, "src", "Makevars")
    }
    if (file.exists(makevars_file)) {
      lines <- readLines(makevars_file)
      lines <- gsub(
        "^(PKG_CXXFLAGS.*)$",
        "\\1 -DARMA_64BIT_WORD",
        lines
      )
      writeLines(lines, makevars_file)
    }

    devtools::install(pkg_path, build_vignettes = FALSE)

see `DESCRIPTION` for more details
