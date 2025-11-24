#' @title Selection of parameter alpha and alpha_2
#'
#' @description
#' Performs cross-validation to select optimal regularization parameters
#' (alpha_1 and alpha_2) for matrix factorization models. This function
#' evaluates a grid of parameter combinations and selects the ones that
#' maximize the cross-validation performance metric.
#'
#' @param Object A structured object containing data and fixed matrices for
#'               the model. Expected to have components:
#'               \itemize{
#'                 \item \code{X}: The main data matrix
#'                 \item \code{A}: Fixed matrix A
#'                 \item \code{L}: Fixed matrix L
#'                 \item \code{D}: Fixed matrix D
#'               }
#' @param method Type of phenotype
#' @param K The rank for matrix factorization.
#' @param cross_k Number of folds for k-fold cross-validation. Defaults to 5.
#' @param para_1_list Numeric vector of candidate values for the first
#'                    regularization parameter (alpha_1). Defaults to
#'                    \code{c(0.01, 0.005, 0.001)}.
#' @param para_2_list Numeric vector of candidate values for the second
#'                    regularization parameter (alpha_2). Defaults to
#'                    \code{c(0.01, 0.005, 0.001)}.
#' @param seed Random seed for reproducible cross-validation splitting.
#'             Defaults to 0.
#' @param parallel Logical indicating whether to use parallel processing for
#'                 parameter evaluation. Defaults to \code{FALSE}.
#' @param workers Number of parallel workers to use if \code{parallel = TRUE}.
#'                If \code{NULL}, uses available cores minus one. Defaults to
#'                \code{NULL}.
#' @param verbose Logical indicating whether to print progress messages.
#'                Defaults to \code{TRUE}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{para}: A list with selected parameters:
#'       \itemize{
#'         \item \code{alpha_1}: The optimal first regularization parameter
#'         \item \code{alpha_2}: The optimal second regularization parameter
#'         \item \code{result_cv}: Matrix of cross-validation results for all parameter combinations
#'       }
#'   }
#'
#' @note
#' This function supports both parallel and sequential evaluation to balance
#' computational efficiency and resource usage.
#'
#'
#' @export
#' @family scAB
#' @family scAB_optimal_param
#'
#'
select_alpha.optimized <- function(
    Object,
    method = c("binary", "survival"),
    K,
    cross_k = 5,
    para_1_list = c(0.01, 0.005, 0.001),
    para_2_list = c(0.01, 0.005, 0.001),
    seed = SigBridgeRUtils::getFuncOption("seed"),
    parallel = SigBridgeRUtils::getFuncOption("parallel"),
    workers = SigBridgeRUtils::getFuncOption("workers"),
    verbose = SigBridgeRUtils::getFuncOption("verbose")
) {
    if (verbose) {
        ts_cli$cli_alert_info(
            "Selecting optimal {.arg alpha} and {.arg alpha_2}, this would take a while"
        )
    }

    train_phenotype <- PreparePheno(Object)
    train_data <- Object$X

    cvlist <- CVgroup2(k = cross_k, datasize = nrow(train_data), seed = seed)

    train_data_norm <- train_data / Matrix::norm(train_data, "F")
    fixed_matrices <- list(A = Object$A, L = Object$L, D = Object$D)

    param_grid <- expand.grid(
        para_1 = para_1_list,
        para_2 = para_2_list,
        stringsAsFactors = FALSE
    )

    cv_results <- if (parallel) {
        ParallelEvaluate(
            method = method,
            param_grid = param_grid,
            train_data = train_data_norm,
            train_phenotype = train_phenotype,
            cvlist = cvlist,
            fixed_matrices = fixed_matrices,
            K = K,
            cross_k = cross_k,
            workers = workers,
            verbose = verbose
        )
    } else {
        SequentialEvaluate(
            method = method,
            param_grid = param_grid,
            train_data = train_data_norm,
            train_phenotype = train_phenotype,
            cvlist = cvlist,
            fixed_matrices = fixed_matrices,
            K = K,
            cross_k = cross_k,
            verbose = verbose
        )
    }

    result_cv <- matrix(
        cv_results,
        nrow = length(para_1_list),
        ncol = length(para_2_list)
    )

    best_idx <- which(result_cv == max(result_cv), arr.ind = TRUE)[1, ]

    alpha_1 <- para_1_list[best_idx[1]]
    alpha_2 <- para_2_list[best_idx[2]]

    if (verbose) {
        ts_cli$cli_alert_success(
            "Best {.arg alpha} and {.arg alpha_2} are {.val {alpha_1}} and {.val {alpha_2}}"
        )
    }

    list(
        para = list(
            alpha_1 = alpha_1,
            alpha_2 = alpha_2,
            result_cv = result_cv
        )
    )
}

#' @title Prepare Phenotype Data
#' @description
#' Prepares phenotype data in the appropriate format for different analysis methods.
#' For survival analysis, returns the phenotype as-is. For other methods, converts
#' to a data frame with status and time columns.
#'
#' @param Object A list containing analysis method and phenotype data
#'
#' @return For survival method: returns Object$phenotype as-is.
#' For other methods: returns a data frame with status and time columns.
#'
#' @family scAB_optimal_param
#' @export
PreparePheno <- function(Object) {
    if (Object$method == "survival") {
        return(Object$phenotype)
    }

    data.frame(
        status = as.integer(Object$phenotype),
        time = ifelse(Object$phenotype, 1, 100),
        row.names = rownames(Object$X)
    )
}

#' @title Sequential Parameter Evaluation
#' @description
#' Evaluates parameter combinations sequentially using cross-validation.
#' For each parameter set in the grid, computes average cross-validation score.
#'
#' @param method Analysis method ("binary" or "survival")
#' @param param_grid Data frame of parameter combinations to evaluate
#' @param train_data Training data matrix
#' @param train_phenotype Phenotype data for training
#' @param cvlist List of cross-validation indices
#' @param fixed_matrices Precomputed matrices (A, L, D)
#' @param K Number of components
#' @param cross_k Number of cross-validation folds
#' @param verbose Whether to show progress
#'
#' @return Numeric vector of average CV scores for each parameter combination
#'
#' @family scAB_optimal_param
#' @export
SequentialEvaluate <- function(
    method = c("binary", "survival"),
    param_grid,
    train_data,
    train_phenotype,
    cvlist,
    fixed_matrices,
    K,
    cross_k,
    verbose = TRUE
) {
    cv_cache <- vector("list", cross_k)

    if (verbose) {
        cli::cli_progress_bar(name = 'Preparing CV cache', total = cross_k)
    }

    for (cv_idx in seq_len(cross_k)) {
        test_idx <- cvlist[[cv_idx]]
        train_pheno <- train_phenotype[-test_idx, , drop = FALSE]

        ss <- guanrank2(train_pheno[, c("time", "status")])
        S_matrix <- diag(1 - ss[rownames(train_pheno), 3])

        cv_cache[[cv_idx]] <- list(
            train_subset = train_data[-test_idx, , drop = FALSE],
            test_subset = train_data[test_idx, , drop = FALSE],
            train_pheno = train_pheno,
            test_pheno = train_phenotype[test_idx, , drop = FALSE],
            S_matrix = S_matrix,
            Object_template = list(
                phenotype = train_pheno,
                A = fixed_matrices$A,
                L = fixed_matrices$L,
                D = fixed_matrices$D,
                method = method
            )
        )
        if (verbose) {
            cli::cli_progress_update()
        }
    }
    if (verbose) {
        cli::cli_progress_done()
    }
    n_params <- nrow(param_grid)

    cv_scores_matrix <- matrix(0, nrow = n_params, ncol = cross_k)

    for (param_idx in seq_len(n_params)) {
        para_1 <- param_grid$para_1[param_idx]
        para_2 <- param_grid$para_2[param_idx]

        if (verbose) {
            cli::cli_progress_bar(
                name = sprintf(
                    'Alpha = %g, alpha_2 = %g',
                    para_1,
                    para_2
                ),
                total = cross_k
            )
        }
        for (cv_idx in seq_len(cross_k)) {
            cache <- cv_cache[[cv_idx]]

            Object_cv <- c(
                list(
                    X = cache$train_subset,
                    S = cache$S_matrix
                ),
                cache$Object_template
            )
            class(Object_cv) <- "scAB_data"

            s_res <- scAB.optimized(
                Object = Object_cv,
                K = K,
                alpha = para_1,
                alpha_2 = para_2,
                maxiter = 2000
            )

            ginvH <- SigBridgeRUtils::ginv2(s_res$H)
            new_W <- cache$test_subset %*% ginvH

            W_matrix <- as.matrix(s_res$W)

            n_train <- nrow(W_matrix)
            n_features <- ncol(W_matrix)

            clin_data <- cbind(
                cache$train_pheno$time,
                cache$train_pheno$status,
                W_matrix
            )

            col_names <- c("time", "status", paste0("V", seq_len(n_features)))
            clin_km <- as.data.frame(clin_data)
            colnames(clin_km) <- col_names

            res.cox <- survival::coxph(
                survival::Surv(time, status) ~ .,
                data = clin_km,
                x = FALSE,
                y = FALSE
            )

            new_W_df <- as.data.frame(new_W)
            colnames(new_W_df) <- paste0("V", seq_len(n_features))

            pre_test <- stats::predict(res.cox, newdata = new_W_df, type = "lp")

            cv_scores_matrix[param_idx, cv_idx] <- survival::concordance(
                survival::Surv(
                    cache$test_pheno$time,
                    cache$test_pheno$status
                ) ~ pre_test
            )$concordance

            if (verbose) {
                cli::cli_progress_update()
            }
        }

        if (verbose) {
            cli::cli_progress_done()
        }
    }

    rowMeans(cv_scores_matrix)
}


#' @title Parallel Parameter Evaluation
#' @description
#' Evaluates parameter combinations in parallel using cross-validation.
#' Uses future framework for parallel processing across multiple workers.
#'
#' @param method Analysis method ("binary" or "survival")
#' @param param_grid Data frame of parameter combinations to evaluate
#' @param train_data Training data matrix
#' @param train_phenotype Phenotype data for training
#' @param cvlist List of cross-validation indices
#' @param fixed_matrices Precomputed matrices (A, L, D)
#' @param K Number of components
#' @param cross_k Number of cross-validation folds
#' @param workers Number of parallel workers
#' @param verbose Whether to show progress
#' @param parallel_type Type of parallel backend to use
#' @param seed Random seed for reproducibility
#'
#' @return Numeric vector of average CV scores for each parameter combination
#'
#' @family scAB_optimal_param
#' @export
ParallelEvaluate <- function(
    method = c("binary", "survival"),
    param_grid,
    train_data,
    train_phenotype,
    cvlist,
    fixed_matrices,
    K,
    cross_k,
    workers = SigBridgeRUtils::getFuncOption("workers"),
    verbose = SigBridgeRUtils::getFuncOption("verbose"),
    parallel_type = SigBridgeRUtils::getFuncOption("parallel.type"),
    seed = SigBridgeRUtils::getFuncOption("seed")
) {
    SigBridgeRUtils::plan(parallel_type, workers = workers)
    on.exit(SigBridgeRUtils::plan("sequential"))

    if (verbose) {
        ts_cli$cli_alert_info(sprintf(
            "Using parallel processing with %d workers",
            workers
        ))
    }

    if (verbose) {
        cli::cli_progress_bar(name = 'Preparing CV cache', total = cross_k)
    }

    cv_cache <- vector("list", cross_k)

    for (cv_idx in seq_len(cross_k)) {
        test_idx <- cvlist[[cv_idx]]
        train_pheno <- train_phenotype[-test_idx, , drop = FALSE]

        ss <- guanrank2(train_pheno[, c("time", "status")])
        S_matrix <- Matrix::diag(1 - ss[rownames(train_pheno), 3])

        cv_cache[[cv_idx]] <- list(
            train_subset = train_data[-test_idx, , drop = FALSE],
            test_subset = train_data[test_idx, , drop = FALSE],
            train_pheno = train_pheno,
            test_pheno = train_phenotype[test_idx, , drop = FALSE],
            S_matrix = S_matrix,
            Object_template = list(
                phenotype = train_pheno,
                A = fixed_matrices$A,
                L = fixed_matrices$L,
                D = fixed_matrices$D,
                method = method
            )
        )
        if (verbose) {
            cli::cli_progress_update()
        }
    }
    if (verbose) {
        cli::cli_progress_done()
    }
    # === 并行评估函数（内联版本）===
    ScoringAll <- function(param_idx) {
        para_1 <- param_grid$para_1[param_idx]
        para_2 <- param_grid$para_2[param_idx]

        cv_scores <- vapply(
            seq_len(cross_k),
            function(cv_idx) {
                cache <- cv_cache[[cv_idx]]

                Object_cv <- cache$Object_template
                Object_cv$X <- cache$train_subset
                Object_cv$S <- cache$S_matrix
                class(Object_cv) <- "scAB_data"

                s_res <- scAB.optimized(
                    Object = Object_cv,
                    K = K,
                    alpha = para_1,
                    alpha_2 = para_2,
                    maxiter = 2000L
                )

                ginvH <- SigBridgeRUtils::ginv2(s_res$H)
                new_W <- cache$test_subset %*% ginvH
                W_matrix <- as.matrix(s_res$W)

                n_features <- ncol(W_matrix)
                clin_data <- cbind(
                    cache$train_pheno$time,
                    cache$train_pheno$status,
                    W_matrix
                )

                col_names <- c(
                    "time",
                    "status",
                    paste0("V", seq_len(n_features))
                )
                clin_km <- as.data.frame(clin_data)
                colnames(clin_km) <- col_names

                res.cox <- survival::coxph(
                    survival::Surv(time, status) ~ .,
                    data = clin_km,
                    x = FALSE,
                    y = FALSE
                )

                new_W_df <- as.data.frame(new_W)
                colnames(new_W_df) <- paste0("V", seq_len(n_features))
                pre_test <- stats::predict(
                    res.cox,
                    newdata = new_W_df,
                    type = "lp"
                )

                survival::concordance(
                    survival::Surv(
                        cache$test_pheno$time,
                        cache$test_pheno$status
                    ) ~ pre_test
                )$concordance
            },
            numeric(1)
        )

        mean(cv_scores)
    }

    res <- SigBridgeRUtils::future_map_dbl(
        seq_len(nrow(param_grid)),
        ScoringAll,
        .progress = verbose,
        .options = SigBridgeRUtils::furrr_options(
            seed = seed,
            packages = c(
                "survival",
                "SigBridgeRUtils",
                "Matrix",
                "Rcpp",
                "scAB"
            ),
            globals = list(
                cv_cache = cv_cache,
                K = K,
                method = method,
                cvlist = cvlist,
                cross_k = cross_k
            )
        )
    )

    res
}
