## name: guanrank.R
## date: 02/20/2018

## This is the R version of guanrank, KM-curve and guanrank-elf

#' @title Kaplan-Meier Estimator
#' @description
#' Calculates Kaplan-Meier survival curve from time-to-event data.
#'
#' @param mat Matrix or data frame with time and status columns. First two columns
#' must be time and status; rownames are used as identifiers.
#'
#' @return Matrix with three columns: time, status, and survival_rate.
#'
#' @examples
#' \dontrun{
#' # Create sample survival data
#' surv_data <- data.frame(
#' time = c(5, 8, 12, 15, 20, 25),
#' status = c(1, 0, 1, 1, 0, 1)
#' )
#' rownames(surv_data) <- paste0("patient", 1:6)
#'
#' # Calculate KM curve
#' km_result <- km_curve(surv_data)
#' print(km_result)
#' }
#'
#' @export
km_curve <- function(mat) {
    # 1. preprocess
    mat = as.matrix(mat)
    mat = mat[, 1:2]
    colnames(mat) = c("time", "status")
    storage.mode(mat) <- "numeric" # convert: numeric matrix
    mat = mat[order(mat[, "time"]), ] # order by time
    mat_curve = cbind(mat, survival_rate = rep(1, nrow(mat)))

    # 2. calculation
    unique_time = unique(mat[, "time"]) # calculate sr at each unique time
    sr = 1 # initial survival rate
    num_current = nrow(mat) # number of current/uncensored/at-risk sample
    for (i in unique_time) {
        ind = mat[, "time"] == i
        censor = mat[ind, "status"]
        sr = sr * (1 - sum(censor) / num_current) # KM estimator
        mat_curve[ind, "survival_rate"] = sr
        num_current = num_current - sum(ind) # update number of samples at risk
    }

    return(mat_curve)
}

# expected future lifetime
# This function calculates the area under KM curve
# Input: A data.frame consists of three columns - time, status, survival_rate; it's also the return of the function km_curve
# Output: A matrix consists of four columns. The last column is the accumulated area under KM curve called future_time.
# The accumulation is calculated in the reverse order so that early cencoring time has more "future time" and vise versa.
exp_future_time <- function(mat, unit_time = FALSE) {
    # 1. preprocess
    mat = as.matrix(mat)
    storage.mode(mat) <- "numeric" # convert: numeric matrix
    if (unit_time) {
        # if unit_time, ignore length of time interval and only consider order
        mat[, "time"] = order(mat[, "time"])
    }
    mat = mat[order(mat[, "time"]), ] # order by time
    mat_future = cbind(mat, future_time = rep(1, nrow(mat)))

    # 2. calculation
    accumulated_time = mat_future[nrow(mat), "future_time"] = 0 # end point is set to 0
    for (i in (nrow(mat) - 1):1) {
        accumulated_time = accumulated_time +
            (mat[i + 1, "time"] - mat[i, "time"]) * mat[i, "survival_rate"] # use i survival rate
        mat_future[i, "future_time"] = accumulated_time
    }

    return(mat_future)
}

# guanrank_elf: Expected LiFetime
# This function calculates guanrank_elf based on the input
# Input: The first two columns are time and status. The rowname is patient id
# Output: A matrix consists of multiple columns.
# Arguments:
# (1) unit_time: if TRUE, the lengths of time intervals are converted into 1
# (2) complete: if TRUE, a complete table is returned; otherwise only three columns are returned: time, status and guanrank
guanrank_elf <- function(mat, unit_time = FALSE, complete = FALSE) {
    # 1. KM_curve & expected future time
    mat = as.matrix(mat)
    mat = mat[, 1:2]
    colnames(mat) = c("time", "status")
    storage.mode(mat) <- "numeric" # convert: numeric matrix
    mat = mat[order(mat[, "time"]), ] # order by time
    mat_curve = km_curve(mat)
    mat_future = exp_future_time(mat_curve, unit_time = unit_time)
    mat_elf = cbind(mat_future, exp_time = rep(1, nrow(mat)))

    # 2. calculation
    l = nrow(mat)
    mat_elf[l, "exp_time"] = mat_elf[l, "time"]
    for (i in (l - 1):1) {
        if (mat_elf[i, "status"]) {
            mat_elf[i, "exp_time"] = mat_elf[i, "time"]
        } else {
            tmp = mat_elf[i, "future_time"] / mat_elf[i, "survival_rate"]
            if (is.na(tmp)) {
                # for end point whose survival rate is 0, estimate it using sr=0.5
                tmp = (mat_elf[i + 1, "time"] - mat_elf[i, "time"]) * 0.5
            }
            mat_elf[i, "exp_time"] = mat_elf[i, "time"] + tmp
        }
    }

    # 3. scale by (x-min)/(max-min)
    val = mat_elf[, "exp_time"]
    val_max = max(val)
    val_min = min(val)
    mat_elf = cbind(
        mat_elf,
        scaled_time = (val - val_min) / (val_max - val_min)
    )

    # 4. reverse rank (max corresponds to 0, not 1)
    mat_elf = cbind(mat_elf, rank = 1 - mat_elf[, "scaled_time"])
    if (!complete) {
        mat_elf = mat_elf[, c("time", "status", "rank")]
    }

    return(mat_elf)
}


#' @title Guanrank for Survival Data
#'
#' @param mat A matrix, data frame, or data table containing survival data.
#'            Must contain at least two columns representing time and status.
#'            The function will use the first two columns as time and status
#'            respectively.
#' @param complete Logical indicating whether to include the survival rate
#'                 column in the output. If \code{FALSE} (default), returns
#'                 only time, status, and rank columns.
#'
#' @return A matrix with the following columns:
#'   \itemize{
#'     \item \code{time}: Original survival times (sorted)
#'     \item \code{status}: Original event indicators (0 = censored, 1 = event)
#'     \item \code{rank}: Computed Guan's rank values (normalized to \[0,1\])
#'     \item \code{survival_rate}: (if \code{complete = TRUE}) Kaplan-Meier survival probabilities from \code{km_curve}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Create sample survival data (matrix format)
#' survival_data <- matrix(c(
#'   10, 1,   # Event at time 10
#'   15, 0,   # Censored at time 15
#'   20, 1,   # Event at time 20
#'   25, 0,   # Censored at time 25
#'   30, 1    # Event at time 30
#' ), ncol = 2, byrow = TRUE)
#'
#' # Compute Guan's rank
#' result <- guanrank(survival_data)
#' print(result)
#'
#' # Get complete output with survival rates
#' complete_result <- guanrank(survival_data, complete = TRUE)
#' print(complete_result)
#'
#' # Using data frame input
#' df_data <- data.frame(
#'   time = c(10, 15, 20, 25, 30),
#'   status = c(1, 0, 1, 0, 1),
#'   other_col = c("A", "B", "C", "D", "E")  # This column will be ignored
#' )
#' df_result <- guanrank(df_data)
#' print(df_result)
#' }
#'
#' @family scAB
#' @export
#'
guanrank2 <- function(mat, complete = FALSE) {
    mat <- as.matrix(mat)
    mat <- mat[, 1:2]
    colnames(mat) = c("time", "status")
    storage.mode(mat) = "numeric"
    mat <- mat[order(mat[, "time"]), ]

    mat_curve <- scAB::km_curve(mat)
    n <- nrow(mat_curve)

    mat_guanrank <- cbind(mat_curve, rank = rep(0, n))
    vect <- mat_guanrank[, "time"]
    vecs <- mat_guanrank[, "status"]
    # vectorized pipeline is not faster than the loop as I expected.
    for (i in seq_len(n)) {
        tA <- mat_guanrank[i, "time"]
        rA <- mat_guanrank[i, "survival_rate"]
        sA <- mat_guanrank[i, "status"]
        if (sA == 1) {
            tBgttA <- mat_guanrank[vect > tA, "survival_rate"]
            tBletA_sBeq0 <- mat_guanrank[
                vect <= tA & vecs == 0,
                "survival_rate"
            ]
            tBeqtA_sBeq1 <- mat_guanrank[
                vect == tA & vecs == 1,
                "survival_rate"
            ]
            mat_guanrank[i, "rank"] = ifelse(
                length(tBgttA) == 0,
                0,
                1 * length(tBgttA)
            ) +
                ifelse(length(tBletA_sBeq0) == 0, 0, sum(rA / tBletA_sBeq0)) +
                ifelse(length(tBeqtA_sBeq1) == 0, 0, 0.5 * length(tBeqtA_sBeq1))
        }
        if (sA == 0) {
            tBgetA_sBeq0 <- mat_guanrank[
                vect >= tA & vecs == 0,
                "survival_rate"
            ]
            tBgetA_sBeq1 <- mat_guanrank[
                vect >= tA & vecs == 1,
                "survival_rate"
            ]
            tBlttA_sBeq0 <- mat_guanrank[vect < tA & vecs == 0, "survival_rate"]
            mat_guanrank[i, "rank"] = ifelse(
                length(tBgetA_sBeq0) == 0,
                0,
                sum(1 - 0.5 * tBgetA_sBeq0 / rA)
            ) +
                ifelse(
                    length(tBgetA_sBeq1) == 0,
                    0,
                    sum(1 - tBgetA_sBeq1 / rA)
                ) +
                ifelse(
                    length(tBlttA_sBeq0) == 0,
                    0,
                    sum(0.5 * rA / tBlttA_sBeq0)
                )
        }
    }
    rank <- mat_guanrank[, "rank"]
    # 0.5 is the correction for self-comparison
    # normalization to [0,1]
    mat_guanrank[, "rank"] <- (rank - 0.5) / max(rank)
    if (!complete) {
        mat_guanrank <- mat_guanrank[, c("time", "status", "rank")]
    }

    mat_guanrank
}

## example:
#tbl=data.frame(time=c(430,257,185,298,506),status=c(0,1,0,1,0))
#rownames(tbl)=paste0("patient",1:5)
#guanrank(tbl)
#guanrank_elf(tbl)
