// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List NMF_optimized(const arma::mat &X, int K, int maxiter = 2000,
                         double tol = 1e-5) {

  const double eps = 2.2204e-256;

  const int nr = X.n_rows;
  const int nc = X.n_cols;

  // Input validation
  if (K <= 0 || K >= std::min(nr, nc)) {
    Rcpp::stop("K must be positive and smaller than both dimensions of X");
  }

  if (any(vectorise(X < 0))) {
    Rcpp::stop("X must be non-negative");
  }

  // Initialize W and H with random uniform values
  arma::mat W = arma::randu<arma::mat>(nr, K);
  arma::mat H = arma::randu<arma::mat>(K, nc);

  double old_eucl = arma::datum::inf;
  double eucl_dist = 0.0;
  int iter = 1;

  // Pre-allocate
  arma::mat WtW(K, K);
  arma::mat WtX(K, nc);
  arma::mat HHt(K, K);
  arma::mat XHt(nr, K);
  arma::mat WH(nr, nc);

  // Main optimization loop
  for (iter = 1; iter <= maxiter; iter++) {
    // Update H
    WtW = W.t() * W; // K x K
    WtX = W.t() * X; // K x nc
    H %= WtX / (WtW * H + eps);

    // Update W
    HHt = H * H.t(); // K x K
    XHt = X * H.t(); // nr x K
    W %= XHt / (W * HHt + eps);

    // Compute loss (skip first iteration)
    if (iter > 1) {
      arma::mat residual = X - W * H;
      eucl_dist = arma::norm(residual, "fro") * arma::norm(residual, "fro");

      double d_eucl = std::abs(eucl_dist - old_eucl);

      // Check convergence
      if (d_eucl < tol) {
        break;
      }

      old_eucl = eucl_dist;
    }
  }

  // Return results as a list
  return Rcpp::List::create(Rcpp::Named("W") = W, Rcpp::Named("H") = H,
                            Rcpp::Named("iter") = iter - 1,
                            Rcpp::Named("loss") = eucl_dist);
}

// [[Rcpp::export]]
List select_K_optimized(const arma::mat &X, int K_max = 20,
                        int repeat_times = 10, int maxiter = 2000,
                        bool verbose = true) {

  // K values to test: 2 to K_max
  int n_K = K_max - 1; // Number of K values (2, 3, ..., K_max)

  // Matrix to store reconstruction losses for each K and repeat
  arma::mat dist_K(n_K, repeat_times);
  dist_K.fill(arma::datum::nan);
  // Vectors for statistics
  arma::vec eii(n_K);
  eii.fill(arma::datum::nan); // 初始化为NA
  arma::vec row_means(n_K);
  row_means.fill(arma::datum::nan); // 初始化为NA

  int optimal_K = 2;
  bool found_optimal = false;

  // Main loop over K values
  for (int Ki = 2; Ki <= K_max; Ki++) {
    int Ki_idx = Ki - 2; // 转换为0-based索引: K=2 -> idx=0, K=3 -> idx=1, etc.

    // Run NMF multiple times with different initializations
    for (int Kj = 0; Kj < repeat_times; Kj++) {
      List res_ij = NMF_optimized(X, Ki, maxiter, 1e-5);

      arma::mat W = Rcpp::as<arma::mat>(res_ij["W"]);
      arma::mat H = Rcpp::as<arma::mat>(res_ij["H"]);

      // Compute reconstruction error: ||X - WH||_F^2
      arma::mat diff_matrix = X - W * H;
      dist_K(Ki_idx, Kj) = std::pow(arma::norm(diff_matrix, "fro"), 2);
    }

    // Compute mean loss for this K
    row_means(Ki_idx) = arma::mean(dist_K.row(Ki_idx));

    if (verbose) {
      Rcout << "loss of " << Ki << ": " << row_means(Ki_idx) << std::endl;
    }

    // Skip calculation for K = 2 (first iteration)
    if (Ki == 2) {
      continue;
    }

    // Calculate elbow improvement index
    // Ki_idx对应当前K，Ki_idx-1对应前一个K
    double current_loss = row_means(Ki_idx);
    double numerator = row_means(Ki_idx - 1) - current_loss;
    double denominator = row_means(0) - current_loss; // row_means(0)对应K=2

    // Check for valid denominator
    eii(Ki_idx) = numerator / denominator;

    // Check stopping criteria
    if (numerator <= 0) {
      optimal_K = Ki - 1; // 选择前一个K作为最优
      found_optimal = true;
      break;
    }

    if (eii(Ki_idx) < 0.05) {
      optimal_K = Ki - 1; // 选择前一个K作为最优
      found_optimal = true;
      break;
    }

    // 如果没有break，当前K就是最优（但会在下次循环被覆盖）
    optimal_K = Ki;
  }

  // 如果没有找到最优K（没有触发break），使用最后一个K
  if (!found_optimal) {
    optimal_K = K_max;
  }

  // Create K_all vector for output consistency
  arma::vec K_all_vec = arma::regspace(2, K_max);

  // Return results
  return List::create(Named("K") = optimal_K, Named("eii") = eii,
                      Named("mean_loss") = row_means,
                      Named("all_losses") = dist_K, Named("K_all") = K_all_vec);
}

double compute_loss(const arma::sp_mat &X, const arma::mat &W,
                    const arma::mat &H, const arma::sp_mat &S,
                    const arma::sp_mat &L, double alpha, double alpha_2) {
  // Loss 1: ||X - W*H||_F^2
  const arma::mat WH = W * H;
  // ||X - WH||_F^2 = ||X||_F^2 + ||WH||_F^2 - 2*trace(X' * WH)
  double x_norm_sq = arma::norm(X, "fro");
  x_norm_sq *= x_norm_sq;

  double wh_norm_sq = arma::norm(WH, "fro");
  wh_norm_sq *= wh_norm_sq;

  // trace(X' * WH) = sum(X .* WH)
  double trace_term = 0.0;
  arma::sp_mat::const_iterator it = X.begin();
  arma::sp_mat::const_iterator it_end = X.end();
  for (; it != it_end; ++it) {
    trace_term += (*it) * WH(it.row(), it.col());
  }

  double loss1 = x_norm_sq + wh_norm_sq - 2.0 * trace_term;

  // Loss 2: alpha * ||S * W||_F^2
  const arma::mat SW = S * W; // nr x K
  double loss2_norm = arma::norm(SW, "fro");
  double loss2 = alpha * loss2_norm * loss2_norm;

  // Loss 3: alpha_2 * sum(H .* (H * L))
  // H * L is K x nc (dense * sparse)
  const arma::mat HL = H * L; // K x nc
  // Element-wise multiplication H .* HL and sum all elements
  double loss3 = alpha_2 * arma::accu(H % HL); // Sum of element-wise product

  return loss1 + loss2 + loss3;
}

// [[Rcpp::export]]
List scAB_inner(const arma::sp_mat &X, const arma::sp_mat &A,
                const arma::sp_mat &D, const arma::sp_mat &L,
                const arma::sp_mat &S, int K, double alpha = 0.005,
                double alpha_2 = 0.005, int maxiter = 2000,
                double convergence_threshold = 1e-5) {

  // Constants
  const double eps = 2.2204e-256;
  const int nr = X.n_rows;
  const int nc = X.n_cols;

  // Initialize W and H with random uniform values
  arma::mat W = arma::randu<arma::mat>(nr, K);
  arma::mat H = arma::randu<arma::mat>(K, nc);

  // Pre-compute
  arma::sp_mat SS = S * S;
  // arma::mat A_mat = arma::mat(A);
  // arma::mat D_mat = arma::mat(D);
  // arma::mat L_mat = arma::mat(L);
  arma::mat WtW, HHt, WtX, XHt, SW, WHHt, SSW, HAt, HDt;

  // Variables for convergence check
  double old_eucl = 0.0;
  double eucl_dist = 0.0;
  int final_iter = 0;

  // Main iteration loop
  for (int iter = 1; iter <= maxiter; iter++) {
    // Update H
    WtW = W.t() * W;
    WtX = W.t() * X;

    HAt = (A * H.t()).t(); // K x nc
    HDt = (D * H.t()).t(); // K x nc

    arma::mat h_numerator = WtX + alpha_2 * HAt;             // K x nc
    arma::mat h_denominator = WtW * H + alpha_2 * HDt + eps; // K x nc

    H %= h_numerator / h_denominator;

    // Update W
    HHt = H * H.t();
    SW = S * W;
    XHt = X * H.t();
    WHHt = W * HHt;
    SSW = SS * W;

    arma::mat w_denominator = WHHt + alpha * SSW + eps; // nr x K

    W %= XHt / w_denominator;

    // Check convergence
    if (iter > 1) {
      eucl_dist = compute_loss(X, W, H, S, L, alpha, alpha_2);
      double d_eucl = std::abs(eucl_dist - old_eucl);

      if (d_eucl < convergence_threshold) {
        final_iter = iter;
        old_eucl = eucl_dist;
        break;
      }
      old_eucl = eucl_dist;
    } else {
      old_eucl = compute_loss(X, W, H, S, L, alpha, alpha_2);
    }

    final_iter = iter;
  }

  // Return results as a list
  return Rcpp::List::create(Rcpp::Named("W") = W, Rcpp::Named("H") = H,
                            Rcpp::Named("iter") = final_iter,
                            Rcpp::Named("loss") = old_eucl);
}
