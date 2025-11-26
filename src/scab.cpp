// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List NMF_optimized(const arma::mat& X, 
                         int K, 
                         int maxiter = 2000, 
                         double tol = 1e-5) {
    
    const double eps = 2.2204e-256;
    
    int nr = X.n_rows;
    int nc = X.n_cols;
    
    // Input validation
    if (K <= 0 || K >= std::min(nr, nc)) {
        Rcpp::stop("K must be positive and smaller than both dimensions of X");
    }
    
    if (X.min() < 0) {
        Rcpp::stop("X must be non-negative");
    }
    
    // Initialize W and H with random uniform values
    arma::mat W = arma::randu<arma::mat>(nr, K);
    arma::mat H = arma::randu<arma::mat>(K, nc);
    
    double old_eucl = arma::datum::inf;
    double eucl_dist = 0.0;
    int iter;
    
    // Main optimization loop
    for (iter = 1; iter <= maxiter; iter++) {
        // Update H
        arma::mat WtW = W.t() * W;           // K x K
        arma::mat WtX = W.t() * X;           // K x nc
        H = H % (WtX / (WtW * H + eps));
        
        // Update W
        arma::mat HHt = H * H.t();           // K x K
        arma::mat XHt = X * H.t();           // nr x K
        W = W % (XHt / (W * HHt + eps));
        
        // Compute loss (skip first iteration)
        if (iter > 1) {
            arma::mat diff = X - W * H;
            eucl_dist = arma::accu(arma::square(diff));
            
            double d_eucl = std::abs(eucl_dist - old_eucl);
            
            // Check convergence
            if (d_eucl < tol) {
                break;
            }
            
            old_eucl = eucl_dist;
        } 
    }
    
    // Return results as a list
    return Rcpp::List::create(
        Rcpp::Named("W") = W,
        Rcpp::Named("H") = H,
        Rcpp::Named("iter") = iter - 1,
        Rcpp::Named("loss") = eucl_dist
    );
}

// [[Rcpp::export]]
List select_K_optimized(const arma::mat& X,
                        int K_max = 20,
                        int repeat_times = 10,
                        int maxiter = 2000,
                        bool verbose = true) {    
    // K values to test: 2 to K_max
    int n_K = K_max - 1;  // Number of K values (2, 3, ..., K_max)
    
    // Matrix to store reconstruction losses for each K and repeat
    mat dist_K(n_K, repeat_times);
    dist_K.fill(datum::nan);
    
    // Vectors for statistics
    vec eii(n_K);
    eii.fill(datum::nan);  // 初始化为NA
    vec row_means(n_K);
    row_means.fill(datum::nan);  // 初始化为NA
    
    int optimal_K = 2;
    bool found_optimal = false;
    
    // Main loop over K values
    for (int Ki = 2; Ki <= K_max; Ki++) {
        int Ki_idx = Ki - 2;  // 转换为0-based索引: K=2 -> idx=0, K=3 -> idx=1, etc.
        
        // Run NMF multiple times with different initializations
        for (int Kj = 0; Kj < repeat_times; Kj++) {
            List res_ij = NMF_optimized(X, Ki, maxiter, 1e-5);
            
            // Extract W and H
            mat W = as<mat>(res_ij["W"]);
            mat H = as<mat>(res_ij["H"]);
            
            // Compute reconstruction error: ||X - WH||_F^2
            mat diff_matrix = X - W * H;
            dist_K(Ki_idx, Kj) = arma::accu(arma::square(diff_matrix));
        }
        
        // Compute mean loss for this K
        row_means(Ki_idx) = arma::mean(dist_K.row(Ki_idx));
        
        if (verbose) {
            Rcout << "loss of " << Ki << ": " 
                  << row_means(Ki_idx) << std::endl;
        }
        
        // Skip calculation for K = 2 (first iteration)
        if (Ki == 2) {
            continue;
        }
        
        // Calculate elbow improvement index
        // Ki_idx对应当前K，Ki_idx-1对应前一个K
        double numerator = row_means(Ki_idx - 1) - row_means(Ki_idx);
        double denominator = row_means(0) - row_means(Ki_idx);  // row_means(0)对应K=2
        
        // Check for valid denominator
        if (std::abs(denominator) < 1e-10 || !arma::is_finite(denominator)) {
            eii(Ki_idx) = 0.0;
        } else {
            eii(Ki_idx) = numerator / denominator;
        }
        
        // Check stopping criteria
        if (numerator <= 0) {
            optimal_K = Ki - 1;  // 选择前一个K作为最优
            found_optimal = true;
            break;
        }
        
        if (eii(Ki_idx) < 0.05) {
            optimal_K = Ki - 1;  // 选择前一个K作为最优
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
    vec K_all_vec = arma::regspace(2, K_max);
    
    // Return results
    return List::create(
        Named("K") = optimal_K,
        Named("eii") = eii,
        Named("mean_loss") = row_means,
        Named("all_losses") = dist_K,
        Named("K_all") = K_all_vec
    );
}



// [[Rcpp::export]]
List scAB_inner(
    const arma::mat& X,  
    const arma::sp_mat& A,
    const arma::sp_mat& D,
    const arma::sp_mat& L,
    const arma::sp_mat& S,
    int K,
    double alpha = 0.005,
    double alpha_2 = 0.005,
    int maxiter = 2000,
    double convergence_threshold = 1e-5
) {
    // Constants
    const double eps = 2.2204e-256;
    const int nr = X.n_rows;
    const int nc = X.n_cols;
    
    // Initialize W and H with random uniform values
    arma::mat W = arma::randu<arma::mat>(nr, K) * std::sqrt(arma::mean(arma::mean(X)) / K);
    arma::mat H = arma::randu<arma::mat>(K, nc) * std::sqrt(arma::mean(arma::mean(X)) / K);
    
    // Pre-compute 
    arma::sp_mat SS = S * S;
    arma::mat A_mat = arma::mat(A);
    arma::mat D_mat = arma::mat(D);
    arma::mat L_mat = arma::mat(L);

    arma::mat WtW(K, K);
    arma::mat HHt(K, K);
    arma::mat WtX(K, nc);
    arma::mat X_Ht(nr, K);
    arma::mat SW(nr, K);

    // Variables for convergence check
    double old_eucl = 0.0;
    double eucl_dist = 0.0;
    int final_iter = 0;
    
    // Loss function (inline computation)
    auto compute_loss = [&]() -> double {
        double loss1 = arma::accu(arma::square(X - W * H)); // Frobenius norm squared
        
        arma::mat SW = arma::mat(S * W);
        double loss2 = alpha * arma::accu(arma::square(SW));
        
        double loss3 = alpha_2 * arma::accu(H % (H * L_mat));
        
        return loss1 + loss2 + loss3;
    };

    // Main iteration loop
    for (int iter = 1; iter <= maxiter; iter++) {
        // Update H
        WtW = W.t() * W;
        WtX = W.t() * X;

        H %= (WtX + alpha_2 * (H * A_mat.t())) / (WtW * H + alpha_2 * (H * D_mat.t()) + eps);
        
        HHt = H * H.t();
        // Update W
        SW = arma::mat(S * W);
        X_Ht = X * H.t();
        arma::mat W_HHt = W * (H * H.t());
        
        W %= X_Ht / (W * HHt + alpha * (arma::mat(SS * W)) + eps);
        
        // Check convergence
        if (iter > 1) {
            eucl_dist = compute_loss();
            double d_eucl = std::abs(eucl_dist - old_eucl);
            
            if (d_eucl < convergence_threshold) {
                final_iter = iter;
                old_eucl = eucl_dist;
                break;
            }
            old_eucl = eucl_dist;
        } else {
            old_eucl = compute_loss();
        }
        
        final_iter = iter;
    }
    
    // Return results as a list
    return List::create(
        Named("W") = W,
        Named("H") = H,
        Named("iter") = final_iter,
        Named("loss") = old_eucl
    );
}

// Handle dgeMatrix
// [[Rcpp::export]]
List scAB_inner_wrapper(
    SEXP X_sexp,              // Matrix of type dgeMatrix
    const arma::sp_mat& A,
    const arma::sp_mat& D,
    const arma::sp_mat& L,
    const arma::sp_mat& S,
    int K,
    double alpha = 0.005,
    double alpha_2 = 0.005,
    int maxiter = 2000,
    double convergence_threshold = 1e-5
) {
    arma::mat X = Rcpp::as<arma::mat>(X_sexp);
    
    return scAB_inner(X, A, D, L, S, K, alpha, alpha_2, maxiter, convergence_threshold);
}