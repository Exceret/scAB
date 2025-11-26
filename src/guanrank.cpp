#include <Rcpp.h>
using namespace Rcpp;

//' @title Guanrank for Survival Data (C++ Implementation)
//'
//' @description Fast C++ implementation of Guan's rank for survival data
//'
//' @param mat_curve A numeric matrix from km_curve output with columns:
//'                  time, status, survival_rate
//'
//' @return A numeric vector of Guan's rank values (normalized to [0,1])
//'
//' @export
// [[Rcpp::export]]
NumericVector guanrank_cpp(NumericMatrix mat_curve) {
    int n = mat_curve.nrow();
    
    // Extract columns
    NumericVector time = mat_curve(_, 0);
    NumericVector status = mat_curve(_, 1);
    NumericVector survival_rate = mat_curve(_, 2);
    
    // Initialize rank vector
    NumericVector rank(n, 0.0);
    
    // Compute ranks
    for (int i = 0; i < n; i++) {
        double tA = time[i];
        double rA = survival_rate[i];
        int sA = status[i];
        
        if (sA == 1) {
            // Event case
            double sum_tBgttA = 0.0;
            double sum_tBletA_sBeq0 = 0.0;
            double sum_tBeqtA_sBeq1 = 0.0;
            
            for (int j = 0; j < n; j++) {
                double tB = time[j];
                int sB = status[j];
                double rB = survival_rate[j];
                
                // tB > tA
                if (tB > tA) {
                    sum_tBgttA += 1.0;
                }
                
                // tB <= tA and sB == 0
                if (tB <= tA && sB == 0) {
                    sum_tBletA_sBeq0 += rA / rB;
                }
                
                // tB == tA and sB == 1
                if (tB == tA && sB == 1) {
                    sum_tBeqtA_sBeq1 += 0.5;
                }
            }
            
            rank[i] = sum_tBgttA + sum_tBletA_sBeq0 + sum_tBeqtA_sBeq1;
            
        } else {
            // Censored case
            double sum_tBgetA_sBeq0 = 0.0;
            double sum_tBgetA_sBeq1 = 0.0;
            double sum_tBlttA_sBeq0 = 0.0;
            
            for (int j = 0; j < n; j++) {
                double tB = time[j];
                int sB = status[j];
                double rB = survival_rate[j];
                
                // tB >= tA and sB == 0
                if (tB >= tA && sB == 0) {
                    sum_tBgetA_sBeq0 += 1.0 - 0.5 * rB / rA;
                }
                
                // tB >= tA and sB == 1
                if (tB >= tA && sB == 1) {
                    sum_tBgetA_sBeq1 += 1.0 - rB / rA;
                }
                
                // tB < tA and sB == 0
                if (tB < tA && sB == 0) {
                    sum_tBlttA_sBeq0 += 0.5 * rA / rB;
                }
            }
            
            rank[i] = sum_tBgetA_sBeq0 + sum_tBgetA_sBeq1 + sum_tBlttA_sBeq0;
        }
    }
    
    // Normalize to [0,1]
    double max_rank = max(rank);
    for (int i = 0; i < n; i++) {
        rank[i] = (rank[i] - 0.5) / max_rank;
    }
    
    return rank;
}


//' @title Complete Guanrank Function (C++ Backend)
//'
//' @param mat A matrix, data frame, or data table containing survival data
//' @param complete Logical indicating whether to include survival_rate column
//'
//' @return A numeric matrix with time, status, rank (and optionally survival_rate)
//'
//' @export
// [[Rcpp::export]]
NumericMatrix guanrank_complete_cpp(NumericMatrix mat, bool complete = false) {
    // This function is called from R after km_curve preprocessing
    // mat should already be sorted and contain time, status, survival_rate
    
    int n = mat.nrow();
    NumericVector rank_values = guanrank_cpp(mat);
    
    NumericMatrix result;
    
    if (complete) {
        result = NumericMatrix(n, 4);
        for (int i = 0; i < n; i++) {
            result(i, 0) = mat(i, 0);  // time
            result(i, 1) = mat(i, 1);  // status
            result(i, 2) = rank_values[i];  // rank
            result(i, 3) = mat(i, 2);  // survival_rate
        }
        colnames(result) = CharacterVector::create("time", "status", "rank", "survival_rate");
    } else {
        result = NumericMatrix(n, 3);
        for (int i = 0; i < n; i++) {
            result(i, 0) = mat(i, 0);  // time
            result(i, 1) = mat(i, 1);  // status
            result(i, 2) = rank_values[i];  // rank
        }
        colnames(result) = CharacterVector::create("time", "status", "rank");
    }
    
    return result;
}