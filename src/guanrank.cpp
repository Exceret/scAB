#include <Rcpp.h>
#include <cstring>
using namespace Rcpp;

// Guanrank for Survival Data (C++ Implementation)
// [[Rcpp::export]]
NumericVector guanrank_cpp(NumericMatrix mat_curve) {
  const int n = mat_curve.nrow();

  const double *time = &mat_curve(0, 0);
  const double *status = &mat_curve(0, 1);
  const double *survival_rate = &mat_curve(0, 2);

  NumericVector rank(n);
  double *rank_ptr = &rank[0];

  for (int i = 0; i < n; i++) {
    const double tA = time[i];
    const double rA = survival_rate[i];
    const int sA = static_cast<int>(status[i]);

    double sum = 0.0;

    if (sA == 1) {
      for (int j = 0; j < n; j++) {
        const double tB = time[j];
        const int sB = static_cast<int>(status[j]);
        const double rB = survival_rate[j];

        if (tB > tA) {
          sum += 1.0;
        } else if (tB == tA && sB == 1) {
          sum += 0.5;
        } else if (tB <= tA && sB == 0) {
          sum += rA / rB;
        }
      }
    } else {
      for (int j = 0; j < n; j++) {
        const double tB = time[j];
        const int sB = static_cast<int>(status[j]);
        const double rB = survival_rate[j];

        if (tB >= tA) {
          sum += (sB == 0) ? (1.0 - 0.5 * rB / rA) : (1.0 - rB / rA);
        } else if (sB == 0) {
          sum += 0.5 * rA / rB;
        }
      }
    }

    rank_ptr[i] = sum;
  }

  // Find max (serial, but fast)
  double max_rank = rank_ptr[0];
  for (int i = 1; i < n; i++) {
    if (rank_ptr[i] > max_rank)
      max_rank = rank_ptr[i];
  }

  for (int i = 0; i < n; i++) {
    rank_ptr[i] = (rank_ptr[i] - 0.5) / max_rank;
  }

  return rank;
}

// Complete Guanrank Function (C++ Backend)
// [[Rcpp::export]]
NumericMatrix guanrank_complete_cpp(NumericMatrix mat, bool complete = false) {
  // This function is called from R after km_curve preprocessing
  // mat should already be sorted and contain time, status, survival_rate

  const int n = mat.nrow();
  const int out_cols = complete ? 4 : 3;

  NumericMatrix result(n, out_cols);

  // Choose computation method based on size
  NumericVector ranks = guanrank_cpp(mat);

  // Assemble result
  if (complete) {
    std::memcpy(&result(0, 0), &mat(0, 0), n * 3 * sizeof(double));
    std::memcpy(&result(0, 3), &ranks[0], n * sizeof(double));
    colnames(result) =
        CharacterVector::create("time", "status", "survival_rate", "rank");
  } else {
    const double *mat_ptr = &mat(0, 0);
    double *result_ptr = &result(0, 0);
    const double *rank_ptr = &ranks[0];

    for (int i = 0; i < n; i++) {
      result_ptr[i * 3] = mat_ptr[i * 3];
      result_ptr[i * 3 + 1] = mat_ptr[i * 3 + 1];
      result_ptr[i * 3 + 2] = rank_ptr[i];
    }
    colnames(result) = CharacterVector::create("time", "status", "rank");
  }

  return result;
}