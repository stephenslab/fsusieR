#include <Rcpp.h>
using namespace Rcpp;

// FUNCTION DECLARATIONS
// ---------------------
NumericVector colSumsCpp (NumericMatrix x);
NumericVector colVarsCpp (NumericMatrix x);
NumericMatrix covCpp     (NumericMatrix x);

// FUNCTION DEFINITIONS
// --------------------
// [[Rcpp::export]]
NumericVector colSumsCpp(NumericMatrix x) {
  int n = x.nrow(), p = x.ncol();
  NumericVector sums(p);
  for (int j = 0; j < p; j++) {
    double s = 0;
    for (int i = 0; i < n; i++) s += x(i, j);
    sums[j] = s;
  }
  return sums;
}

// [[Rcpp::export]]
NumericVector colVarsCpp(NumericMatrix x) {
  int n = x.nrow(), p = x.ncol();
  NumericVector vars(p);
  for (int j = 0; j < p; j++) {
    double mean = 0.0, M2 = 0.0;
    for (int i = 0; i < n; i++) {
      double val = x(i, j);
      double delta = val - mean;
      mean += delta / (i + 1);
      M2 += delta * (val - mean);
    }
    vars[j] = M2 / (n - 1);
  }
  return vars;
}

// [[Rcpp::export]]
NumericMatrix covCpp(NumericMatrix x) {
  int n = x.nrow(), p = x.ncol();
  NumericVector means(p);
  NumericMatrix covmat(p, p);
  for (int j = 0; j < p; j++) {
    double s = 0;
    for (int i = 0; i < n; i++) s += x(i, j);
    means[j] = s / n;
  }
  for (int j = 0; j < p; j++) {
    for (int k = j; k < p; k++) {
      double s = 0;
      for (int i = 0; i < n; i++) {
        s += (x(i, j) - means[j]) * (x(i, k) - means[k]);
      }
      double cov = s / (n - 1);
      covmat(j, k) = cov;
      if (j != k) covmat(k, j) = cov;
    }
  }
  return covmat;
}
