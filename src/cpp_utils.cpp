#include <RcppArmadillo.h>
// [[Rcpp::depends(BH)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nu_k_root_finder_equation_cpp(arma::cube data, int p, arma::vec z_k, int n_k, double nu_k, arma::mat Sigma_inv_k, int N) {

  // Sum of determinant terms
  double sum_log_det = 0;
  for (int i = 0; i < N; i++) {
    double sign = 0;
    double log_det_val = 0;
    arma::log_det(log_det_val, sign, data.slice(i) * Sigma_inv_k / 2.0);
    sum_log_det += z_k(i) * log_det_val;
  }

  // Sum of digamma terms
  double sum_digamma = 0;
  for (int j = 0; j< p; j++) {
    sum_digamma += R::digamma((nu_k - (j + 1) + 1) / 2.0);
  }
  sum_digamma *= n_k;

  // Final result
  double result = sum_log_det - sum_digamma;

  return result;
}
