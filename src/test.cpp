#include <RcppArmadillo.h>

//' @export
// [[Rcpp::export]]
arma::mat get_submatrix(const arma::mat& X_, const arma::uvec& which_rows)
{
    // arma::mat X_ = Rcpp::as<arma::mat>(X);
    // arma::uvec which_rows = Rcpp::as<arma::uvec>(rows);
    return X_.rows(which_rows);
}

//' @export
// [[Rcpp::export]]
arma::vec get_submatrix_mult(const arma::mat& X_,
        const arma::uvec& which_rows,
        const arma::vec& v)
{
    // arma::mat X_ = Rcpp::as<arma::mat>(X);
    // arma::uvec which_rows = Rcpp::as<arma::uvec>(rows);
    return X_.rows(which_rows) * v;
}

// arma::mat get_submatrix(Rcpp::NumericMatrix X, Rcpp::NumericVector rows)
