#include <RcppArmadillo.h>

typedef std::shared_ptr< arma::uvec > arma_uvec_ptr;

//' @export
// [[Rcpp::export]]
arma::vec zero_mat(arma::mat& X)
{
    X.fill(0);
    return X.col(0);
}
