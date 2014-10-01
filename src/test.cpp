#include <RcppArmadillo.h>

typedef std::shared_ptr< arma::uvec > arma_uvec_ptr;

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
    return X_.rows(which_rows) * v;
}

//' @export
// [[Rcpp::export]]
arma::vec get_submatrix_mult_ptr(const arma::mat& X,
        const arma::uvec& which_rows,
        const arma::vec& v)
{
    arma::mat Y = X.rows(which_rows);
    return  Y * v;
}
// arma::mat get_submatrix(Rcpp::NumericMatrix X, Rcpp::NumericVector rows)
//


//' @export
// [[Rcpp::export]]
Rcpp::NumericVector practice_arma(arma::vec & x)
{
    arma_uvec_ptr z(new arma::uvec(10));
    std::cout << "sup doodie ";
    std::cout << z->n_rows << std::endl;
    (*z)[0] = 2.0;
    std::cout << "hilowww " << z->at(0) << '\t' << z->at(2) << std::endl;

    Rcpp::NumericVector y;
    y.push_back(3.4);

    return y;
}
