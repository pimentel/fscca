#include "fwd.h"

#include <RcppArmadillo.h>

//' @useDynLib sccaf


//' NIPALS CCA algorithm
//'
//' @param x a matrix x that has been centered and scaled
//' @param y a matrix y that has been centered and scaled
//' @return a list containing a1 and b1
//' @export
// [[Rcpp::export]]
Rcpp::List nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr) 
{
    double eps = 1.0;

    arma::mat X = Rcpp::as<arma::mat>(Xr);
    arma::mat Y = Rcpp::as<arma::mat>(Yr);

    arma::vec v1(Y.begin(), Yr.nrow(), true);

    arma::vec a1;
    arma::vec v1_2(v1.n_elem);

    double v1_2_sum = 0.0;
    // while (eps > MIN_EPS)
    // {
    v1_2 = v1 % v1;
    a1 = arma::trans(X) * v1;
    v1_2_sum = arma::sum(v1_2);
    //     eps = 1e-10;
    // }

    Rcpp::List c = Rcpp::List::create( v1, a1, v1_2, v1_2_sum);

    return c;
}
