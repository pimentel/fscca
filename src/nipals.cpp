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

    arma::vec a1();
    arma::vec b1;

    arma::vec u1;
    arma::vec v1_old;

    while (eps > MIN_EPS)
    {
        a1 = arma::trans(X) * v1;
        a1 = a1 / arma::sum( v1 % v1 );
        a1 = a1 / sqrt( arma::sum( a1 % a1 ) );

        u1 = X * a1;

        b1 = arma::trans(Y) * u1;
        b1 = b1 / arma::sum( u1 % u1 );
        b1 = b1 / sqrt( arma::sum( b1 % b1 ) );

        v1_old = v1;

        v1 = Y * b1;
        eps = max(abs(v1 - v1_old));
    }

    arma::vec rho1 = arma::trans(u1) * v1;

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho1") = rho1,
            Rcpp::Named("u1") = u1,
            Rcpp::Named("v1") = v1,
            Rcpp::Named("a1") = a1,
            Rcpp::Named("b1") = b1
            );

    return result;
}
