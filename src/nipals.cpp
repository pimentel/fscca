#include <RcppArmadillo.h>

#include "fwd.hpp"
#include "penalties.hpp"

//' @useDynLib fscca


//' NIPALS CCA algorithm
//'
//' @param x a matrix x that has been centered and scaled
//' @param y a matrix y that has been centered and scaled
//' @return a list containing a1 and b1
//' @export
// [[Rcpp::export]]
Rcpp::List nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr) 
{

    if (Xr.nrow() != Yr.nrow())
    {
        forward_exception_to_r(
                std::runtime_error("nrows of X and Y must be equal!")
                );
    }

    double eps = 1.0;

    arma::mat X = Rcpp::as<arma::mat>(Xr);
    arma::mat Y = Rcpp::as<arma::mat>(Yr);

    arma::vec v1(Y.begin(), Y.n_rows, true);

    arma::vec a1(X.n_cols);
    arma::vec b1(Y.n_cols);

    arma::vec u1(X.n_rows);
    arma::vec v1_old(Y.n_rows);

    while (eps > NIPALS_EPS_CONVERGE)
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


//' Sparse NIPALS CCA algorithm
//'
//' @param x a matrix x that has been centered and scaled
//' @param y a matrix y that has been centered and scaled
//' @param lamx a positive enalty on 'a'
//' @param lamy a positive penalty on 'b'
//' @return a list containing a1 and b1
//' @export
// [[Rcpp::export]]
Rcpp::List sparse_nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr,
        double lamx, double lamy)
{
    if (lamx < 0.0 || lamy < 0.0)
    {
        forward_exception_to_r(
                std::runtime_error("lamx and lamy must be at least 0.0")
                );
    }

    // nipals will check the dimensions of X and Y
    Rcpp::List init_res = nipals(Xr, Yr);

    size_t p = Xr.ncol();
    size_t q = Yr.ncol();
    size_t n_iter = 0;

    double eps = 1;

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < 100) )
    {
        // TODO: write me
        ++n_iter;
        n_iter = 1000;
    }

    return init_res;
}

void iterate_sparse_nipals(arma::mat &Z, arma::vec &c)
{
    double eps = 1;
    size_t n_iter = 0;

    arma::vec abs_c;

    LassoPenalty np(3);

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < 50) )
    {
        // for now only compute LASSO
        // TODO: make this a call a functor for penalty derivative

        break;
    }

}


