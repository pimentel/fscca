#include <RcppArmadillo.h>

#include "fwd.hpp"
#include "penalties.hpp"

#include "nipals.h"

//' @useDynLib fscca

// Requires a, b, u, and v to be allocated
void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v)
{
    double eps = 1.0;

    v = Y.col(0);

    arma::vec v_prev(Y.n_rows);

    while (eps > NIPALS_EPS_CONVERGE)
    {
        a = arma::trans(X) * v;
        a = a / l2_norm_sq( v );
        a = a / l2_norm( a ) ;

        u = X * a;

        b = arma::trans(Y) * u;
        b = b / l2_norm_sq( u );
        b = b / l2_norm( b );

        v_prev = v;

        v = Y * b;
        eps = arma::max(abs(v - v_prev));
    }
}


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

    arma::mat X = Rcpp::as<arma::mat>(Xr);
    arma::mat Y = Rcpp::as<arma::mat>(Yr);

    // arma::vec v(Y.begin(), Y.n_rows, true);

    arma::vec a(X.n_cols);
    arma::vec b(Y.n_cols);

    arma::vec u(X.n_rows);
    arma::vec v(Y.n_rows);

    nipals_(X, Y, a, b, u, v);

    arma::mat rho = arma::trans(u) * v;

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho") = rho,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("a") = a,
            Rcpp::Named("b") = b
            );

    return result;
}


size_t count_zeros(const arma::vec &x)
{
    size_t count = 0;

    for (arma::vec::const_iterator it = x.begin(); it != x.end(); ++it )
    {
        if (*it < ROUND_ZERO)
            ++count;
    }

    return count;
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

    arma::mat X = Rcpp::as<arma::mat>(Xr);
    arma::mat Y = Rcpp::as<arma::mat>(Yr);

    arma::vec a(X.n_cols), b(Y.n_cols), u(X.n_rows), v(Y.n_rows);

    // nipals will check the dimensions of X and Y
    sparse_nipals_(X, Y, a, b, u, v, lamx, lamy);

    arma::vec rho = arma::trans(u) * v;

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho") = rho,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("a") = a,
            Rcpp::Named("b") = b
            );

    return result;
}

void sparse_nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v,
        double lamx, double lamy)
{
    size_t p = X.n_cols;
    size_t q = Y.n_cols;
    size_t n_iter = 0;

    double eps = 1.0;

    v = Y.col(0);

    arma::vec v_prev(Y.n_rows);

    LassoPenalty penalty_x(lamx);
    LassoPenalty penalty_y(lamy);

    size_t a_zeros;
    size_t b_zeros;

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < 100) )
    {
        iterate_sparse_nipals(X, a, u, v, penalty_x);
        u = X * a;

        iterate_sparse_nipals(Y, b, v, u, penalty_y);
        v_prev = v;
        v = Y * b;

        eps = arma::max(abs( v - v_prev ));

        a_zeros = count_zeros( a );
        b_zeros = count_zeros( b );

        if (a_zeros >= (p - 1) && b_zeros >= (q - 2))
            break;

        ++n_iter;
    }

    if ( n_iter == 100 )
    {
        Rcpp::Rcout << "No convergence" << std::endl;
    }

    a = a / l2_norm( a );
    b = b / l2_norm( b );

    u = X * a;
    v = Y * b;

}


void iterate_sparse_nipals(const arma::mat &Z, arma::vec &coef,
        arma::vec &u, const arma::vec &v, const NipalsPenalty& np)
{
    // Z: The data matrix
    // coef: 'a' in Xa
    // u: 'u' in u = Xa
    // v: 'v' in t(X) %*% v
    double eps = 1;
    size_t n_iter = 0;

    arma::vec w(coef.n_elem);

    arma::vec coef_prev(coef.n_elem);
    arma::vec Ztv(v.n_elem);

    double v_norm = 0;

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < 50) )
    {
        // for now only compute LASSO
        // TODO: make this a call a functor for penalty derivative
        np.w( coef, w );

        Ztv = arma::trans(Z) * v;

        v_norm = l2_norm_sq( v );

        // TODO: make sure this is a deep copy
        coef_prev = coef;

        for (size_t i = 0; i < coef.n_elem; ++i)
        {
            coef[i] = Ztv[i] / (v_norm + np.lambda() * w[i]);
        }

        // TODO: is max the correct thing to do here?
        eps = arma::max( abs(coef - coef_prev) );
    }

    coef = coef / l2_norm( coef );
    // u = Z * coef;
}
