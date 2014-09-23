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
        a1 = a1 / l2_norm_sq( v1 );
        a1 = a1 / l2_norm( a1 ) ;

        u1 = X * a1;

        b1 = arma::trans(Y) * u1;
        b1 = b1 / l2_norm_sq( u1 );
        b1 = b1 / l2_norm( b1 );

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

void iterate_sparse_nipals(const arma::mat &Z, arma::vec &coef,
        arma::vec &u, const arma::vec &v, const NipalsPenalty& np);

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

    // nipals will check the dimensions of X and Y
    Rcpp::List init_res = nipals(Xr, Yr);

    size_t p = Xr.ncol();
    size_t q = Yr.ncol();
    size_t n_iter = 0;

    double eps = 1;

    arma::mat X = Rcpp::as<arma::mat>(Xr);
    arma::mat Y = Rcpp::as<arma::mat>(Yr);

    arma::vec u(X.n_rows);
    arma::vec v(Y.begin(), Y.n_rows, true);
    arma::vec v_prev;

    arma::vec a = Rcpp::as<arma::vec>(init_res["a1"]);
    arma::vec b = Rcpp::as<arma::vec>(init_res["b1"]);

    LassoPenalty penalty_x(lamx);
    LassoPenalty penalty_y(lamy);

    size_t a_zeros;
    size_t b_zeros;

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < 100) )
    {
        // TODO: write me
        iterate_sparse_nipals(X, a, u, v, penalty_x);
        u = X * a;

        iterate_sparse_nipals(Y, b, v, u, penalty_y);
        v_prev = v;
        v = Y * b;

        eps = max(abs( v - v_prev ));

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
    arma::vec rho = arma::trans(u) * v;

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho") = rho,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("a") = a,
            Rcpp::Named("b") = b,
            Rcpp::Named("niter") = n_iter
            );

    return result;
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
        eps = max( abs(coef - coef_prev) );
    }

    coef = coef / l2_norm( coef );
    // u = Z * coef;
}


