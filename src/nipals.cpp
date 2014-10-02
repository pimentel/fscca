#include "nipals.h"

//' @useDynLib fscca

// [[Rcpp::depends(RcppArmadillo)]]

//' NIPALS CCA algorithm
//'
//' @param X a matrix X that has been centered and scaled
//' @param Y a matrix Y that has been centered and scaled
//' @return a list containing a, b, u, v, and rho
//' @export
// [[Rcpp::export]]
Rcpp::List nipals(const arma::mat& X, const arma::mat& Y)
{
    arma::vec a(X.n_cols);
    arma::vec b(Y.n_cols);

    arma::vec u(X.n_rows);
    arma::vec v(Y.n_rows);

    nipals_(X, Y, a, b, u, v);

    double rho = nipals_cov(u, v);

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho") = rho,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("a") = a,
            Rcpp::Named("b") = b
            );

    return result;
}

// Requires a, b, u, and v to be allocated
void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v)
{
    if (X.n_rows != Y.n_rows)
    {
        forward_exception_to_r(
                std::runtime_error("nrows of X and Y must be equal!")
                );
    }

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
//' @param X a matrix X that has been centered and scaled
//' @param Y a matrix Y that has been centered and scaled
//' @param penalty_x A character string of type "lasso"
//' @param penalty_y A character string of type "lasso"
//' @param lamx a positive penalty on 'a'
//' @param lamy a positive penalty on 'b'
//' @return a list containing a, b, u, v, and rho (covariance)
//' @export
// [[Rcpp::export]]
Rcpp::List sparse_nipals(const arma::mat& X, const arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        double lamx, double lamy)
{

    std::unique_ptr< NipalsPenalty > np_x = 
        PenaltyFactory::make_penalty( penalty_x, lamx );
    std::unique_ptr< NipalsPenalty > np_y = 
        PenaltyFactory::make_penalty( penalty_y, lamy );

    arma::vec a(X.n_cols), b(Y.n_cols), u(X.n_rows), v(Y.n_rows);

    // nipals will check the dimensions of X_ and Y_
    sparse_nipals_(X, Y, a, b, u, v, *np_x, *np_y);

    double rho = nipals_cov(u, v);

    Rcpp::List result = Rcpp::List::create(Rcpp::Named("rho") = rho,
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("a") = a,
            Rcpp::Named("b") = b
            );

    return result;
}

// In: X, Y, penalty_x, penalty_y
// Out: a, b, u, v. Requires these variables to be preallocated
void sparse_nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v,
        const NipalsPenalty &penalty_x,
        const NipalsPenalty &penalty_y)
{
    size_t p = X.n_cols;
    size_t q = Y.n_cols;
    size_t n_iter = 0;
    const size_t MAX_ITER = 100;

    double eps = 1.0;

    nipals_(X, Y, a, b, u, v);

    arma::vec v_prev(Y.n_rows);

    size_t a_zeros;
    size_t b_zeros;

    while ( (eps > S_NIPALS_EPS_CONVERGE) &&
            (n_iter < MAX_ITER) )
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

    if ( n_iter == MAX_ITER )
    {
        Rcpp::Rcout << "No convergence (sparse_nipals)." << std::endl;
    }

    a = a / l2_norm( a );
    b = b / l2_norm( b );

    u = X * a;
    v = Y * b;

}


// XXX: *IMPORTANT* Does not update 'u' at the end. Must compute u = Z * coef
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
        np.w( coef, w );

        Ztv = arma::trans(Z) * v;

        v_norm = l2_norm_sq( v );

        // TODO: make sure this is a deep copy
        coef_prev = coef;

        for (size_t i = 0; i < coef.n_elem; ++i)
        {
            coef[i] = Ztv[i] / (v_norm + np.lambda() * w[i]);
        }

        eps = arma::max( abs(coef - coef_prev) );
    }

    coef = coef / l2_norm( coef );
}

double nipals_cov(const arma::vec& u, const arma::vec& v)
{
    arma::mat res = arma::trans( u ) * v;
    return res(0, 0) / (u.n_rows - 1);
}


double nipals_cov(const arma::mat& X, const arma::vec& a,
        const arma::mat& Y, const arma::vec& b)
{
    arma::vec u = X * a;
    arma::vec v = Y * b;
    return nipals_cov( u, v );
}
