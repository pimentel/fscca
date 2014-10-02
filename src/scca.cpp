#include "scca.h"

Rcpp::List cv_alt_wrapper(arma::mat& X, arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        size_t k_folds,
        const arma::vec& lamx, const arma::vec& lamy)
{
    unsigned int best_x;
    unsigned int best_y;
    double best_avg_cv;

    cross_validation_alt(X, Y, penalty_x, penalty_y, k_folds, lamx, lamy,
            best_x, best_y, best_avg_cv);

    // TODO: once the 'best' have been found, use nipals_sparse to find the
    // solution

    return Rcpp::List::create(
            Rcpp::Named("lamx", lamx[best_x]),
            Rcpp::Named("lamy", lamy[best_y]),
            Rcpp::Named("cv_cor", best_avg_cv));
}

//' Fast sparse canonical correlation
//'
//' This is a fast implementation of sparse canonical correlation analysis.
//'
//' @param X a matrix of dimension n x p
//' @param Y a matrix of dimension n x q
//' @param penalty_x a string indicating the penalty function to use on matrix
//' X. Currently only "lasso" is implemented.
//' @param penalty_y a string indicating the penalty function to use on matrix
//' Y. Currently only "lasso" is implemented.
//' @param lam_x a numeric vector of tuning parameters on X
//' @param lam_y a numeric vector of tuning parameters on Y
//' @param k_folds an integer denoting the number of folds to use in cros-validation
//' @param n_components the number of components to compute
//' @param center center the columns to mean zero?
//' @param scale scale the columns to standard deviation one?
//' @return A list with matrices:
//' \describe{
//'     \item{A}{ A matrix of dimension p x n_components of the canonical
//' vector }
//'     \item{B}{ A matrix of dimension q x n_components of the canonical
//' vector }
//'     \item{U}{ A matrix of dimension n x n_components of the X * a }
//'     \item{V}{ A matrix of dimension n x n_components of the Y * b }
//'     \item{lambda}{ A matrix of dimension n_components x n of the optimal
//' tuning parameters}
//'     \item{covar}{ The average cross-validated covariance }
//' }
//' @export
// [[Rcpp::export]]
Rcpp::List fscca(arma::mat X, arma::mat Y,
        const std::string& penalty_x, const std::string& penalty_y,
        const arma::vec& lam_x, const arma::vec& lam_y,
        size_t k_folds = 5, size_t n_components = 1,
        bool center = true, bool scale = false)
{
    if (X.n_rows != Y.n_rows)
    {
        forward_exception_to_r(
                std::runtime_error("nrows of X and Y must be equal!")
                );
    }

    scale_in_place(X, center, scale);
    scale_in_place(Y, center, scale);

    arma::mat U(X.n_rows, n_components);
    arma::mat V(X.n_rows, n_components);
    arma::mat A(X.n_cols, n_components);
    arma::mat B(Y.n_cols, n_components);
    arma::mat L(n_components, 2);
    arma::vec cov(n_components);

    unsigned int best_x;
    unsigned int best_y;
    double best_avg_cv;

    arma::vec a(X.n_cols), b(Y.n_cols), u(X.n_rows), v(X.n_rows);

    for (size_t c = 0; c < n_components; ++c)
    {
        Rcpp::Rcout << "Computing component: " << c + 1 << std::endl;

        if (c > 0)
        {
            X = X - ( u * arma::trans(a) );
            Y = Y - ( v * arma::trans(b) );
        }

        cross_validation_alt(X, Y, penalty_x, penalty_y, k_folds,
                lam_x, lam_y, best_x, best_y, best_avg_cv);

        std::unique_ptr<NipalsPenalty> pen_x =
            PenaltyFactory::make_penalty( penalty_x, lam_x[best_x] );
        std::unique_ptr<NipalsPenalty> pen_y =
            PenaltyFactory::make_penalty( penalty_y, lam_y[best_y] );

        sparse_nipals_(X, Y, a, b, u, v, *pen_x, *pen_y);

        arma::rowvec l(2);
        l(0) = lam_x[best_x];
        l(1) = lam_y[best_y];

        A.col(c) = a;
        B.col(c) = b;
        U.col(c) = u;
        V.col(c) = v;
        L.row(c) = l;
        cov(c) = best_avg_cv;

    }

    Rcpp::Environment base("package:base");
    Rcpp::Function round = base["round"];

    return Rcpp::List::create(
            Rcpp::Named("A", round(A, 4)),
            Rcpp::Named("B", round(B, 4)),
            Rcpp::Named("U", round(U, 4)),
            Rcpp::Named("V", round(V, 4)),
            Rcpp::Named("lambda", L),
            Rcpp::Named("covar", cov)
            );
}
