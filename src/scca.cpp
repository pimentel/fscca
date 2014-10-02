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

//' @export
// [[Rcpp::export]]
Rcpp::List fscca(arma::mat& X, arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        const arma::vec& lamx, const arma::vec& lamy,
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
        Rcpp::Rcout << "Computing component: " << c << std::endl;

        cross_validation_alt(X, Y, penalty_x, penalty_y, k_folds,
                lamx, lamy, best_x, best_y, best_avg_cv);

        std::unique_ptr<NipalsPenalty> pen_x =
            PenaltyFactory::make_penalty( penalty_x, lamx[best_x] );
        std::unique_ptr<NipalsPenalty> pen_y =
            PenaltyFactory::make_penalty( penalty_y, lamy[best_y] );

        sparse_nipals_(X, Y, a, b, u, v, *pen_x, *pen_y);

        arma::rowvec l(2);
        l(0) = lamx[best_x];
        l(1) = lamy[best_y];

        X = X - ( u * arma::trans(a) );
        Y = Y - ( v * arma::trans(b) );

        A.col(c) = a;
        B.col(c) = b;
        U.col(c) = u;
        V.col(c) = v;
        L.row(c) = l;
        cov(c) = nipals_cov(u, v);

    }

    return Rcpp::List::create(
            Rcpp::Named("A", A),
            Rcpp::Named("B", B),
            Rcpp::Named("U", U),
            Rcpp::Named("V", V),
            Rcpp::Named("L", L),
            Rcpp::Named("cov", cov)
            );
}
