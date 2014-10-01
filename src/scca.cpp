#include "scca.h"

//' @export
// [[Rcpp::export]]
Rcpp::List fscca(arma::mat& X, arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        size_t k_folds,
        const arma::vec& lamx, const arma::vec& lamy)
{
    unsigned int best_x;
    unsigned int best_y;
    double best_avg_cv;

    cross_validation_alt(X, Y, penalty_x, penalty_y, k_folds, lamx, lamy,
            best_x, best_y, best_avg_cv);

    return Rcpp::List::create(Rcpp::Named("lamx", lamx[best_x]),
            Rcpp::Named("lamy", lamy[best_y]),
            Rcpp::Named("cv_cor", best_avg_cv));
}
