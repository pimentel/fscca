// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// split_in_groups
arma::uvec split_in_groups(size_t length, size_t k);
RcppExport SEXP fscca_split_in_groups(SEXP lengthSEXP, SEXP kSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< size_t >::type length(lengthSEXP );
        Rcpp::traits::input_parameter< size_t >::type k(kSEXP );
        arma::uvec __result = split_in_groups(length, k);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// groups_to_rows
Rcpp::List groups_to_rows(const arma::uvec& x, size_t k);
RcppExport SEXP fscca_groups_to_rows(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::uvec& >::type x(xSEXP );
        Rcpp::traits::input_parameter< size_t >::type k(kSEXP );
        Rcpp::List __result = groups_to_rows(x, k);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// split_cv
Rcpp::List split_cv(size_t n_rows, size_t k);
RcppExport SEXP fscca_split_cv(SEXP n_rowsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< size_t >::type n_rows(n_rowsSEXP );
        Rcpp::traits::input_parameter< size_t >::type k(kSEXP );
        Rcpp::List __result = split_cv(n_rows, k);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// nipals
Rcpp::List nipals(const arma::mat& X, const arma::mat& Y);
RcppExport SEXP fscca_nipals(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP );
        Rcpp::List __result = nipals(X, Y);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// sparse_nipals
Rcpp::List sparse_nipals(const arma::mat& X, const arma::mat& Y, const std::string& penalty_x, const std::string& penalty_y, double lamx, double lamy);
RcppExport SEXP fscca_sparse_nipals(SEXP XSEXP, SEXP YSEXP, SEXP penalty_xSEXP, SEXP penalty_ySEXP, SEXP lamxSEXP, SEXP lamySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP );
        Rcpp::traits::input_parameter< const std::string& >::type penalty_x(penalty_xSEXP );
        Rcpp::traits::input_parameter< const std::string& >::type penalty_y(penalty_ySEXP );
        Rcpp::traits::input_parameter< double >::type lamx(lamxSEXP );
        Rcpp::traits::input_parameter< double >::type lamy(lamySEXP );
        Rcpp::List __result = sparse_nipals(X, Y, penalty_x, penalty_y, lamx, lamy);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// fscca
Rcpp::List fscca(arma::mat& X, arma::mat& Y, const std::string& penalty_x, const std::string& penalty_y, size_t k_folds, const arma::vec& lamx, const arma::vec& lamy);
RcppExport SEXP fscca_fscca(SEXP XSEXP, SEXP YSEXP, SEXP penalty_xSEXP, SEXP penalty_ySEXP, SEXP k_foldsSEXP, SEXP lamxSEXP, SEXP lamySEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< arma::mat& >::type Y(YSEXP );
        Rcpp::traits::input_parameter< const std::string& >::type penalty_x(penalty_xSEXP );
        Rcpp::traits::input_parameter< const std::string& >::type penalty_y(penalty_ySEXP );
        Rcpp::traits::input_parameter< size_t >::type k_folds(k_foldsSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type lamx(lamxSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type lamy(lamySEXP );
        Rcpp::List __result = fscca(X, Y, penalty_x, penalty_y, k_folds, lamx, lamy);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// get_submatrix
arma::mat get_submatrix(const arma::mat& X_, const arma::uvec& which_rows);
RcppExport SEXP fscca_get_submatrix(SEXP X_SEXP, SEXP which_rowsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< const arma::uvec& >::type which_rows(which_rowsSEXP );
        arma::mat __result = get_submatrix(X_, which_rows);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// get_submatrix_mult
arma::vec get_submatrix_mult(const arma::mat& X_, const arma::uvec& which_rows, const arma::vec& v);
RcppExport SEXP fscca_get_submatrix_mult(SEXP X_SEXP, SEXP which_rowsSEXP, SEXP vSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X_(X_SEXP );
        Rcpp::traits::input_parameter< const arma::uvec& >::type which_rows(which_rowsSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP );
        arma::vec __result = get_submatrix_mult(X_, which_rows, v);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// get_submatrix_mult_ptr
arma::vec get_submatrix_mult_ptr(const arma::mat& X, const arma::uvec& which_rows, const arma::vec& v);
RcppExport SEXP fscca_get_submatrix_mult_ptr(SEXP XSEXP, SEXP which_rowsSEXP, SEXP vSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP );
        Rcpp::traits::input_parameter< const arma::uvec& >::type which_rows(which_rowsSEXP );
        Rcpp::traits::input_parameter< const arma::vec& >::type v(vSEXP );
        arma::vec __result = get_submatrix_mult_ptr(X, which_rows, v);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// practice_arma
Rcpp::NumericVector practice_arma(arma::vec& x);
RcppExport SEXP fscca_practice_arma(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP );
        Rcpp::NumericVector __result = practice_arma(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
