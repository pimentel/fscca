#ifndef SCCA_H
#define SCCA_H

#include <RcppArmadillo.h>

#include "fwd.h"

#include "cv.h"
#include "nipals.h"
#include "penalties.h"

Rcpp::List fscca(arma::mat& X, arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        size_t k_folds,
        const arma::vec& lamx, const arma::vec& lamy);

#endif // SCCA_H
