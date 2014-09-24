#ifndef NIPALS_H
#define NIPALS_H

#include <RcppArmadillo.h>

#include "fwd.hpp"
#include "penalties.hpp"


Rcpp::List nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr);

void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v);

#endif // NIPALS_H
