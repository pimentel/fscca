#ifndef NIPALS_H
#define NIPALS_H

#include <RcppArmadillo.h>

#include "fwd.hpp"
#include "penalties.hpp"

Rcpp::List nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr);

void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v);

Rcpp::List sparse_nipals(Rcpp::NumericMatrix Xr, Rcpp::NumericMatrix Yr,
        double lamx, double lamy);

void sparse_nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v,
        double lamx, double lamy);

void iterate_sparse_nipals(const arma::mat &Z, arma::vec &coef,
        arma::vec &u, const arma::vec &v, const NipalsPenalty& np);

size_t count_zeros(const arma::vec &x);


#endif // NIPALS_H
