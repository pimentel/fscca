#ifndef NIPALS_H
#define NIPALS_H

#include <string>

#include <RcppArmadillo.h>

#include "fwd.hpp"
#include "penalties.hpp"

Rcpp::List nipals(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y);

void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v);

Rcpp::List sparse_nipals(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y,
        std::string penalty_x, std::string penalty_y,
        double lamx, double lamy);

void sparse_nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v,
        const NipalsPenalty &np_x,
        const NipalsPenalty &np_y
        );

void iterate_sparse_nipals(const arma::mat &Z, arma::vec &coef,
        arma::vec &u, const arma::vec &v, const NipalsPenalty &np);

size_t count_zeros(const arma::vec &x);


#endif // NIPALS_H
