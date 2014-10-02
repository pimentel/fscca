#ifndef NIPALS_H
#define NIPALS_H

#include <string>

#include <RcppArmadillo.h>

#include "fwd.h"
#include "penalties.h"

Rcpp::List nipals(const arma::mat &X, const arma::mat &Y);

void nipals_(const arma::mat &X, const arma::mat &Y,
        arma::vec &a, arma::vec &b,
        arma::vec &u, arma::vec &v);

Rcpp::List sparse_nipals(const arma::mat& X, const arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
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

double nipals_cov(const arma::mat& X, const arma::vec& a,
        const arma::mat& Y, const arma::vec& b);
double nipals_cov(const arma::vec& u, const arma::vec& v);

#endif // NIPALS_H
