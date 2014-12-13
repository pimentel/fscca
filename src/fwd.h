#ifndef FSCCA_FWD_H
#define FSCCA_FWD_H

// #if !defined(ARMA_64BIT_WORD)
//   #define ARMA_64BIT_WORD
// #endif

#include <memory>

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]

#define NIPALS_EPS_CONVERGE 1e-6
#define S_NIPALS_EPS_CONVERGE 1e-3

#define ROUND_ZERO 5e-5

#define DELTA 1.0e-08

double square(double x);

double l2_norm_sq(const arma::vec &x);

double l2_norm(const arma::vec &x);

void scale_in_place(arma::mat& X, bool center, bool scale);

void round_in_place(arma::vec& x, short precision);
void round_in_place(arma::mat& x, short precision);
double d_round(double value, short precision);

#endif // FSCCA_FWD_H
