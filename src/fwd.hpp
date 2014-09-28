#ifndef FSCCA_FWD_H
#define FSCCA_FWD_H

#include <RcppArmadillo.h>

#define NIPALS_EPS_CONVERGE 1e-6
#define S_NIPALS_EPS_CONVERGE 1e-3

#define ROUND_ZERO 5e-5

#define DELTA 1.0e-08

double square(double x)
{
    return x * x;
}

double l2_norm_sq(const arma::vec &x)
{
    double res = 0.0;

    for (arma::vec::const_iterator i = x.begin(); i != x.end(); ++i)
    {
        res += square(*i);
    }

    return res;
}

double l2_norm(const arma::vec &x)
{
    return sqrt( l2_norm_sq(x) );
}

#endif // FSCCA_FWD_H
