#include "fwd.h"

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
