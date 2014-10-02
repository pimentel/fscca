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

void scale_in_place(arma::mat& X, bool center, bool scale)
{
    if (center)
    {
        arma::rowvec mu = arma::mean(X, 0);
        for (size_t col = 0; col < mu.n_cols; ++col)
        {
            for (size_t row = 0; row < X.n_rows; ++row)
            {
                X(row, col) = X(row, col) - mu(col);
            }
        }
    }

    if (scale)
    {
        arma::rowvec sd = arma::stddev(X, 0, 0);
        for (size_t col = 0; col < sd.n_cols; ++col)
        {
            for (size_t row = 0; row < X.n_rows; ++row)
            {
                X(row, col) = X(row, col) / sd(col);
            }
        }
    }
}
