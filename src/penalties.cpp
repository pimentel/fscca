#include "penalties.h"


LassoPenalty::LassoPenalty(double lam) :
    lambda_(lam)
{
    if (lam < 0.0)
    {
        forward_exception_to_r(
                std::runtime_error("Lasso lambda must be >= 0.0")
                );
    }
    lambda_ = lam;
}


void LassoPenalty::w(const arma::vec &x, arma::vec &result) const
{
    result = 1 / (arma::abs(x) + DELTA);
}

double LassoPenalty::lambda() const
{
    return lambda_;
}

