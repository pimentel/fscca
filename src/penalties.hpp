#ifndef PENALTIES_HPP
#define PENALTIES_HPP

#include <RcppArmadillo.h>

// TODO: Expose C++ classes to R and unit test

class NipalsPenalty {
    public:
        // Given a vector, returns the W(x) function 
        virtual void w(const arma::vec &x, arma::vec &result) const = 0;
        virtual double lambda() const = 0;
};

class LassoPenalty : public NipalsPenalty {
    public:
        LassoPenalty(double );
        virtual void w(const arma::vec &x, arma::vec &result) const;
        virtual double lambda() const;
    protected:
        double lambda_;
};

LassoPenalty::LassoPenalty(double lam) :
    lambda_(lam)
{}

void LassoPenalty::w(const arma::vec &x, arma::vec &result) const
{
    result = abs(x) + 1.0e-8;
}

double LassoPenalty::lambda() const
{
    return lambda_;
}

#endif // PENALTIES_HPP
