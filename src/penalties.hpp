#ifndef PENALTIES_HPP
#define PENALTIES_HPP

#include <RcppArmadillo.h>

// TODO: Expose C++ classes to R and unit test

class NipalsPenalty {
    public:
        // Given a vector, returns the W(x) function 
        virtual arma::vec w(arma::vec &) const = 0;
};

class LassoPenalty : NipalsPenalty {
    public:
        LassoPenalty(double );
        virtual arma::vec w(arma::vec&) const;
    protected:
        double lambda;
};

LassoPenalty::LassoPenalty(double lam) :
    lambda(lam)
{}

arma::vec LassoPenalty::w(arma::vec& x) const
{
    arma::vec res;
    res = abs(x) + 1.0e-8;


    return res;
}

#endif // PENALTIES_HPP
