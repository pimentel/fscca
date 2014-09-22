#ifndef PENALTIES_HPP
#define PENALTIES_HPP

class NipalsPenalty {
    public:
        // Given a vector, returns the derivative of the penalty on that
        // function
        virtual double dpenalty(arma::vec &) const = 0;
};

class LassoPenalty : NipalsPenalty {
    public:
        LassoPenalty(double );
        virtual double dpenalty(arma::vec&) const;
    protected:
        double lambda;
};

LassoPenalty::LassoPenalty(double lam) :
    lambda(lam)
{}

double LassoPenalty::dpenalty(arma::vec& x) const
{
    // TODO: compute penalty

    return 42.0;
}

#endif // PENALTIES_HPP
