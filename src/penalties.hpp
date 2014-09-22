class NipalsPenalty {
    public:
        // Given a vector, returns the derivative of the penalty on that
        // function
        virtual double dpenalty(arma::vec &) = 0;
};

class LassoPenalty {
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
    return 42.0;
}

