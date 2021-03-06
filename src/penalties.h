#ifndef PENALTIES_HPP
#define PENALTIES_HPP

#include <memory>

#include <RcppArmadillo.h>

#include "fwd.h"

// TODO: Unit test C++
// TODO: Include constructor with 'X' and 'a' for H30 and SCAD

class NipalsPenalty {
    public:
        // Given a vector, returns the W(x) function 
        virtual void w(const arma::vec &x, arma::vec &result) const = 0;
        virtual double lambda() const = 0;
        virtual ~NipalsPenalty(){};
};

class LassoPenalty : public NipalsPenalty {
    public:
        LassoPenalty(double );
        ~LassoPenalty() {}
        virtual void w(const arma::vec &x, arma::vec &result) const;
        virtual double lambda() const;
    protected:
        double lambda_;
};


class HLPenalty : public NipalsPenalty {
    public:
        HLPenalty(double );
        ~HLPenalty() {}
    protected:
        double lambda_;
};

class PenaltyFactory {
    public:
        static std::unique_ptr<NipalsPenalty> make_penalty(
                const std::string &type,
                double lam)
        {
            if (type.compare( "lasso" ) == 0)
            {
                return std::unique_ptr<NipalsPenalty>( new LassoPenalty(lam) );
            }

            forward_exception_to_r(
                    std::runtime_error("Invalid 'type' of penalty function.")
                    );
            return std::unique_ptr<NipalsPenalty>();
        }
};

#endif // PENALTIES_HPP
