#ifndef CV_H
#define CV_H

#include "fwd.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include <RcppArmadillo.h>

#include "nipals.h"
#include "penalties.h"

typedef std::unique_ptr< arma::uvec > arma_uvec_ptr;

void groups_to_row_ptr(arma::uvec& x, size_t k,
        std::vector< arma_uvec_ptr >& fit, std::vector< arma_uvec_ptr >& test);

void cross_validation_alt(const arma::mat& X, const arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        size_t k_folds,
        const arma::vec& lamx, const arma::vec& lamy,
        unsigned int& best_lamx_idx, unsigned int& best_lamy_idx,
        double& best_avg_cv);

#endif // CV_H
