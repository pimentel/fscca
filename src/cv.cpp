#include "cv.h"

//' @export
// [[Rcpp::export]]
arma::uvec split_in_groups(size_t length, size_t k)
{
    std::unique_ptr<size_t[]> groups( new size_t[ k ] );
    size_t n_group = length / k;
    size_t last_group = length - (n_group * (k-1));

    for (size_t i = 0; i < k; ++i)
        groups[i] = n_group;

    // add leftovers in one-by-one
    for (size_t i = 0; i < last_group - n_group; ++i)
        ++groups[i];

    arma::uvec x(length);

    size_t n_samp = 0;
    while (n_samp < length)
    {
        size_t samp = std::rand() % k;
        if ( groups[samp] > 0 )
        {
            --groups[samp];
            x[n_samp] = samp;
            ++n_samp;
        }
    }

    return x;
}

void groups_to_row_ptr(arma::uvec& x, size_t k,
        std::vector< arma_uvec_ptr >& fit, std::vector< arma_uvec_ptr >& test)
{

    std::unique_ptr<size_t[]> groups( new size_t[ k ] );

    for (size_t i = 0; i < x.n_rows; ++i)
        ++groups[ x[i] ];

    for (size_t grp = 0; grp < k; ++grp)
    {
        arma_uvec_ptr fit_vec( new arma::uvec( x.n_rows - groups[grp] ) );
        arma_uvec_ptr test_vec( new arma::uvec( groups[grp] ) );
        size_t fit_count = 0;
        size_t test_count = 0;
        for (size_t j = 0; j < x.n_rows; ++j)
        {
            if ( x[j] == grp)
            {
                (*test_vec)[test_count] = j;
                ++test_count;
            }
            else
            {
                (*fit_vec)[fit_count] = j;
                ++fit_count;
            }
        }
        fit.push_back( fit_vec );
        test.push_back( test_vec );
    }

}

// This function does the bulk of the CV
void cross_validation_alt(const arma::mat& X, const arma::mat& Y,
        const std::string& penalty_x, const std::string& penalty_y,
        size_t k_folds,
        const arma::vec& lamx, const arma::vec& lamy,
        unsigned int& best_lamx_idx, unsigned int& best_lamy_idx,
        double& best_avg_cv)
{

    // TODO: add seed for randomness
    
    // get rows for cross validation
    arma::uvec groups = split_in_groups(X.n_rows, k_folds);
    
    std::vector< arma_uvec_ptr > fit;
    std::vector< arma_uvec_ptr > test;
    
    // fit and test now contain k_folds arma_uvec_ptrs of rows for fit and test
    // set for each k
    groups_to_row_ptr(groups, k_folds, fit, test);

    size_t it = 0;
    double eps = 1.0;

    // initialize all the penalties
    std::vector< std::shared_ptr<NipalsPenalty> > pen_x;
    std::vector< std::shared_ptr<NipalsPenalty> > pen_y;

    for (size_t i = 0; i < lamx.n_rows; ++i)
        pen_x.push_back( PenaltyFactory::make_penalty(penalty_x, lamx[i]) );
    for (size_t i = 0; i < lamx.n_rows; ++i)
        pen_y.push_back( PenaltyFactory::make_penalty(penalty_y, lamy[i]) );

    unsigned int opt_x = std::rand() % pen_x.size();
    unsigned int opt_y = 0;

    NipalsPenalty* cur_pen_x;
    NipalsPenalty* cur_pen_y;
    cur_pen_x = pen_x[opt_x].get();

    arma::mat avg_cv_cor( pen_x.size(), pen_y.size(), arma::fill::zeros );

    double cur_max = 0.0, prev_max = 0.0;

    while (eps > 1e-02 & it < 20)
    {
        for (size_t i_y = 0; i_y < pen_y.size(); ++i_y)
        {
            cur_pen_y = pen_y[i_y].get();
            double cur_cv_cor = 0.0;
            for (size_t k = 0; k < k_folds; ++k)
            {
                arma::mat X_ = X.rows( *(fit[k]) );
                arma::mat Y_ = Y.rows( *(fit[k]) );
                arma::vec a(X_.n_cols), b(Y_.n_cols), u(X_.n_rows), v(X_.n_rows);
                sparse_nipals_(X_, Y_, a, b, u, v, *cur_pen_x, *cur_pen_y);
                arma::mat X_t = X.rows( *(test[k]) );
                arma::mat Y_t = Y.rows( *(test[k]) );

                cur_cv_cor += abs(nipals_cov(X_t, a, Y_t, b));
            }
            cur_cv_cor /= k_folds;

            avg_cv_cor(opt_x, i_y) = cur_cv_cor;
        }

        // gets the index of the optimal y
        avg_cv_cor.row( opt_x ).max( opt_y );

        cur_pen_y = pen_y[opt_y].get();

        for (size_t i_x = 0; i_x < pen_x.size(); ++i_x)
        {
            cur_pen_x = pen_x[i_x].get();
            double cur_cv_cor = 0.0;
            for (size_t k = 0; k < k_folds; ++k)
            {
                arma::mat X_ = X.rows( *(fit[k]) );
                arma::mat Y_ = Y.rows( *(fit[k]) );
                arma::vec a(X_.n_cols), b(Y_.n_cols), u(X_.n_rows), v(X_.n_rows);

                sparse_nipals_(X_, Y_, a, b, u, v, *cur_pen_x, *cur_pen_y);

                arma::mat X_t = X.rows( *(test[k]) );
                arma::mat Y_t = Y.rows( *(test[k]) );
                cur_cv_cor += abs(nipals_cov(X_t, a, Y_t, b));
            }

            cur_cv_cor /= k_folds;
            avg_cv_cor(i_x, opt_y) = cur_cv_cor;
        }

        // should have the optimum x now
        avg_cv_cor.col( opt_y ).max( opt_x );

        cur_max = avg_cv_cor(opt_x, opt_y);
        eps = abs( cur_max - prev_max );
        prev_max = cur_max;
        ++it;
    }

    best_lamx_idx = opt_x;
    best_lamy_idx = opt_y;
    best_avg_cv = avg_cv_cor(opt_x, opt_y);
}
