#include "cv.h"

//' @export
// [[Rcpp::export]]
arma::vec split_in_groups(size_t length, size_t k )
{
    std::unique_ptr<size_t[]> groups( new size_t[ k ] );
    size_t n_group = length / k;
    size_t last_group = length - (n_group * (k-1));

    for (size_t i = 0; i < k; ++i)
        groups[i] = n_group;

    // add leftovers in one-by-one
    for (size_t i = 0; i < last_group - n_group; ++i)
        ++groups[i];

    arma::vec x(length);

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

void split_matrix(const arma::mat& X, int k )
{

}
