#include "cv.h"

//' @export
// [[Rcpp::export]]
arma::vec split_in_groups(size_t length, size_t k)
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

//' @export
// [[Rcpp::export]]
Rcpp::List groups_to_rows(const arma::vec& x, size_t k)
{

    Rcpp::List fit, test;

    for (size_t grp = 0; grp < k; ++grp)
    {
        Rcpp::NumericVector fit_vec;
        Rcpp::NumericVector test_vec;
        for (size_t j = 0; j < x.n_rows; ++j)
        {
            if ( x[j] == grp)
                test_vec.push_back(j);
            else
                fit_vec.push_back(j);
        }
        fit.push_back( fit_vec );
        test.push_back( test_vec );
    }

    return Rcpp::List::create( Rcpp::Named("fit", fit),
            Rcpp::Named("test", test));
}

//' @export
// [[Rcpp::export]]
Rcpp::List ret_list(Rcpp::NumericVector x)
{
    Rcpp::List a_list;
    for (size_t i = 0; i < 5; ++i)
        a_list.push_back( x );
    return a_list;
}

void split_matrix(const arma::mat& X, int k )
{

}
