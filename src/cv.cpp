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

typedef std::shared_ptr< arma::vec > arma_vec_ptr;

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector practice_arma(arma::vec & x)
{
    arma_vec_ptr z(new arma::vec(10));
    std::cout << "sup doodie ";
    std::cout << z->n_rows << std::endl;
    (*z)[0] = 2.0;
    std::cout << "hilowww " << z->at(0) << '\t' << z->at(2) << std::endl;

    Rcpp::NumericVector y;
    y.push_back(3.4);

    return y;
}

void groups_to_row_ptr(arma::vec& x, size_t k,
        std::vector< arma_vec_ptr >& fit, std::vector< arma_vec_ptr >& test)
{

    std::unique_ptr<size_t[]> groups( new size_t[ k ] );

    for (size_t i = 0; i < x.n_rows; ++i)
        ++groups[ x[i] ];

    for (size_t grp = 0; grp < k; ++grp)
    {
        arma_vec_ptr fit_vec( new arma::vec( x.n_rows - groups[grp] ) );
        arma_vec_ptr test_vec( new arma::vec( groups[grp] ) );
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

//' @export
// [[Rcpp::export]]
Rcpp::List split_cv(size_t n_rows, size_t k)
{
    arma::vec groups = split_in_groups(n_rows, k);
    std::vector< arma_vec_ptr > fit;
    std::vector< arma_vec_ptr > test;

    groups_to_row_ptr(groups, k, fit, test);

    Rcpp::List fit_list;
    Rcpp::List test_list;
    for (size_t i = 0; i < fit.size(); ++i)
        fit_list.push_back( *fit[i] );

    for (size_t i = 0; i < fit.size(); ++i)
        test_list.push_back( *test[i] );

    return Rcpp::List::create( Rcpp::Named("fit", fit_list),
            Rcpp::Named("test", test_list));
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
