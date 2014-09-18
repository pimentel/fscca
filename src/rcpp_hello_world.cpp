#include <Rcpp.h>

//' Hello world
//'
//' @return a list with foo bar
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_hello_world() {

    Rcpp::CharacterVector x = Rcpp::CharacterVector::create( "foo", "bar" )  ;
    Rcpp::NumericVector y   = Rcpp::NumericVector::create( 0.0, 1.0 ) ;
    Rcpp::List z            = Rcpp::List::create( x, y ) ;

    return z ;
}

//' First column
//'
//' @param x a matrix
//' @return the first column of the matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector firstColumn(Rcpp::NumericMatrix nm, int c)
{
    return nm.column(c);
}
