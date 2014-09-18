
#include <Rcpp.h>
using namespace Rcpp;


//' @useDynLib sccaf


//' Hello world
//'
//' @return a list with foo bar
//' @export
// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}
