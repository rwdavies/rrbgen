#include <Rcpp.h>
using namespace Rcpp;


//' @export
// [[Rcpp::export]]
int rcpp_return_same_int(const int hello) {
    // do something here
    int a = hello;
    return(a);
}

