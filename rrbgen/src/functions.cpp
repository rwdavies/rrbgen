#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;


//' @export
// [[Rcpp::export]]
int rcpp_return_same_int(const int hello) {
    // do something here
    int a = hello;
    return(a);
}



// https://codereview.stackexchange.com/questions/30593/split-a-long-integer-into-eight-4-bit-values

unsigned short get_nybble( std::uint32_t number, const unsigned short part )
{
  //if( part > 1 ) {
  //  std::cout << "You have selected " << part << std::endl;    
  //     throw std::out_of_range("'part' must be a number between 0 and 7");
  //}
    return (number >> (8 * part)) & 255;
}

//' @export
// [[Rcpp::export]]
Rcpp::RawVector rcpp_make_raw_data_vector_for_probabilities(Rcpp::NumericMatrix gp_sub, int B_bit_prob = 16)
{
  int N = gp_sub.nrow();
  int B_bit_prob_divide_8 = (B_bit_prob / 8);
  int B_bit_prob_divide_8_minus_1 = (B_bit_prob_divide_8 - 1);
  int const_2_bit = pow(2, B_bit_prob) - 1;
  Rcpp::RawVector output(2 * N * B_bit_prob_divide_8);
  int c = 0;
  for(int i_sample =0; i_sample < N; i_sample++) {
    for(int i_col =0; i_col < 2; i_col++) {
      std::uint32_t number = gp_sub(i_sample, i_col) * const_2_bit;
      std::cout << "i_sample=" << i_sample << ", i_col=" << i_col << ", gp_sub[i_sample, i_col]=" << gp_sub(i_sample, i_col) << ", number=" << number << std::endl;
      for( short i = B_bit_prob_divide_8_minus_1; i >= 0; i--) {
	// use c as counter as this is in order. do not change for loop order
        output[c] = get_nybble(number, i);
	c++;
      }
    }
  }
  return(output);
}
