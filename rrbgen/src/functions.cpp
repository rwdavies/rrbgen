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
  std::uint32_t const_2_bit = pow(2, B_bit_prob) - 1;
  Rcpp::RawVector output(2 * N * B_bit_prob_divide_8);
  int c = 0;
  Rcpp::NumericVector v(3);
  Rcpp::NumericVector v2(3);
  std::uint32_t v3[3];
  double F = 0;
  int i_sample, i_col, F_int;
  IntegerVector idx(3);
  short i;
  F_int = 0;    
  for(i_sample =0; i_sample < N; i_sample++) {
    if (NumericVector::is_na(gp_sub(i_sample, 0))) {
      for(i_col = 0; i_col <= 2; i_col++) {
	v3[i_col] = 0;
      }
    } else {
      //
      // need to normalize gp_sub into integers here
      //
      F = 0;
      for(i_col =0; i_col <= 2; i_col++) {    
	v(i_col) = gp_sub(i_sample, i_col) * const_2_bit;
	v2(i_col) = v(i_col) - floor(v(i_col));
        F = F + v2(i_col);
      }
      F_int = 0.1 + std::round(F); // this is insane. dammit C++ sometimes. so round will get it almost there, then add 0.1 and floor
      //std::cout << "TEST1 - F=" << F << ",F_int=" << F_int << ",std::round(F)=" << std::round(F) << std::endl;            
      //
      if ((F_int == 0) || (F_int == 3)) {
	for(i_col =0; i_col <= 2; i_col++) {
	  v3[i_col] = v[i_col];
	}
      }
      if (F_int > 0) {
	idx = seq_along(v2) - 1;
	std::sort(idx.begin(), idx.end(), [&](int i, int j){return v2[i] <= v2[j];});
	for(i_col =0; i_col <= 2; i_col++) {
	  if (i_col >= (3 - F_int)) {
	    //std::cout << "CEILING i_col=" << i_col << ",idx[i_col]=" << idx[i_col] << std::endl;
	    v3[idx[i_col]] = ceil(v[idx[i_col]]);
	  } else {
	    //std::cout << "FLOOR i_col=" << i_col << ",idx[i_col]=" << idx[i_col] << std::endl;	    
	    v3[idx[i_col]] = floor(v[idx[i_col]]);
	  }
	}
      }
    }
    //
    // argh this was ugly
    //std::cout << "i_sample=" << i_sample << "gp_sub[i_sample, ]=" << gp_sub(i_sample, 0) << "," <<  gp_sub(i_sample, 1) << "," << gp_sub(i_sample, 2) << ",F=" << F << ",F_int=" << F_int << ", v=" << v << ", v2=" << v2 << ", v3=" << v3[0] << "," << v3[1] << "," << v3[2] << ",idx=" << idx << std::endl;
    //
    for(i_col =0; i_col < 2; i_col++) {
      std::uint32_t number = v3[i_col];
      //std::cout << "i_sample=" << i_sample << ",i_col=" << i_col << ",number=" << number << std::endl;
      for(i = B_bit_prob_divide_8_minus_1; i >= 0; i--) {
	// use c as counter as this is in order. do not change for loop order
	output[c + i] = get_nybble(number, i);
      }
      c = c + B_bit_prob_divide_8;
    }
  }
  return(output);
}
