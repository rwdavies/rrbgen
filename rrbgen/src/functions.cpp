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
Rcpp::RawVector rcpp_make_raw_data_vector_for_probabilities(Rcpp::NumericMatrix& gp_sub, int B_bit_prob = 16)
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



// basically, a copied and pasted of the above
// but, operate on a different ordering of things

//' @export
// [[Rcpp::export]]
void rcpp_place_gp_t_into_output(
    const Rcpp::NumericMatrix& gp_t,
    Rcpp::RawMatrix& to_out,
    const int i_sample,
    const int nSNPs,
    const int B_bit_prob = 16
) {
  int B_bit_prob_divide_8 = (B_bit_prob / 8);
  int B_bit_prob_divide_8_minus_1 = (B_bit_prob_divide_8 - 1);
  std::uint32_t const_2_bit = pow(2, B_bit_prob) - 1;
  int c = 0;
  Rcpp::NumericVector v(3);
  Rcpp::NumericVector v2(3);
  std::uint32_t v3[3];
  double F = 0;
  int i_SNP, i_col, F_int;
  IntegerVector idx(3);
  short i;
  F_int = 0;
  // now, for offset into to_out:
  //   i_sample is 1-based
  //   each sample takes up 2 * B_bit_prob_divide_8
  int offset = 2 * B_bit_prob_divide_8 * (i_sample - 1);
  for(i_SNP = 0; i_SNP < nSNPs; i_SNP++) {
    if (NumericVector::is_na(gp_t(0, i_SNP))) {
      for(i_col = 0; i_col <= 2; i_col++) {
	v3[i_col] = 0;
      }
    } else {
      //
      // need to normalize gp_sub into integers here
      //
      F = 0;
      for(i_col =0; i_col <= 2; i_col++) {    
	v(i_col) = gp_t(i_col, i_SNP) * const_2_bit;
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
    //std::cout << "i_SNP=" << i_SNP << "gp_sub[i_SNP, ]=" << gp_sub(i_SNP, 0) << "," <<  gp_sub(i_SNP, 1) << "," << gp_sub(i_SNP, 2) << ",F=" << F << ",F_int=" << F_int << ", v=" << v << ", v2=" << v2 << ", v3=" << v3[0] << "," << v3[1] << "," << v3[2] << ",idx=" << idx << std::endl;
    //
    c = 0; // re-set counter here
    for(i_col =0; i_col < 2; i_col++) {
        std::uint32_t number = v3[i_col];
	//std::cout << "i_SNP=" << i_SNP << ",i_col=" << i_col << ",number=" << number << std::endl;
	for(i = B_bit_prob_divide_8_minus_1; i >= 0; i--) {
	  // use c as counter as this is in order. do not change for loop order
	  to_out(offset + c + i, i_SNP) = get_nybble(number, i);
	}
	c = c + B_bit_prob_divide_8;	
    }
  }
  return;
}



// see http://www.cplusplus.com/forum/general/45432/


//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_convert_raw_probabilities_to_double_probabilities(const Rcpp::RawVector& data_raw_for_probs, int N, int B_bit_prob, Rcpp::LogicalVector& is_missing) {
  Rcpp::NumericMatrix gen_probs(N, 3);
  int B_bit_prob_divide_8 = B_bit_prob / 8;
  char hom_ref_vals[] =  {0, 0, 0, 0}; //some char/byte array
  char het_vals[] =  {0, 0, 0, 0}; //some char/byte array  
  //std::uint32_t denom = pow(2, B_bit_prob) - 1;
  double denomD = pow(2, B_bit_prob) - 1;  
  std::uint32_t x_hom_ref, x_het;
  int last_byte_used = 0;
  int i_sample;
  int i = 0;
  double p_hom_ref, p_het, p_hom_alt;
  for(i_sample = 0; i_sample < N; i_sample++) {
    if (is_missing[i_sample] == false) {
      for(i=0; i < B_bit_prob_divide_8; i++) {
	hom_ref_vals[i]=data_raw_for_probs[last_byte_used + i];
	het_vals[i]=data_raw_for_probs[last_byte_used + i + B_bit_prob_divide_8];
      }
      // 
      std::uint32_t *plong_hom_ref = (std::uint32_t*) hom_ref_vals;
      std::uint32_t *plong_het = (std::uint32_t*) het_vals;
      x_hom_ref = plong_hom_ref[0];
      p_hom_ref = x_hom_ref / denomD;
      x_het = plong_het[0];
      p_het = x_het / denomD;
      p_hom_alt = 1 - p_hom_ref - p_het;
      gen_probs(i_sample, 0) = p_hom_ref;
      gen_probs(i_sample, 1) = p_het;
      gen_probs(i_sample, 2) = p_hom_alt;    
      //std::cout << "i_sample=" << i_sample << ",x_hom_ref=" << x_hom_ref << ",p_hom_ref=" << p_hom_ref << std::endl;
      //std::cout << "i_sample=" << i_sample << ",x_het=" << x_het << ",p_het=" << p_het << std::endl;    
      //std::cout << "hom_ref_vals:" << hom_ref_vals[0] << std::endl;
    } else {
      gen_probs(i_sample, 0) = NA_REAL;
      gen_probs(i_sample, 1) = NA_REAL;
      gen_probs(i_sample, 2) = NA_REAL;
    }
    last_byte_used = last_byte_used + B_bit_prob_divide_8 * 2;
  }
  return(gen_probs);
}




  
