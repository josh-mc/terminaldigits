#include <Rcpp.h>
#include "terminaldigits.h"

using namespace Rcpp;

// [[Rcpp::export]]

double average_fre(IntegerVector x) {
  int vec_length = x.size();
  IntegerVector tab = table(x);
  IntegerVector tab_sq = tab * tab;
  double tab_sm = sum(tab_sq);
  double out = tab_sm / vec_length;
  return out;
}

// [[Rcpp::export]]

double average_fre2(IntegerVector x,
                    int n) {

  double a = sum(x * x);
  double out = a / n;

  return out;
}


// [[Rcpp::export]]

double ft_stat(int new_n,
               NumericVector vec_1_frac,
               NumericVector vec_2_frac)  {

  NumericVector x = pow(pow(vec_1_frac, 0.5) - pow(vec_2_frac, 0.5), 2);

  // Taking the sum of x
  double y = 0;

  int n = x.size();

  for(int i = 0; i < n; i++){
    y += x[i];
  }

  double ft = 4 * new_n * y;

  return ft;

}


// [[Rcpp::export]]

double rms_stat(int new_n,
                NumericVector vec_1_frac,
                NumericVector vec_2_frac)  {

  NumericVector x = pow((pow(new_n, 0.5) * (vec_1_frac - vec_2_frac)), 2);

  // Taking the sum of x
  double rms = 0;

  int n = x.size();

  for(int i = 0; i < n; i++){
    rms += x[i];
  }

  return rms;
}


// [[Rcpp::export]]

double chisq_stat(int draws,
                  NumericVector vec_1_frac,
                  NumericVector vec_2_frac)  {

  NumericVector x0 = pow((vec_1_frac - vec_2_frac), 2) / vec_2_frac;

  NumericVector x = x0[!is_na(x0)];

  // Taking the sum of x
  double y = 0;

  int n = x.size();

  for(int i = 0; i < n; i++){
    y += x[i];
  }

  double chi = y * draws;

  return chi;

}



// [[Rcpp::export]]

double g2_stat(int draws,
               NumericVector vec_1_frac,
               NumericVector vec_2_frac)  {

  NumericVector x0 = vec_1_frac * log(vec_1_frac / vec_2_frac);

  NumericVector x = x0[!is_na(x0)];

  // Taking the sum of x
  double y = 0;

  int n = x.size();

  for(int i = 0; i < n; i++){
    y += x[i];
  }

  double g2 = 2 * draws * y;

  return g2;

}




