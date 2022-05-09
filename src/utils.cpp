#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

List int_dec(NumericVector s,
             int decimals) {

  // Truncating

  NumericVector s_1 = s * pow(10, decimals);

  IntegerVector s_f = as<IntegerVector>(s_1);

  // Pulling out the "integers" (everything except the last digit)


  IntegerVector int_1 = s_f / pow(10, 1);

  IntegerVector int_f = int_1 * 10;


  // Subtracting integers from sample to get decimals

  NumericVector dec_f = Rcpp::abs(s_f - int_f);
  IntegerVector dec_g = as<IntegerVector>(dec_f);

  List out = List::create(Named("s_f") = s_f,
                          Named("int_f") = int_1,
                          Named("dec_f") = dec_g);

  return out;

}

// [[Rcpp::export]]

IntegerVector observed_vec(IntegerVector u_int,
                           IntegerVector u_dec,
                           IntegerVector u_sam,
                           IntegerVector tab_sam) {

  IntegerVector neg_dec = Rcpp::clone(u_dec).sort(true);

  int n_int = u_int.size();
  int n_dec = u_dec.size();
  int n_sam = u_sam.size();
  int n = n_int * n_dec;
  IntegerVector out(n);

  int count = 0;
  int cc = 0;

  for(int i = 0; i < n_int; ++i)  {

    for(int j = 0; j < n_dec; ++j)  {

      int k = 0;

      if(u_int(i) >= 0)  {

        k = (u_int(i) * 10) + u_dec(j);

      }

      // So that with negative integers -1 and the decimal 2 leads to -12, not -8

      else {

        k = (u_int(i) * 10) - neg_dec(j);

      }

      if (cc >= n_sam) {

        out(count) = 0;
      }

      else if (k == u_sam(cc)) {

        out(count) = tab_sam(cc);

        cc = cc + 1;

      }

      else {

        out(count) = 0;

      }

      count = count + 1;

      }
  }


    return out;
  }






// [[Rcpp::export]]

IntegerVector full_vec(IntegerVector int_1,
                       IntegerVector dec_1)  {

  IntegerVector int_uni_x = sort_unique(int_1);
  IntegerVector dec_uni_x = sort_unique(dec_1);
  IntegerVector int_vec_x = rep_each(int_uni_x, dec_uni_x.size());
  IntegerVector dec_vec_x = rep(dec_uni_x, int_uni_x.size());
  IntegerVector full_vec = int_vec_x + dec_vec_x;

  return(full_vec);

}


// [[Rcpp::export]]

NumericVector actual_frac(IntegerVector int_full,
                          IntegerVector dec,
                          int new_n)  {

  IntegerVector uni_1 = sort_unique(dec);
  IntegerVector mat_0 = match(int_full, uni_1);
  IntegerVector mat_1 = ifelse(is_na(mat_0), 0, mat_0);

  IntegerVector tab_1 = table(dec);
  NumericVector vector_1_bins(int_full.size());

  for(int i = 0; i < int_full.size(); ++i)  {

    if(mat_1(i) > 0) {

      vector_1_bins(i) = tab_1(mat_1(i) - 1);

    }

    else {

      vector_1_bins(i) = 0;

    }

  }

  NumericVector vec_1_frac = vector_1_bins / new_n;

  return vec_1_frac;

}

// [[Rcpp::export]]

IntegerVector tab_it(IntegerVector x, int bins, int a, int b) {

  IntegerVector cts(bins);

  int t;

  for(int i = a; i < b; i++) {

    t = x[i] - 1;

    if (0 <= t && t < bins)

      cts[t]++;
  }

  return cts;
}

// [[Rcpp::export]]

IntegerVector perm_vector(IntegerVector v, IntegerVector r_sum, IntegerVector c_sum)  {

  IntegerVector s_v = sample(v, v.size());

  IntegerVector q(r_sum.size() * c_sum.size());

  int count = 0;
  int count_b = 0;

  for(int i = 0; i < r_sum.size(); i++) {

    int l = r_sum(i);

    IntegerVector out = tab_it(s_v, c_sum.size(), count, count + l);

    for(int j = 0; j < out.size(); j++) {

      q(j + count_b) = out(j);

    }

    count += l;
    count_b += out.size();

  }

  return q;

}


// [[Rcpp::export]]

NumericVector expected_cells(NumericVector r_frac,
                             NumericVector c_frac)  {

  NumericVector int_vec = rep_each(r_frac, c_frac.size());
  NumericVector dec_vec = rep(c_frac, r_frac.size());

  return int_vec * dec_vec;

}


// [[Rcpp::export]]

IntegerVector out_vector_cpp(IntegerVector c_sums)  {

  std::vector<int> mat_draw;

  mat_draw.reserve(sum(c_sums));

  for(int i = 0; i < c_sums.size(); ++i){

    for(int j = 0; j < c_sums[i]; ++ j){

      //This for some reason would not take a (). Thus [] above.

      mat_draw.push_back(i);

    }
  }

  return wrap(mat_draw);

}



// [[Rcpp::export]]

List RCONT(int n,
           IntegerVector r_sum,
           IntegerVector c_sum)  {

  List return_list(n);

  IntegerVector v = out_vector_cpp(c_sum) + 1;

  for(int i = 0; i < n; ++i)  {

    IntegerVector x = perm_vector(v,
                                  r_sum,
                                  c_sum);

    return_list[i] = x;

  }

  return(return_list);

}

