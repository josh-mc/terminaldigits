
#include <Rcpp.h>
#include "terminaldigits.h"

using namespace Rcpp;

//[[Rcpp::export]]

DataFrame perm_basic(int distribution,
                     int duplicates,
                     int set_n,
                     double set_mean,
                     double set_sd,
                     int decimals,
                     int reps,
                     int times,
                     double tolerance)  {

  NumericVector d_perm_p(times);
  NumericVector d_av_fre_p(times);
  NumericVector d_chi_p(times);
  NumericVector d_g2_p(times);
  NumericVector d_ft_p(times);
  NumericVector d_rms_p(times);

  for(int i = 0; i < times; ++i)  {

    NumericVector sample_0 = data_in(distribution, set_n, set_mean, set_sd, duplicates);

    int n = sample_0.size();

    // Truncating

    List a = int_dec(sample_0, decimals);

    IntegerVector sample_1 = a["s_f"];
    IntegerVector int_1 = a["int_f"];
    IntegerVector dec_1 = a["dec_f"];

    // Taking the unique length for the draw

    IntegerVector uni_d = unique(sample_1);

    int d_length = uni_d.size();

    // And the average frequency for the draw

    double av_fre_d = average_fre(sample_1);

    // Here we generate the expected fractions.
    // Row sums are integers. Column sums are decimals.\

    IntegerVector r_sums_t = table(int_1);
    IntegerVector c_sums_t = table(dec_1);

    NumericVector r_frac = as<NumericVector>(r_sums_t) / sum(r_sums_t);
    NumericVector c_frac = as<NumericVector>(c_sums_t) / sum(c_sums_t);

    // Calculating expected fraction.

    NumericVector expected_frac = expected_cells(r_frac, c_frac);

    int n_bins = expected_frac.size();

    // Take the full vector, i.e. all possible bins for each integer

    IntegerVector full_vec_1 = full_vec(int_1, dec_1);

    //Now take the actual fraction, a vector of the same length as the full vector
    //where the actual counts for each occupied bin are displayed as fractions.

    NumericVector actual_frac_d = actual_frac(full_vec_1, sample_1, n);


    double chi_gof_d = chisq_stat(n,
                                  actual_frac_d,
                                  expected_frac);


    double g2_gof_d = g2_stat(n,
                              actual_frac_d,
                              expected_frac);

    double ft_gof_d = ft_stat(n,
                              actual_frac_d,
                              expected_frac);


    double rms_gof_d = rms_stat(n,
                                actual_frac_d,
                                expected_frac);

    // Now setting up loops for permutation test

    double total_actual = 0;
    double av_fre_draw = 0;
    double chi_draw = 0;
    double g2_draw = 0;
    double ft_draw = 0;
    double rms_draw = 0;

    //Generating vector from which to draw.

    //Adding one below b/c tab_it doesn't count zeros.

    IntegerVector to_draw = out_vector_cpp(c_sums_t) + 1;

    for(int i = 0; i < reps; ++i)  {

      IntegerVector shuffled_num = perm_vector(to_draw, r_sums_t, c_sums_t);

      NumericVector actual_frac_x = as<NumericVector>(shuffled_num) / n;

      //O

      IntegerVector uni_perm = shuffled_num[shuffled_num > 0];

      int out = uni_perm.size();

      //AF

      double out_av = average_fre2(shuffled_num, n);

      //FT GoF

      double out_chi = chisq_stat(n,
                                  actual_frac_x,
                                  expected_frac);

      double out_g2 = g2_stat(n,
                              actual_frac_x,
                              expected_frac);

      double out_ft = ft_stat(n,
                              actual_frac_x,
                              expected_frac);

      double out_rms = rms_stat(n,
                                actual_frac_x,
                                expected_frac);



      if (out * (1 - tolerance) <= d_length )  {
        total_actual += 1;
      }


      if (out_av >= (1 - tolerance) * av_fre_d)  {
        av_fre_draw += 1;
      }


      //Chi

      if (out_chi >= (1 - tolerance) * chi_gof_d)  {
        chi_draw += 1;
      }


      //G2

      if (out_g2 >= (1 - tolerance) * g2_gof_d)  {
        g2_draw += 1;
      }


      //FT

      if (out_ft >= (1 - tolerance) * ft_gof_d)  {
        ft_draw += 1;
      }


      //RMS

      if (out_rms >= (1 - tolerance) * rms_gof_d)  {
        rms_draw += 1;
      }

    }

    d_perm_p(i) = (total_actual + 1) / (reps + 1);

    d_av_fre_p(i) = (av_fre_draw + 1) / (reps + 1);
    d_chi_p(i) = (chi_draw + 1) / (reps + 1);
    d_g2_p(i) = (g2_draw + 1) / (reps + 1);
    d_ft_p(i) = (ft_draw + 1) / (reps + 1);

    d_rms_p(i) = (rms_draw + 1) / (reps + 1);

    //This checks if in R the process has been stopped.

    Rcpp::checkUserInterrupt();

  }

  DataFrame out = DataFrame::create(Named("d_perm_p") = d_perm_p,
                                    Named("d_av_fre_p") = d_av_fre_p,
                                    Named("d_chi_p") = d_chi_p,
                                    Named("d_g2_p") = d_g2_p,
                                    Named("d_ft_p") = d_ft_p,
                                    Named("d_rms_p") = d_rms_p);

  return out;

}




//[[Rcpp::export]]

List terminal_independence(NumericVector x,
                           int decimals,
                           double reps,
                           double tolerance,
                           int type)  {

  int n = x.size();

  // Identifying terminal and preceding digits

  List a = int_dec(x, decimals);

  IntegerVector sample_1 = a["s_f"];
  IntegerVector int_1 = a["int_f"];
  IntegerVector dec_1 = a["dec_f"];

  // Here we generate the expected fractions.
  // Row sums are preceding digits. Column sums are terminal digits.

  IntegerVector r_sums_t = table(int_1);
  IntegerVector c_sums_t = table(dec_1);

  NumericVector r_frac = as<NumericVector>(r_sums_t) / n;
  NumericVector c_frac = as<NumericVector>(c_sums_t) / n;

  // Calculating expected fraction.

  NumericVector expected_frac = expected_cells(r_frac, c_frac);

  // Take the full vector, i.e. all possible bins for each preceding digit

  IntegerVector full_vec_1 = full_vec(int_1, dec_1);

  //Now take the actual fraction, a vector of the same length as the full vector
  //where the actual counts for each occupied bin are displayed as fractions.

  NumericVector actual_frac_d = actual_frac(full_vec_1, sample_1, n);


  auto stat_fun = [&](int type,
                      NumericVector actual_frac_d,
                      NumericVector expected_frac,
                      int n)  {

    switch(type) {

    case 2:

      return g2_stat(n,
                     actual_frac_d,
                     expected_frac);

      break;

    case 3:

      return ft_stat(n,
                     actual_frac_d,
                     expected_frac);

      break;

    case 4:

      return rms_stat(n,
                      actual_frac_d,
                      expected_frac);

    default:

      return chisq_stat(n,
                        actual_frac_d,
                        expected_frac);


    }



  };

  double STATISTIC = stat_fun(type,
                              actual_frac_d,
                              expected_frac,
                              n);

  // Now setting up loops for simulations

  double greater_equal = 0;

  //Generating vector from which to draw.

  //Adding one below b/c tab_it doesn't count zeros.

  IntegerVector to_draw = out_vector_cpp(c_sums_t) + 1;

  for(int i = 0; i < reps; ++i)  {

    if (i % 1000 == 0)  {

      Rcpp::checkUserInterrupt();

    }

    IntegerVector shuffled_num = perm_vector(to_draw, r_sums_t, c_sums_t);

    NumericVector actual_frac_x = as<NumericVector>(shuffled_num) / n;

    double simulated_statistic = stat_fun(type,
                                          actual_frac_x,
                                          expected_frac,
                                          n);

    if (simulated_statistic >= (1 - tolerance) * STATISTIC)  {
      greater_equal += 1;
    }


  }

  double P_VALUE = (greater_equal + 1) / (reps + 1);

  List out = List::create(Named("statistic") = STATISTIC,
                          Named("p_value") = P_VALUE);

  return out;

}

