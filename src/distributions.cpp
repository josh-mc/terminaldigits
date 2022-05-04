#include <Rcpp.h>
using namespace Rcpp;

//First choose the distribution.

//[[Rcpp::export]]

NumericVector dist_in(int distribution,
                      int set_n,
                      double set_mean,
                      double set_sd)  {

  if(distribution == 1)  {

    return rnorm(set_n, set_mean, set_sd);

  }

  if(distribution == 2)  {

    return runif(set_n, set_mean, set_sd);

  }

  //3 == exponential

  else {

    return rexp(set_n, set_mean);

  }

}



//[[Rcpp::export]]

NumericVector data_in(int distribution,
                      int set_n,
                      double set_mean,
                      double set_sd,
                      int duplicates)  {

  NumericVector sample_0 = dist_in(distribution, set_n, set_mean, set_sd);

  //Setting up the addition of repititions.

  if(duplicates == 0)  {

    return sample_0;

  }

  else  {

    NumericVector sample_0_0 = rep_len(sample_0, duplicates);

    NumericVector sample_0_1(set_n + duplicates);

    //Combining the two vectors.

    for(int i = 0; i < set_n; ++i)  {

      sample_0_1(i) = sample_0(i);
    }

    for(int i = 0; i < duplicates; ++i)  {

      sample_0_1(i + set_n) = sample_0_0(i);
    }

    return sample_0_1;

  }

}



//Violation type

//Currently not used since focusing on simple violation (type 1 below)

//[[Rcpp::export]]

NumericVector violation(NumericVector sample_0,
                        int duplicates,
                        int type)  {

  if(type == 0)  {

    NumericVector sample_0_0 = rep(sample_0(0), duplicates);

    return sample_0_0;

  }

  else if(type == 1)  {

    NumericVector sample_0_0 = rep_len(sample_0, duplicates);

    return sample_0_0;

  }

  else  {

    NumericVector sample_0_0 = rpois(duplicates, 2);
    NumericVector sample_0_1(duplicates);

    for(int i = 0; i < duplicates; ++i)  {

      sample_0_1(i) = sample_0(sample_0_0(i));

    }

    return sample_0_1;

  }
}
