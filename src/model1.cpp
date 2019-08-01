//
// Author: Anita Nandi
// Date: 2019-07-24
//
// Simple example of a TMB model - linear regression
//
// Data: prevalence survey data and covariate data
//
// The model: prev = invlogit(intercept + x * slope)
//

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // ------------------------------------------------------------------------ //
  // Input data
  // ------------------------------------------------------------------------ //
  
  // Covariates: Environmental data
  DATA_MATRIX(x);
  
  // Response: Prevalence point data
  DATA_VECTOR(positive_cases);
  DATA_VECTOR(examined_cases);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  Type priormean_intercept = -4.0;
  Type priorsd_intercept = 2.0; //priormean_intercept from data entry
  Type priormean_slope = 0.0;
  Type priorsd_slope = 0.5;
  
  // Priors for beta liklihood
  PARAMETER(log_point_sd);
  Type point_sd_mean = 0.1;
  Type point_sd_sd = 0.1;
  
  // Number of prev points
  int n_points = positive_cases.size();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  nll -= dnorm(exp(log_point_sd), point_sd_mean, point_sd_sd, true);
  
  Type nllpriors = nll;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_points);
  vector<Type> pixel_pred(n_points);
  vector<Type> reportnll(n_points);

  pixel_linear_pred = intercept + x * slope;
  pixel_pred = invlogit(pixel_linear_pred);
  
  for (int point = 0; point < n_points; point++) {
    nll -= dbinom(positive_cases[point], examined_cases[point], pixel_pred[point], true);
    reportnll[point] = -dbinom(positive_cases[point], examined_cases[point], pixel_pred[point], true);
  }

  REPORT(pixel_pred);
  REPORT(reportnll);
  REPORT(positive_cases);
  REPORT(examined_cases);
  REPORT(nllpriors);
  REPORT(nll);
  
  return nll;
  }
