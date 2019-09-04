//
// Author: Anita Nandi
// Date: 2019-07-24
//
// Simple example of a TMB model - linear regression with field
//
// Data: prevalence survey data and covariate data
//
// The model: prev = invlogit(intercept + x * slope + logit_prevalence_field.array())
//

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // ------------------------------------------------------------------------ //
  // Spatial field data
  // ------------------------------------------------------------------------ //
  
  // The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
  DATA_SPARSE_MATRIX(Apixel);
  DATA_STRUCT(spde, spde_t);
  
  // ------------------------------------------------------------------------ //
  // Input data
  // ------------------------------------------------------------------------ //
  
  // Covariates: Environmental data
  DATA_MATRIX(x);
  
  // Response: Prevalence point data
  DATA_VECTOR(positive_cases);
  DATA_VECTOR(examined_cases);
  
  DATA_VECTOR(incidence_cases);
  DATA_VECTOR(incidence_pop);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // Prevalence to incidence relationship parameters
  DATA_VECTOR(prev_inc_par);
  PARAMETER(log_prev_inc_extra_slope);
  DATA_SCALAR(priormean_log_prev_inc_extra_slope); 
  DATA_SCALAR(priorsd_log_prev_inc_extra_slope);
  
  Type prev_inc_extra_slope = exp(log_prev_inc_extra_slope);
  

  // spde hyperparameters
  // sigma defines strength of random field. 
  // rho defines range (corr = 0.1 after rho space). 
  //   Might need to switch to log.
  PARAMETER(log_sigma);
  PARAMETER(log_rho);
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);


  // Priors on spde hyperparameters
  //   kappa -- i.e. exp(priormean_log_kappa) -- set as approximately the width of the region being studied.
  //   This implies prior belief in a fairly flat field.
  //   tau -- exp(priormean_log_tau) -- set to close to zero. Betas on regression coefficients have priors of 0 so this is reasonable.
  //Type priormean_log_kappa = -3;
  //Type priorsd_log_kappa   = 0.5;
  //Type priormean_log_tau   = -0.50;
  //Type priorsd_log_tau     = 2.0;
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);

  // Convert hyperparameters to natural scale
  // todo
  Type kappa = sqrt(8) / rho;
  Type nu = 1;
  // nu = 1
  // gamma(nu + d / 2) = gamma(2) = 1
  // gamma(nu) = 1
  // So gamma(nu + d / 2)(4pi)^{d/2) / gamma(nu) is just 4 * pi
  Type tau = sigma * pow(kappa, nu) * sqrt(4 * M_PI);


  PARAMETER_VECTOR(nodemean);
  
  // Number of prev and inc points
  int n_prev_points = positive_cases.size();
  int n_inc_points = incidence_cases.size();
  int n_points = n_prev_points + n_inc_points;
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  nll -= dnorm(log_prev_inc_extra_slope, priormean_log_prev_inc_extra_slope, priorsd_log_prev_inc_extra_slope, true);
  
  // Likelihood of hyperparameter of polygon iid random effect.
  Type lambda = -log(prior_iideffect_sd_prob) / prior_iideffect_sd_max;
  Type pcdensityiid = lambda / 2 * pow(iideffect_tau, -3/2) * exp( - lambda * pow(iideffect_tau, -1/2));
  nll -= log(pcdensityiid);
  
  // Build spde matrix
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  // Likelihood of the random field.
  nll += SCALE(GMRF(Q), 1.0 / tau)(nodemean);
  
  Type nllpriors = nll;

  // ------------------------------------------------------------------------ //
  // Calculate random field effects
  // ------------------------------------------------------------------------ //
  
  // Calculate field for pixel data
  vector<Type> logit_prevalence_field;
  logit_prevalence_field = Apixel * nodemean;

  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //

  vector<Type> pixel_linear_pred(n_points);
  vector<Type> pixel_pred(n_points);
  vector<Type> pixel_inc(n_points);
  vector<Type> reportnll(n_points);

  pixel_linear_pred = intercept + x * slope  + logit_prevalence_field.array();
  pixel_pred = invlogit(pixel_linear_pred);

  pixel_inc = pixel_pred * prev_inc_par[0] +
                    pixel_pred.pow(2) * prev_inc_par[1] +
                    pixel_pred.pow(3) * prev_inc_par[2];
  
  pixel_inc = prev_inc_extra_slope * pixel_inc;
  
  for (int prev_point = 0; prev_point < n_prev_points; prev_point++) {
    nll -= dbinom(positive_cases[prev_point], examined_cases[prev_point], pixel_pred[prev_point], true);
    reportnll[prev_point] = -dbinom(positive_cases[prev_point], examined_cases[prev_point], pixel_pred[prev_point], true);
  }

  int point_index;
  vector<Type> pixel_inc_count(n_inc_points);
  for (int inc_point = 0; inc_point < n_inc_points; inc_point++) {
    point_index = inc_point + n_prev_points;
    pixel_inc_count[inc_point] = pixel_inc[point_index] * incidence_pop[inc_point];
    nll -= dpois(incidence_cases[inc_point], pixel_inc_count[inc_point], true);
    reportnll[point_index] = -dpois(incidence_cases[inc_point], pixel_inc_count[inc_point], true);
  }

  REPORT(pixel_pred);
  REPORT(pixel_inc);
  REPORT(pixel_inc_count)
  REPORT(reportnll);
  REPORT(positive_cases);
  REPORT(examined_cases);
  REPORT(incidence_cases)
  REPORT(nllpriors);
  REPORT(nll);
  
  return nll;
  }
