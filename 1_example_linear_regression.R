#
# Implements a linear regression model in TMB
#
# Inputs: malaria prevalence survey points
#         covariate rasters: accessibility, elevation, temperature suitability
#
# Anita Nandi
# 01/08/19
#

library(dplyr)
library(sp)
library(raster)
library(TMB)
library(stats)

compile('src/model1.cpp')
dyn.load('src/model1')

source('prepare_data.R')

head(response_reduced)
head(cov_matrix)

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_point_sd = -2.3)

input_data <- list(x = cov_matrix, 
                   positive_cases = response_reduced$pf_pos,
                   examined_cases = response_reduced$examined)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  DLL = "model1")

its <- 10
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))

sd_out <- sdreport(obj, getJointPrecision = TRUE)

report <- obj$report()
plot(report$positive_cases/report$examined_cases, report$pixel_pred)
