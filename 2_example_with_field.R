#
# Implements a linear model with gaussian field in TMB
#
# Inputs: malaria prevalence survey points
#         covariate rasters: accessibility, elevation, temperature suitability
#         gaussian field
#
# Anita Nandi
# 01/08/19
#

library(dplyr)
library(sp)
library(raster)
library(TMB)
library(stats)

compile('src/model2.cpp')
dyn.load('src/model2')

source('prepare_data.R')

x <- runif(nrow(response_reduced), survey_loc@bbox[1, 1], survey_loc@bbox[1, 2])
y <- runif(nrow(response_reduced), survey_loc@bbox[2, 1], survey_loc@bbox[2, 2])
coords <- cbind(x, y) %>% as.matrix()

mesh <- INLA::inla.mesh.create(response_reduced[ , 2:1], extend = list(offset = 4))
spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
Apix <- INLA::inla.mesh.project(mesh, loc = coords)$A
n_s <- nrow(spde$M0)

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_point_sd = -2.3,
                   log_kappa = -3,
                   log_tau = -0.5,
                   nodemean = rep(0, n_s))

input_data <- list(x = cov_matrix, 
                   Apixel = Apix,
                   spde = spde, 
                   positive_cases = response_reduced$pf_pos,
                   examined_cases = response_reduced$examined)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean'),
  DLL = "model2")

its <- 10
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))

sd_out <- sdreport(obj, getJointPrecision = TRUE)

report <- obj$report()

# In sample performance
pred_df <- data.frame(obs = report$positive_cases/report$examined_cases, pred = report$pixel_pred)
insample_plot <- ggplot(pred_df, aes(obs, pred)) + geom_point() + geom_abline(intercept = 0, slope = 1, color = 'red')
