#
# Implements a linear model with gaussian field and optimised link function in TMB
#
# Inputs: malaria prevalence survey points
#         malaria incidence survey points
#         covariate rasters: accessibility, elevation, temperature suitability
#         gaussian field
#
# Anita Nandi
# 01/08/19
#

library(malariaAtlas)
library(dplyr)
library(sp)
library(raster)
library(TMB)
library(stats)

compile('src/model3.cpp')
dyn.load('src/model3')

############
# Data preparation
############

# We now want to use incidence data
inc <- TRUE

source('prepare_data.R')

x <- runif(nrow(cov_matrix), survey_loc@bbox[1, 1], survey_loc@bbox[1, 2])
y <- runif(nrow(cov_matrix), survey_loc@bbox[2, 1], survey_loc@bbox[2, 2])
coords <- cbind(x, y) %>% as.matrix()

mesh <- INLA::inla.mesh.create(survey_loc_df, extend = list(offset = 4))
spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
Apix <- INLA::inla.mesh.project(mesh, loc = coords)$A
n_s <- nrow(spde$M0)

#-------------- MAIN MODEL FITTING CODE ------------

############
# Model fitting
############

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_point_sd = -2.3,
                   log_kappa = -3,
                   log_tau = -0.5,
                   nodemean = rep(0, n_s))

input_data <- list(x = cov_matrix, 
                   Apixel = Apix,
                   spde = spde,
                   prev_inc_par = c(2.616, -3.596, 1.594), 
                   positive_cases = prev_data$positive,
                   examined_cases = prev_data$examined,
                   incidence_cases = inc_data$incidence,
                   incidence_pop = inc_data$population)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean'),
  DLL = "model3")

its <- 10
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))

sd_out <- sdreport(obj, getJointPrecision = TRUE)

report <- obj$report()

#---------------------------------------------

########
## Model prediction
########

# In sample performance
pred_prev_df <- data.frame(obs = report$positive_cases/report$examined_cases, 
                           pred = report$pixel_pred[1:length(report$examined_cases)])
insample_prev_plot <- ggplot(pred_prev_df, aes(obs, pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red')

pred_inc_df <- data.frame(obs = report$incidence_cases, 
                          pred = report$pixel_inc_count)
insample_inc_plot <- ggplot(pred_inc_df, aes(obs, pred)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red')


print(cowplot::plot_grid(plotlist = list(insample_prev_plot, insample_inc_plot)))
