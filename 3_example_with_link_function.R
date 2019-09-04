#
# Implements a linear model with gaussian field and optimised link function in TMB
#
# Inputs: malaria prevalence survey points
#         malaria incidence survey points
#         covariate rasters: temperature, temperature seasonality, precipitation
#         gaussian field
#
# Anita Nandi
# 01/08/19
#

list.of.packages <- c("malariaAtlas", "dplyr", "sp", "raster", "TMB", "stats", "ggplot2", "Matrix")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(malariaAtlas)
library(dplyr)
library(sp)
library(raster)
library(TMB)
library(stats)
library(ggplot2)
library(Matrix)
library(INLA)

# Not totally necessary.
library(cowplot)

############
# Data preparation
############

# We now want to use incidence data

# source('prepare_data.R')
source('prepare_incidence_data.R')

head(inc_data)

coords <- as.matrix(rbind(prev_data[, c('longitude', 'latitude')], inc_data[, c('longitude', 'latitude')]))

mesh <- INLA::inla.mesh.2d(prev_data[ , 2:1],
                     cutoff = 0.1,
                     max.edge = c(0.5, 2), 
                     offset = c(1, 3))
spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
Apix <- INLA::inla.mesh.project(mesh, loc = coords)$A
n_s <- nrow(spde$M0)

#-------------- MAIN MODEL FITTING CODE ------------

compile('src/model3.cpp')
dyn.load(dynlib('src/model3'))


############
# Model fitting
############

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_prev_inc_extra_slope = 0,
                   log_kappa = -3,
                   log_tau = -0.5,
                   nodemean = rep(0, n_s))

input_data <- list(x = cov_matrix, 
                   Apixel = Apix,
                   spde = spde,
                   prev_inc_par = c(2.616, -3.596, 1.594),
                   priormean_log_prev_inc_extra_slope = 0,
                   priorsd_log_prev_inc_extra_slope = 0.001,
                   priormean_intercept = -4.0,
                   priorsd_intercept = 2.0,
                   priormean_slope = 0.0,
                   priorsd_slope = 0.5,
                   priormean_log_kappa = -3,
                   priorsd_log_kappa = 0.5,
                   priormean_log_tau = -0.50,
                   priorsd_log_tau = 2.0,
                   positive_cases = prev_data$positive,
                   examined_cases = prev_data$examined,
                   incidence_cases = inc_data$incidence,
                   incidence_pop = inc_data$population)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean'),
  DLL = "model3")

its <- 1000
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, trace = 0))

sd_out <- sdreport(obj, getJointPrecision = TRUE)
summary(sd_out, select = 'fixed')

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


# Mean and uncertainty maps
source('predict.R')
