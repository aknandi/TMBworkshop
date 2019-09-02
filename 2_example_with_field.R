#
# Implements a linear model with gaussian field in TMB
#
# Inputs: malaria prevalence survey points
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

############
# Data preparation
############

# source('prepare_data.R')

# Make coordinates for building the INLA mesh
x <- runif(nrow(prev_data), survey_loc_prev@bbox[1, 1], survey_loc_prev@bbox[1, 2])
y <- runif(nrow(prev_data), survey_loc_prev@bbox[2, 1], survey_loc_prev@bbox[2, 2])
coords <- cbind(x, y) %>% as.matrix()

# Build INLA mesh for the gaussian field
mesh <- INLA::inla.mesh.2d(prev_data[ , 2:1],
                     cutoff = 0.1,
                     max.edge = c(0.5, 2), 
                     offset = c(1, 3))

# Get spde and A matrix for the gaussian field
spde <- (INLA::inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
Apix <- INLA::inla.mesh.project(mesh, loc = coords)$A
n_s <- nrow(spde$M0)

#-------------- MAIN MODEL FITTING CODE ------------


compile('src/model2.cpp')
dyn.load(dynlib('src/model2'))

############
# Model fitting
############

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)),
                   log_kappa = 0,
                   log_tau = 1,
                   nodemean = rep(0, n_s))

input_data <- list(x = cov_matrix, 
                   Apixel = Apix,
                   spde = spde, 
                   positive_cases = prev_data$positive,
                   examined_cases = prev_data$examined)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean'),
  DLL = "model2")

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
pred_df <- data.frame(obs = report$positive_cases/report$examined_cases, pred = report$pixel_pred)
insample_plot <- ggplot(pred_df, aes(obs, pred)) + geom_point() + geom_abline(intercept = 0, slope = 1, color = 'red')

source('predict.R')






