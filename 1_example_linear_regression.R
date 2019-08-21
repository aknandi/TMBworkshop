#
# Implements a linear regression model in TMB
#
# Inputs: malaria prevalence survey points
#         covariate rasters: temperature, temperature seasonality, precipitation
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
library(ggplot2)
library(Matrix)

############
# Data preparation
############

source('prepare_data.R')

# We downloaded and organised some prevalence data for Kenya.
head(prev_data)

# We downloaded some rasters of environmental data.
cov_rasters

# And then organised it into a matrix matching the prevalence data.
head(cov_matrix)

plot(cov_matrix[, 'Precipitation'], prev_data$pr)


#-------------- MAIN MODEL FITTING CODE ------------

compile('src/model1.cpp')
dyn.load(dynlib('src/model1'))

############
# Model fitting
############

parameters <- list(intercept = -5,
                   slope = rep(0, ncol(cov_matrix)))

input_data <- list(x = cov_matrix, 
                   positive_cases = prev_data$positive,
                   examined_cases = prev_data$examined)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  DLL = "model1")

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

# Mean prediction

parameters <- obj$env$last.par.best
parameters <- split(parameters, names(parameters))

covs_by_betas <- list()
for(i in seq_len(nlayers(cov_rasters))){
  covs_by_betas[[i]] <- parameters$slope[i] * cov_rasters[[i]]
}

cov_by_betas <- stack(covs_by_betas)
linear_predictor <- sum(cov_by_betas) + parameters$intercept

mean_prevalence <- 1 / (1 + exp(-1 * linear_predictor))

plot(mean_prevalence)
plot(kenya_shapefile, add = T)





