#
# Predict malaria prevalence for the 3rd model.
#
# Mean and realisations predictions
# Produce mean and uncertainty maps
#

##############
# Mean prediction
##############

parameters <- obj$env$last.par.best
parameters <- split(parameters, names(parameters))

rasterForPoints <- cov_rasters[[1]]
rasterForPoints[is.na(rasterForPoints)] <- -9999
raster_pts <- rasterToPoints(rasterForPoints, spatial = TRUE)
coords <- raster_pts@coords

# Get random field predicted
spde <- (inla.spde2.matern(mesh, alpha = 2)$param.inla)[c("M0", "M1", "M2")]	
n_s <- nrow(spde$M0)						

Amatrix <- inla.mesh.project(mesh, loc = as.matrix(coords))$A

field <- (Amatrix %*% parameters$nodemean)[, 1]
field_ras <- rasterFromXYZ(cbind(coords, field))

covs_by_betas <- list()
for(i in seq_len(nlayers(cov_rasters))){
  covs_by_betas[[i]] <- parameters$slope[i] * cov_rasters[[i]]
}

cov_by_betas <- stack(covs_by_betas)

linear_predictor <- sum(cov_by_betas) + parameters$intercept + field_ras

mean_prevalence <- 1 / (1 + exp(-1 * linear_predictor))
plot(mean_prevalence)
plot(kenya_shapefile, add = T)


##############
# Uncertainty prediction
##############

N <- 100
CI <- 0.95

parameters <- obj$env$last.par.best

ch <- Cholesky(sd_out$jointPrecision)
par_draws <- sparseMVN::rmvn.sparse(N, parameters, ch, prec = TRUE)

prevalence <- list()

for(r in seq_len(N)){
  
  p <- split(par_draws[r, ], names(parameters))
  
  field <- (Amatrix %*% p$nodemean)[, 1]
  field_ras <- rasterFromXYZ(cbind(coords, field))
  
  covs_by_betas <- list()
  for(i in seq_len(nlayers(cov_rasters))){
    covs_by_betas[[i]] <- p$slope[i] * cov_rasters[[i]]
  }
  
  cov_by_betas <- stack(covs_by_betas)
  linear_predictor <- sum(cov_by_betas) + p$intercept + field_ras
  
  prevalence[[r]] <- 1 / (1 + exp(-1 * linear_predictor))
  
}

prevalence <- do.call(stack, prevalence)

probs <- c((1 - CI) / 2, 1 - (1 - CI) / 2)
prevalence_ci <- calc(prevalence, function(x) quantile(x, probs = probs, na.rm = TRUE))

plot(prevalence_ci[[2]] - prevalence_ci[[1]])
plot(kenya_shapefile, add = T)




