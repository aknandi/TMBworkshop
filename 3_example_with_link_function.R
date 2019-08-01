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

library(dplyr)
library(sp)
library(raster)
library(TMB)
library(stats)

response_data_path <- 'Z:/GBD2019/Processing/Stages/04b_PR_DB_Import_Export/Checkpoint_Outputs/pfpr.csv'

accessibility_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif' 
elevation_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/elev_srtm/Elev_5km.tif'
temperature_path <- 'Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'

cov_raster_paths <- c(accessibility_path, elevation_path, temperature_path)

# Get prevelance data
response <- read.csv(response_data_path)
response_reduced <- response[response$country == 'Zambia' & response$year_end == '2012', ] %>% 
  dplyr::select(c(latitude, longitude, examined, pf_pos, pf_pr))
response_reduced <- response_reduced[response_reduced$examined > 20, ]

survey_loc <- SpatialPoints(response_reduced[ , 2:1], 
                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Get covariate data
cov_rasters <- lapply(cov_raster_paths, raster)
cov_rasters <- lapply(cov_rasters, function(x) crop(x, extent(survey_loc)))
cov_rasters <- lapply(cov_rasters, function(x) scale(log(200 + x)))
cov_rasters <- stack(cov_rasters)

cov_df <- data.frame(raster::extract(cov_rasters, survey_loc))

# Remove any points with NA values in covariates
cov_matrix <- cov_df[complete.cases(cov_df), ] %>% as.matrix
response_reduced <- response_reduced[complete.cases(cov_df), ]

compile('src/model3.cpp')
dyn.load('src/model3')

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
                   prev_inc_par = c(2.616, -3.596, 1.594), 
                   positive_cases = response_reduced$pf_pos,
                   examined_cases = response_reduced$examined)

obj <- MakeADFun(
  data = input_data, 
  parameters = parameters,
  random = c('nodemean'),
  DLL = "model3")

its <- 10
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = its, eval.max = 2*its, trace = 0))

sd_out <- sdreport(obj, getJointPrecision = TRUE)

report <- obj$report()
plot(report$positive_cases/report$examined_cases, report$pixel_pred)
