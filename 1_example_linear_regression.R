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

compile('src/model1.cpp')
dyn.load('src/model1')

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
