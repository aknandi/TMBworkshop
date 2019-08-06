

response_data_path <- 'Z:/GBD2019/Processing/Stages/04b_PR_DB_Import_Export/Checkpoint_Outputs/pfpr.csv'

accessibility_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif' 
elevation_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/elev_srtm/Elev_5km.tif'
temperature_path <- 'Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'

cov_raster_paths <- c(accessibility_path, elevation_path, temperature_path)

# Get prevelance data
response <- read.csv(response_data_path)
response_reduced <- response[response$country == 'Madagascar' & response$year_end %in% 2013:2017, ] %>% 
  dplyr::select(c(latitude, longitude, examined, pf_pos, pf_pr))
response_reduced <- response_reduced[response_reduced$examined > 20, ]

survey_loc <- SpatialPointsDataFrame(response_reduced[ , 2:1], 
                                     data = response_reduced[, c("pf_pr","examined")],
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
