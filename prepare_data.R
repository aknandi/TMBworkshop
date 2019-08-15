
accessibility_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif' 
elevation_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/elev_srtm/Elev_5km.tif'
temperature_path <- 'Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'

cov_raster_paths <- c(accessibility_path, elevation_path, temperature_path)

# Get prevelance data
response <- getPR(country = "Madagascar", species = "Pf")
autoplot(response)

response_reduced <- dplyr::select(response, c(latitude, longitude, examined, positive, pr))
response_reduced <- response_reduced[complete.cases(response_reduced), ]

survey_loc <- SpatialPoints(as.data.frame(response_reduced[ , 2:1]), 
                            proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Get covariate data
cov_rasters <- lapply(cov_raster_paths, raster)
cov_rasters <- lapply(cov_rasters, function(x) crop(x, extent(survey_loc) + 2))
cov_rasters <- lapply(cov_rasters, function(x) scale(log(200 + x)))
cov_rasters <- stack(cov_rasters)

cov_df <- data.frame(raster::extract(cov_rasters, survey_loc))

# Remove any points with NA values in covariates
cov_matrix <- cov_df[complete.cases(cov_df), ] %>% as.matrix
response_reduced <- response_reduced[complete.cases(cov_df), ]
