#
# Prepare data for all three models
#
# Prevalence surveys and covariate data extracted at survey locations
#

# Get prevelance data
prev_data_full <- getPR(country = 'Kenya', species = "Pf")
kenya_shapefile <- getShp(ISO = "KEN", admin_level = c("admin0"))
autoplot(prev_data_full)

prev_data <- dplyr::select(prev_data_full, c(latitude, longitude, examined, positive, pr))
prev_data <- prev_data[complete.cases(prev_data), ]

survey_loc_prev <- SpatialPoints(as.data.frame(prev_data[ , 2:1]), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Get covariate data from World Clim (https://www.worldclim.org/bioclim)
r <- getData("worldclim", var = "bio", res = 5)
cov_rasters <- r[[c(1, 4, 12)]]
names(cov_rasters) <- c('Temperature', 'TemperatureSeasonality', 'Precipitation')
cov_rasters <- crop(cov_rasters, extent(kenya_shapefile))
cov_rasters <- scale(log(cov_rasters + 50))

cov_df <- data.frame(raster::extract(cov_rasters, survey_loc_prev))

# Remove any points with NA values in covariates
cov_matrix <- cov_df[complete.cases(cov_df), ] %>% as.matrix
prev_data <- prev_data[complete.cases(cov_df), ]


