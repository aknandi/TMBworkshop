#
# Prepare data for all three models
#
# Prevalence surveys, incidence surveys (if required) and covariate data extracted at survey locations
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

if(inc) {
  # Can be accessed from https://datadryad.org/bitstream/handle/10255/dryad.82194/PfPvAllData01042015_AgeStand.csv?sequence=1

  inc_data_path <- 'https://datadryad.org/bitstream/handle/10255/dryad.82194/PfPvAllData01042015_AgeStand.csv?sequence=1'
  #if(!file.exists(inc_data_path)) stop('Please download data from https://datadryad.org/bitstream/handle/10255/dryad.82194/PfPvAllData01042015_AgeStand.csv?sequence=1')  

  inc_data <- read.csv(inc_data_path)
  inc_data <- inc_data[inc_data$COUNTRY == 'Kenya', ]
  inc_data <- dplyr::select(inc_data, c(LAT, LONG, POP, INC))
  names(inc_data) <- c('latitude', 'longitude', 'population', 'incidence')
  inc_data <- inc_data[complete.cases(inc_data), ]
  inc_data$incidence_rate <- inc_data$incidence / inc_data$population
  
  
  survey_loc_inc <- SpatialPoints(inc_data[ , 2:1], 
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  cov_df_inc <- data.frame(raster::extract(cov_rasters, survey_loc_inc))
  
  # Remove any points with NA values in covariates
  cov_matrix_inc <- cov_df_inc[complete.cases(cov_df_inc), ] %>% as.matrix
  cov_matrix <- rbind(cov_matrix, cov_matrix_inc)
  inc_data <- inc_data[complete.cases(cov_df_inc), ]
  
  # Get combined survey locations for the mesh
  survey_loc_df <- rbind(as.data.frame(prev_data[ , 2:1]), inc_data[ , 2:1])
  survey_loc <- SpatialPoints(survey_loc_df, 
                              proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
}

