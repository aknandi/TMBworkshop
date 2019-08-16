#
# Prepare data for all three models
#
# Prevalence surveys, incidence surveys (if required) and covariate data extracted at survey locations
#

accessibility_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/accessibility/accessibility.5k.MEAN.tif' 
elevation_path <- 'Z:/GBD2017/Processing/Static_Covariates/MAP/other_rasters/elev_srtm/Elev_5km.tif'
temperature_path <- 'Z:/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Synoptic/TSI-Martens2-Pf.Synoptic.Overall.Mean.5km.Data.tif'

cov_raster_paths <- c(accessibility_path, elevation_path, temperature_path)

# Get prevelance data
prev_data <- getPR(country = 'Kenya', species = "Pf")
kenya_shapefile <- getShp(ISO = "KEN", admin_level = c("admin0"))
autoplot(prev_data)

prev_data <- dplyr::select(prev_data, c(latitude, longitude, examined, positive, pr))
prev_data <- prev_data[complete.cases(prev_data), ]

survey_loc_prev <- SpatialPoints(as.data.frame(prev_data[ , 2:1]), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# Get covariate data
cov_rasters <- lapply(cov_raster_paths, raster)
cov_rasters <- lapply(cov_rasters, function(x) crop(x, extent(kenya_shapefile)))
cov_rasters <- lapply(cov_rasters, function(x) scale(log(200 + x)))
cov_rasters <- stack(cov_rasters)

cov_df <- data.frame(raster::extract(cov_rasters, survey_loc_prev))

# Remove any points with NA values in covariates
cov_matrix <- cov_df[complete.cases(cov_df), ] %>% as.matrix
prev_data <- prev_data[complete.cases(cov_df), ]

if(inc) {
  inc_data_path <- 'data/PfPvAllData01042015_AgeStand.csv'
  
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

