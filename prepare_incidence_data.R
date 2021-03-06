#
# Prepare incidence data for model 3
#

inc_data_path <- 
  'https://datadryad.org/bitstream/handle/10255/dryad.82194/PfPvAllData01042015_AgeStand.csv?sequence=1'

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


