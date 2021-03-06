TMB workshop
===========

Presentation to accompany the workshop can be found here: https://docs.google.com/presentation/d/1foFjojH5EN8MyzqNVSeTafMI8m1BKgfrSer0Bku9qwE/edit?usp=sharing

Helpful links for reference:
* Comprehensive TMB documentation: http://kaskr.github.io/adcomp/_book/Introduction.html
* Specific help for coming to C++ from R: http://kaskr.github.io/adcomp/_book/CppTutorial.html#i-know-r-but-not-c
* Github link with some simple TMB examples: https://github.com/kaskr/adcomp/wiki

Overview
-------

### Data used

Malaria prevalence data: Open-access malaria data hosted by the Malaria Atlas Project https://map.ox.ac.uk/ accessed using the malariaAtlas R package

Malaria incidence data: Open-access data from https://datadryad.org/bitstream/handle/10255/dryad.82194/PfPvAllData01042015_AgeStand.csv?sequence=1. Saved in the data folder

Covaraite data: WorldClim data (temperature, temparature variation and precipitation) at 5 minutes resolution from https://www.worldclim.org/bioclim using the getData R function


### Structure of the course

Three models of increasing complexity implemented in TMB to predict malaria prevalence. 

**Model 1:** Linear regression with covariates

**Model 2:** Linear model with Gaussian field

**Model 3:** Two types of response data. Two contributions to the likelihood. Customised link function. Linear model with Gaussian field and optimised link function


### File structure

Three main scripts for 3 different models:
* 1_example_linear_regression.R
* 2_example_with_field.R
* 3_example_with_link_function.R

model specification cpp files (in src folder):
* src/model1.cpp
* src/model2.cpp
* src/model3.cpp

Data preparation:
* **prepare_data.R:** prepares prevalence survey data and covariate data for all 3 models
* **prepare_incidence_data.R:** prepared incidence survey data for 3rd model

Predictions:
* **predict.R:** predict prevalence for model 2 and 3. Plot mean and uncertainty


### Needed libraries

If you are coming to the workshop please install these packages ahead of time.

* malariaAtlas
* dplyr
* sp
* raster
* TMB
* INLA
* Matrix
* ggplot2
* cowplot

INLA must be installed from a separate repository.
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)





