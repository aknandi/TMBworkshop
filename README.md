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
