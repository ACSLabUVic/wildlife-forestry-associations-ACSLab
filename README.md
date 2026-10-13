## Manuscript Title: Wildlife associations at the watershed-scale: Understanding species-habitat relationships for forestry planning in a temperate rainforest
### Description of the data and code files

The following datasets and code have been scrubbed of location data. 

1_covariate_preparation_Khan2026.R: COVARIATE PREPARATION - using shapefiles created in "spatial_data_cleanup.R", upload shapefiles (if not done in last step). The following script will assign covariate values to the site level and watershed level covariate dataframes that were prepared for the GLMM phase of analysis in Khan et al. (2026). 

2_species_matrix_Khan2026.R: SPECIES MATRICES - using the cumulative detection dataset from the HT Camera Project (cleaned and reduced into independent detections during data cleaning phase post-processing) and the project effort/camera site file. This R file creates matrices for each species analyzed in Khan et al. (2026) to represent detection/non-detection.

3_modelling_Khan2026.R: MODELLING - This script contains the candidate model set for each of the species of interest (wolf, sitka black-tailed deer, black bear, and grizzly bear) in Khan et al. (2026), as well as model selection and diagnostics. This script builds off of two previous # scripts: 1) '1_covariate_preparation_Khan2026.R', and 2) '2_species_matrix_Khan2026.R'

4_plotting_Khan2026.R: PLOTTING - The following code was used to produce the figures in Khan et al. (2026). This script builds on the model outputs generated in '3_modelling_Khan2026'. Figures are listed in the order they appear within the manuscript.

covariates_2026-04-22.csv: Output of covariate data attached to sites and watersheds generated in '1_covariate_preparation_Khan2026.R' and required to run '3_modelling_Khan2026.R'.

deer_matrix-2026-04-22.csv: Detection/non-detection matrix for Sitka black-tailed deer used in '3_modelling_Khan2026.R'. 

wolf_matrix-2026-04-22.csv: Detection/non-detection matrix for wolves used in '3_modelling_Khan2026.R'. 

bb_matrix-2026-04-22.csv: Detection/non-detection matrix for black bears used in '3_modelling_Khan2026.R'. 

gb_matrix-2026-04-22.csv: Detection/non-detection matrix for grizzly bears used in '3_modelling_Khan2026.R'. 

### Code/Software
Analyses were conducted in RStudio using version 2024.12.0. 
