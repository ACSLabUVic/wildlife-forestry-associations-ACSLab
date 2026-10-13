# Author: Persia Khan (persia.khan@gmail.com)
# Date Last Modified: 2025-04-22

################################################################################
# 1) COVARIATE PREPARATION - using shapefiles created in "spatial_data_cleanup.R", 
# upload shapefiles (if not done in last step). The following script will assign
# covariate values to the site level and watershed level covariate dataframes that 
# were prepared for the GLMM phase of analysis in Khan et al. (2026). 
################################################################################

# **All spatial data is projected in BC Albers** 

library(readr)
library(tidyverse)
library(sf)
library(dplyr)
library(readxl)
library(vegan)
library(ggplot2)
library(units)
library(tidyr)
library(Hmisc)
library(corrplot)

### FILE UPLOAD (if not in environment) ### ------------------------------------
setwd("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Spatial Data")

camsites <- read_excel("HT_camera_project_site_database.xlsx")

sites <- st_read("HTsites_cleaned.shp")

buffer500 <- st_read("camsites_buff_500m.shp")

wsheds <- st_read("all_wsheds_FID.shp")

forest_age <- st_read("forest_age_master.shp")

coastline <- st_read("coastline_FALB_diss.shp")

salmon_biomass <- read_csv("focal_salm_biomass_per_year_w_density.csv")

################################################################################

options(scipen = 999)

### CREATE OBJECT FOR STORING COVARIATES ### -----------------------------------
# This object will be used to store all covariates used in the analysis
covariates <- buffer500 %>%
  select(SiteID, wshed_name, wshed_id, buff_area) %>%
  st_drop_geometry()

# Add watershed FID for later (to merge salmon data)
covariates <- covariates %>%
  left_join(wsheds %>% select(wshed_id, WTRSHD_FID), by = "wshed_id") %>%
  st_drop_geometry() %>%
  select(-geometry)

### FOREST AGE -> CLASSES ### --------------------------------------------------
# This section will classify stand age data into the following classes as proportions
# and categories: 
#   - Recent Harvest = 0-20yrs
#   - Second Growth = 21-80yrs
#   - Old Growth = 140+ yrs
# This is based on the VRI and Harvest Cutblock layers from the BC Gov, which have 
# been merged during data cleaning (using Rank 1 layer in VRI and prioritizing age
# in Harvest Cutblock as this layer was most accurate for harvest when compared 
# with imagery).

# Clean up object column names 
forest_age <- forest_age %>%
  rename(site_index = SITE_INDEX) %>%
  rename(species_code = SPECIES_CD) 

# Select columns of interest (this is a large dataset, needs to be reduced)
forest_age_select <- forest_age %>%
  select(stand_age, site_index, species_code, geometry)

# Intersect forest data with the stand-scale buffer (500m)
buff500_with_forest <- st_intersection(buffer500, forest_age_select)

# Calculate the area for each polygon in the intersection (there will be multiple
# polygons associated with each Site ID)
buff500_with_forest$stand_area <- as.numeric(st_area(buff500_with_forest))

# Classify using the criteria above 
buff500_with_forest <- buff500_with_forest %>%
  mutate(seral_class = case_when(
    stand_age >= 0 & stand_age <= 20 ~ "recent_harvest", 
    stand_age > 20 & stand_age <= 80 ~ "second_growth",
    stand_age > 80 & stand_age <= 140 ~ "maturing_forest",
    stand_age > 140  ~ "old_growth",
    TRUE ~ NA_character_
  ))

# Group by SiteID and seral_class 
class_area <- buff500_with_forest %>%
  group_by(SiteID, seral_class) %>%
  summarise(area = sum(stand_area, na.rm = TRUE), .groups = "drop") %>%
  st_drop_geometry()

# Pivot wider to create columns for each seral class and a row per site
class_area_summary <- class_area %>%
  pivot_wider(
    names_from = seral_class,                  
    values_from = area,      
    values_fill = list(area = 0),  
    names_glue = "{seral_class}_area"    
  )

# Attach buffer areas for calculations
class_area_summary <- class_area_summary %>%
  left_join(
    covariates %>% select(SiteID, buff_area), 
    by = "SiteID"
  )

# Calculate proportions of class_area / buff_area 
class_area_summary <- class_area_summary %>%
  mutate(
    # Proportion within buffer area
    recent_harvest_buff500 = recent_harvest_area / buff_area,
    second_growth_buff500 = second_growth_area / buff_area,
    maturing_forest_buff500 = maturing_forest_area / buff_area,
    old_growth_buff500 = old_growth_area / buff_area
  )

# Re-attach to covariates object
covariates <- covariates %>%
  left_join(
    class_area_summary %>% select(-buff_area), 
    by = "SiteID"
  )

### DISTANCE TO COASTLINE ### --------------------------------------------------

# Calculate the euclidean distance from the site to the coastline
dist_to_coast <- sites %>%
  mutate(coast_dist_m = st_distance(sites, coastline, which = "Euclidean")) %>%
  select(SiteID, coast_dist_m) %>%
  st_drop_geometry()

dist_to_coast$coast_dist_m <- as.numeric(dist_to_coast$coast_dist_m)

### ELEVATION ### --------------------------------------------------------------
# Add field collected elevation to each site (in meters)
elevation <- camsites %>%
  select(SiteID, elevation) %>%
  rename(elevation_m = elevation)

### MERGE COVARIATES ### -------------------------------------------------------
# Join distance to covariates to the main object
covariates <- covariates %>%
  left_join(dist_to_coast, by = "SiteID") %>%
  left_join(elevation, by = "SiteID") 

### WATERSHED FOREST COVARIATES ### --------------------------------------------
# Calculate area of each watershed
wsheds$wshed_area <- as.numeric(st_area(wsheds))

# Perform spatial intersection between wsheds and forest data
wshed_with_forest <- st_intersection(wsheds, forest_age_select)

# Calculate the area of each polygon in the selected watersheds
wshed_with_forest$stand_area <- as.numeric(st_area(wshed_with_forest))

# Classify using the criteria above 
wshed_with_forest <- wshed_with_forest %>%
  mutate(seral_class = case_when(
    stand_age >= 0 & stand_age <= 20 ~ "recent_harvest", 
    stand_age > 20 & stand_age <= 80 ~ "second_growth",
    stand_age > 80 & stand_age <= 140 ~ "maturing_forest",
    stand_age > 140  ~ "old_growth",
    TRUE ~ NA_character_
  ))

# Group by wshed_name and seral_class 
class_area_wshed <- wshed_with_forest %>%
  group_by(wshed_name, seral_class) %>%
  summarise(area = sum(stand_area, na.rm = TRUE), .groups = "drop") %>%
  st_drop_geometry()

# Pivot wider to create columns for each seral class and a row per site
class_area_summary_wshed <- class_area_wshed %>%
  pivot_wider(
    names_from = seral_class,                  
    values_from = area,      
    values_fill = list(area = 0),  
    names_glue = "{seral_class}_area_wshed"    
  ) 

class_area_summary_wshed <- class_area_summary_wshed %>%
  left_join(
    wsheds %>% select(wshed_name, wshed_area), 
    by = "wshed_name"
  ) %>%
  st_drop_geometry()

# Calculate proportions of class_area / buff_area and class_area / FALB_area
class_area_summary_wshed <- class_area_summary_wshed %>%
  mutate(
    # Proportion within wshed area
    recent_harvest_wshed = recent_harvest_area_wshed / wshed_area,
    second_growth_wshed = second_growth_area_wshed / wshed_area,
    maturing_forest_wshed = maturing_forest_area_wshed / wshed_area,
    old_growth_wshed = old_growth_area_wshed / wshed_area,
 )

# Attach to covariates object
covariates <- covariates %>%
  left_join(class_area_summary_wshed, by = "wshed_name")

### SALMON BY WATERSHED ### ----------------------------------------------------
# Bring in salmon data to calculate decadal biomass average for each watershed
salmon_biomass_avg <- salmon_biomass %>%
  group_by(WTRSHD_FID) %>%
  summarise(avg_sal_biomass = mean(LogBMDensityinFWS, na.rm = TRUE))

# Attach to covariates object
covariates <- covariates %>%
  left_join(salmon_biomass_avg, by = "WTRSHD_FID")

### FINAL CLEAN ### ------------------------------------------------------------
# NA -> 0
covariates[is.na(covariates)] <- 0

# drop any units
drop_units(covariates)

# Change site names to upper case to match detection data
covariates <- covariates %>%
  mutate(SiteID = toupper(SiteID)) %>%
  select(-geometry)

# write csv
write_excel_csv(covariates, 
                "C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Covariate Data/covariates_2026-04-22.csv")

### CORRELATION ### ------------------------------------------------------------
# Select for numeric covariates
numeric_covs <- covariates %>%
  select(where(is.numeric)) %>%
  select(where(~ sum(!is.na(.)) > 0))

numeric_covs_filtered <- numeric_covs %>%
  select(c("recent_harvest_buff500",
           "second_growth_buff500",
           "old_growth_buff500",
           "recent_harvest_wshed",
           "second_growth_wshed",
           "old_growth_wshed",
           "coast_dist_m",
           "elevation_m",
           "avg_sal_biomass"))

rcorr_result <- rcorr(
  as.matrix(numeric_covs_filtered),
  type = "pearson"
)

corrplot(
  rcorr_result$r,
  method = "number",
  type = "upper",
  tl.cex = 0.8,
  number.cex = 0.8,
  diag = FALSE
)
