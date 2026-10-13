# Author: Persia Khan (persia.khan@gmail.com)
# Date Last Modified: 2026-04-22

################################################################################
# 2) SPECIES MATRICES - using the cumulative detection dataset from the HT Camera 
# Project (cleaned and reduced into independent detections during data cleaning 
# phase post-processing) and the project effort/camera site file. This R file
# creates matrices for each species analyzed in Khan et al. (2026) to represent
# detection/non-detection.
################################################################################

### GLMM PREP ### --------------------------------------------------------------

# Load packages
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)

setwd("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Detection Data")

# Read in datasets 
ind_det_data <- read_csv("HT_independent_detection_dataset_cumulative_2025-11-21.csv")
operational_periods <- read_csv("HT_operational_periods_cumulative_2026-02-17.csv")

## FORMAT DETECTION DATA -------------------------------------------------------

# Collapse sites with multiple site ID's 
ind_det_data <- ind_det_data %>%
  mutate(SiteID = dplyr::recode(as.character(SiteID),
                                "HT632A"  = "HT632",
                                "HT2300A" = "HT2300",
                                "HT391A"  = "HT391"))

operational_periods <- operational_periods %>%
  mutate(SiteID = dplyr::recode(as.character(SiteID),
                                "HT632A"  = "HT632",
                                "HT2300A" = "HT2300",
                                "HT391A"  = "HT391"))

## DEER ## ---------------------------------------------------------------------
# Filter for deer independent detections 
deer_det <- ind_det_data %>%
  filter(Species == "sitka black-tailed deer")

# Establish the start and end date of data collection period, then trim independent
# detections using this temporal window
start_date <- as.Date("2023-09-01")
end_date <- as.Date("2025-08-31")

# Trim detections to study period
deer_det_trimmed <- deer_det %>%
  filter(
    ymd_hms(DateTime) >= start_date &
      ymd_hms(DateTime) <= end_date
  )

# Create a data frame storing site–year–month deer detections
deer_monthly_counts <- deer_det_trimmed %>%
  mutate(
    Year  = factor(year(DateTime)),   
    Month = factor(month(DateTime))   
  ) %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Detections = n(), .groups = "drop")

# Create all site–year–month combinations across the study period
month_seq <- seq(
  floor_date(min(deer_det_trimmed$DateTime), "month"),
  floor_date(max(deer_det_trimmed$DateTime), "month"),
  by = "month"
)

all_site_month_combos <- expand.grid(
  SiteID = unique(deer_det_trimmed$SiteID),
  Year   = factor(year(month_seq)),
  Month  = factor(month(month_seq))
) %>%
  distinct()

# Join detections to full site–year–month grid and assign 0s
deer_monthly_counts_full <- all_site_month_combos %>%
  left_join(deer_monthly_counts,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Expand each operational period into daily dates
site_operations_expanded <- operational_periods %>%
  rowwise() %>%
  mutate(
    DailyDates = list(seq(ymd(OperationalStart),
                          ymd(OperationalEnd),
                          by = "day"))
  ) %>%
  unnest(DailyDates)

# Extract year and month from daily effort
site_operations_expanded <- site_operations_expanded %>%
  mutate(
    Year  = factor(year(DailyDates)),  
    Month = factor(month(DailyDates))   
  )

# Calculate monthly camera effort per site
site_month_effort <- site_operations_expanded %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Monthly_cam_effort = n(), .groups = "drop") %>% 
  filter(!(Year == 2023 & Month %in% c(6, 7, 8)))

# Combine effort and detections; keep only months with camera activity
deer_data <- site_month_effort %>%
  left_join(deer_monthly_counts_full,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Double check that all sites are included (n = 100)
unique_site_ids <- unique(deer_data$SiteID)
print(unique_site_ids) # Good to go!

# Double check all detections are preserved
sum(deer_data$Detections) #should equal the trimmed dataset detections

# Double check for duplicates in the Site/Year/Month combo
deer_data %>%
  count(SiteID, Year, Month) %>%
  filter(n > 1) #Should return 0x4

# Double check that there are no detections where effort = 0
deer_data %>%
  filter(Monthly_cam_effort == 0 & Detections > 0) #Should return 0x5

# Save new matrix
write.csv(deer_data, "deer_matrix-2026-04-22.csv")


## WOLF ## ---------------------------------------------------------------------
# Filter for wolf independent detections 
wolf_det <- ind_det_data %>%
  filter(Species == "wolf")

# Establish the start and end date of data collection period, then trim independent
# detections using this temporal window
start_date <- as.Date("2023-09-01")
end_date <- as.Date("2025-08-31")

# Trim detections to study period
wolf_det_trimmed <- wolf_det %>%
  filter(
    ymd_hms(DateTime) >= start_date &
      ymd_hms(DateTime) <= end_date
  )

# Create a data frame storing site–year–month wolf detections
wolf_monthly_counts <- wolf_det_trimmed %>%
  mutate(
    Year  = factor(year(DateTime)),   
    Month = factor(month(DateTime))   
  ) %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Detections = n(), .groups = "drop")

# Create all site–year–month combinations across the study period
month_seq <- seq(
  floor_date(min(wolf_det_trimmed$DateTime), "month"),
  floor_date(max(wolf_det_trimmed$DateTime), "month"),
  by = "month"
)

all_site_month_combos <- expand.grid(
  SiteID = unique(wolf_det_trimmed$SiteID),
  Year   = factor(year(month_seq)),
  Month  = factor(month(month_seq))
) %>%
  distinct()

# Join detections to full site–year–month grid and assign 0s
wolf_monthly_counts_full <- all_site_month_combos %>%
  left_join(wolf_monthly_counts,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Expand each operational period into daily dates
site_operations_expanded <- operational_periods %>%
  rowwise() %>%
  mutate(
    DailyDates = list(seq(ymd(OperationalStart),
                          ymd(OperationalEnd),
                          by = "day"))
  ) %>%
  unnest(DailyDates)

# Extract year and month from daily effort
site_operations_expanded <- site_operations_expanded %>%
  mutate(
    Year  = factor(year(DailyDates)), 
    Month = factor(month(DailyDates))   
  )

# Calculate monthly camera effort per site
site_month_effort <- site_operations_expanded %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Monthly_cam_effort = n(), .groups = "drop") %>% 
  filter(!(Year == 2023 & Month %in% c(6, 7, 8)))

# Combine effort and detections; keep only months with camera activity
wolf_data <- site_month_effort %>%
  left_join(wolf_monthly_counts_full,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Double check that all sites are included (n = 100)
unique_site_ids <- unique(wolf_data$SiteID)
print(unique_site_ids) # Good to go!

# Double check all detections are preserved
sum(wolf_data$Detections) #should equal the trimmed dataset detections

# Double check for duplicates in the Site/Year/Month combo
wolf_data %>%
  count(SiteID, Year, Month) %>%
  filter(n > 1) #Should return 0x4

# Double check that there are no detections where effort = 0
wolf_data %>%
  filter(Monthly_cam_effort == 0 & Detections > 0) #Should return 0x5

# Save new matrix
write.csv(wolf_data, "wolf_matrix-2026-04-22.csv")


## BLACK BEAR ## ---------------------------------------------------------------
# Filter for black bear independent detections 
bb_det <- ind_det_data %>%
  filter(Species == "black bear")

# Establish the start and end date of data collection period, then trim independent
# detections using this temporal window
start_date <- as.Date("2023-09-01")
end_date <- as.Date("2025-08-31")

# Trim detections to study period
bb_det_trimmed <- bb_det %>%
  filter(
    ymd_hms(DateTime) >= start_date &
      ymd_hms(DateTime) <= end_date
  )

# Create a data frame storing site–year–month bb detections
bb_monthly_counts <- bb_det_trimmed %>%
  mutate(
    Year  = factor(year(DateTime)),   
    Month = factor(month(DateTime))   
  ) %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Detections = n(), .groups = "drop")

# Create all site–year–month combinations across the study period
month_seq <- seq(
  floor_date(min(bb_det_trimmed$DateTime), "month"),
  floor_date(max(bb_det_trimmed$DateTime), "month"),
  by = "month"
)

all_site_month_combos <- expand.grid(
  SiteID = unique(bb_det_trimmed$SiteID),
  Year   = factor(year(month_seq)),
  Month  = factor(month(month_seq))
) %>%
  distinct()

# Join detections to full site–year–month grid and assign 0s
bb_monthly_counts_full <- all_site_month_combos %>%
  left_join(bb_monthly_counts,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Expand each operational period into daily dates
site_operations_expanded <- operational_periods %>%
  rowwise() %>%
  mutate(
    DailyDates = list(seq(ymd(OperationalStart),
                          ymd(OperationalEnd),
                          by = "day"))
  ) %>%
  unnest(DailyDates)

# Extract year and month from daily effort
site_operations_expanded <- site_operations_expanded %>%
  mutate(
    Year  = factor(year(DailyDates)),  
    Month = factor(month(DailyDates))  
  )

# Calculate monthly camera effort per site
site_month_effort <- site_operations_expanded %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Monthly_cam_effort = n(), .groups = "drop") %>% 
  filter(!(Year == 2023 & Month %in% c(6, 7, 8)))

# Combine effort and detections; keep only months with camera activity
bb_data <- site_month_effort %>%
  left_join(bb_monthly_counts_full,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Double check that all sites are included (n = 100)
unique_site_ids <- unique(bb_data$SiteID)
print(unique_site_ids) # Good to go!

# Double check all detections are preserved
sum(bb_data$Detections) #should equal the trimmed dataset detections

# Double check for duplicates in the Site/Year/Month combo
bb_data %>%
  count(SiteID, Year, Month) %>%
  filter(n > 1) #Should return 0x4

# Double check that there are no detections where effort = 0
bb_data %>%
  filter(Monthly_cam_effort == 0 & Detections > 0) #Should return 0x5

# Save new matrix
write.csv(bb_data, "bb_matrix-2026-04-22.csv")


## GRIZZLY BEAR ## -------------------------------------------------------------
# Filter for grizzly bear independent detections 
gb_det <- ind_det_data %>%
  filter(Species == "grizzly bear")

# Establish the start and end date of data collection period, then trim independent
# detections using this temporal window
start_date <- as.Date("2023-09-01")
end_date <- as.Date("2025-08-31")

# Trim detections to study period
gb_det_trimmed <- gb_det %>%
  filter(
    ymd_hms(DateTime) >= start_date &
      ymd_hms(DateTime) <= end_date
  )

# Create a data frame storing site–year–month gb detections
gb_monthly_counts <- gb_det_trimmed %>%
  mutate(
    Year  = factor(year(DateTime)),   
    Month = factor(month(DateTime))   
  ) %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Detections = n(), .groups = "drop")

# Create all site–year–month combinations across the study period
month_seq <- seq(
  floor_date(min(gb_det_trimmed$DateTime), "month"),
  floor_date(max(gb_det_trimmed$DateTime), "month"),
  by = "month"
)

all_site_month_combos <- expand.grid(
  SiteID = unique(gb_det_trimmed$SiteID),
  Year   = factor(year(month_seq)),
  Month  = factor(month(month_seq))
) %>%
  distinct()

# Join detections to full site–year–month grid and assign 0s
gb_monthly_counts_full <- all_site_month_combos %>%
  left_join(gb_monthly_counts,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Expand each operational period into daily dates
site_operations_expanded <- operational_periods %>%
  rowwise() %>%
  mutate(
    DailyDates = list(seq(ymd(OperationalStart),
                          ymd(OperationalEnd),
                          by = "day"))
  ) %>%
  unnest(DailyDates)

# Extract year and month from daily effort
site_operations_expanded <- site_operations_expanded %>%
  mutate(
    Year  = factor(year(DailyDates)),   # Extract year
    Month = factor(month(DailyDates))   # Extract month (1–12)
  )

# Calculate monthly camera effort per site
site_month_effort <- site_operations_expanded %>%
  group_by(SiteID, Year, Month) %>%
  summarise(Monthly_cam_effort = n(), .groups = "drop") %>% 
  filter(!(Year == 2023 & Month %in% c(6, 7, 8)))

# Combine effort and detections; keep only months with camera activity
gb_data <- site_month_effort %>%
  left_join(gb_monthly_counts_full,
            by = c("SiteID", "Year", "Month")) %>%
  mutate(Detections = replace_na(Detections, 0))

# Double check that all sites are included (n = 100)
unique_site_ids <- unique(gb_data$SiteID)
print(unique_site_ids) # Good to go!

# Double check all detections are preserved
sum(gb_data$Detections) #should equal the trimmed dataset detections

# Double check for duplicates in the Site/Year/Month combo
gb_data %>%
  count(SiteID, Year, Month) %>%
  filter(n > 1) #Should return 0x4

# Double check that there are no detections where effort = 0
gb_data %>%
  filter(Monthly_cam_effort == 0 & Detections > 0) #Should return 0x5

# Save new matrix
write.csv(gb_data, "gb_matrix-2026-04-22.csv")

