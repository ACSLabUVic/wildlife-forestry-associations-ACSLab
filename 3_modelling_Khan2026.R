# Author: Persia Khan (persia.khan@gmail.com)
# Date Last Modified: 2026-04-22

################################################################################
# 3) MODELLING - This script contains the candidate model set for each of the 
# species of interest (wolf, sitka black-tailed deer, black bear, and grizzly bear) 
# in Khan et al. (2026), as well as model selection and diagnostics. This script 
# builds off of two previous # scripts: 1) '1_covariate_preparation_Khan2026.R', 
# and 2) '2_species_matrix_Khan2026.R'
################################################################################

# Load packages
library(car)
library(effects)
library(MuMIn)
library(DHARMa)
library(broom)
library(broom.mixed)
library(dotwhisker)
library(hms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(glmmTMB)
library(bbmle) 
library(geosphere)
library(performance)
library(formattable)
library(ggeffects)
library(ggformula)
library(gridExtra)
library(cowplot)
library(finalfit)
library(readr)
library(pROC)
library(interactions)
library(sf)
library(emmeans)
library(AICcmodavg)
library(ggpubr)
library(patchwork)

## COVARIATES ## ---------------------------------------------------------------
# Read in covariate data 
all_covs <- read_csv("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Covariate Data/covariates_2026-04-22.csv")

# NA -> 0
all_covs[is.na(all_covs)] <- 0

# Retain mean and standard deviations for each numeric covariate, to be back-transformed 
# later for predictions
scaling_info <- all_covs %>%
  summarise(across(where(is.numeric), list(mean = \(x) mean(x, na.rm = TRUE), 
                                           sd = \(x) sd(x, na.rm = TRUE)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("variable", "stat"),  
    names_pattern = "(.*)_(mean|sd)"  
  ) %>%
  pivot_wider(names_from = stat, values_from = value)

# Center and scale numeric covariates so that the mean value is 0
all_covs <- all_covs %>%
  mutate(across(where(~ is.numeric(.) & length(unique(.)) > 2), scale)) 

summary(all_covs) # check that mean values = ~0

## MERGE SPECIES AND COVARIATE DATA ## -----------------------------------------

## Wolves ## -------------------------------------------------------------------
wolf_data <- read.csv("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Detection Data/wolf_matrix-2026-04-22.csv")

# Create new presence/absence column 
wolf_data$Detections_binom <- ifelse(wolf_data$Detections > 0, 1, 0)

# Check distribution 
table(wolf_data$Detections_binom)

# Then join species matrix to covariate data
wolf_data <- wolf_data %>%
  left_join(all_covs, by = "SiteID") %>%
  dplyr::select(-ends_with("X"))

# Using timestamps, assign seasons based on the following: 
# - Summer = April through August 
# - Spawning Season = September through November 
# - Winter = December through March
wolf_data <- wolf_data %>%
  mutate(
    Season = case_when(
      lubridate::month(Month) %in% c(4, 5, 6, 7, 8) ~ "Summer (Apr-Aug)",
      lubridate::month(Month) %in% c(9, 10, 11) ~ "Spawning Season (Sept-Nov)",
      lubridate::month(Month) %in% c(12, 1, 2, 3) ~ "Winter (Dec-Mar)")) 

as.factor(wolf_data$Season)
as.factor(wolf_data$Year)
as.factor(wolf_data$Month)

# Convert effort measurement to log-scale for offset
wolf_data <- wolf_data %>%
  mutate(
    monthly_effort_log = log(Monthly_cam_effort),
    monthly_effort_log_sc = as.numeric(scale(monthly_effort_log))
  ) 

# Re-level Season factor if relevant
wolf_data$Season <- relevel(factor(wolf_data$Season), ref = "Summer (Apr-Aug)")

## Deer ## ---------------------------------------------------------------------
deer_data <- read.csv("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Detection Data/deer_matrix-2026-04-22.csv")

# Create new presence/absence column 
deer_data$Detections_binom <- ifelse(deer_data$Detections > 0, 1, 0)

# Check distribution 
table(deer_data$Detections_binom)

# Then join species matrix to covariate data
deer_data <- deer_data %>%
  left_join(all_covs, by = "SiteID") %>%
  dplyr::select(-ends_with("X"))

# Using timestamps, assign seasons based on the following: 
# - Summer = April through August 
# - Spawning Season = September through November 
# - Winter = December through March
deer_data <- deer_data %>%
  mutate(
    Season = case_when(
      lubridate::month(Month) %in% c(4, 5, 6, 7, 8) ~ "Summer (Apr-Aug)",
      lubridate::month(Month) %in% c(9, 10, 11) ~ "Spawning Season (Sept-Nov)",
      lubridate::month(Month) %in% c(12, 1, 2, 3) ~ "Winter (Dec-Mar)")) 

as.factor(deer_data$Season)
as.factor(deer_data$Year)
as.factor(deer_data$Month)

# Convert effort measurement to log-scale for offset
deer_data <- deer_data %>%
  mutate(
    monthly_effort_log = log(Monthly_cam_effort),
    monthly_effort_log_sc = as.numeric(scale(monthly_effort_log))
  )

# Re-level Season factor if relevant
deer_data$Season <- relevel(factor(deer_data$Season), ref = "Summer (Apr-Aug)")

## Black Bear ## ---------------------------------------------------------------
bb_data <- read.csv("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Detection Data/bb_matrix-2026-04-22.csv")

# Filter black bear data to remove hibernation period 
bb_data <- bb_data %>%
  filter(!is.na(Month)) %>%  # Remove NA values
  filter(!(month(Month) %in% c(12, 1, 2, 3)))

# Create new presence/absence column 
bb_data$Detections_binom <- ifelse(bb_data$Detections > 0, 1, 0)

# Check distribution 
table(bb_data$Detections_binom)

# Then join species matrix to covariate data
bb_data <- bb_data %>%
  left_join(all_covs, by = "SiteID") %>%
  dplyr::select(-ends_with("X"))

# Using timestamps, assign seasons based on the following: 
# - Summer = April through August 
# - Spawning Season = September through November 
# - Winter = December through March
bb_data <- bb_data %>%
  mutate(
    Season = case_when(
      lubridate::month(Month) %in% c(4, 5, 6, 7, 8) ~ "Summer (Apr-Aug)",
      lubridate::month(Month) %in% c(9, 10, 11) ~ "Spawning Season (Sept-Nov)",
      lubridate::month(Month) %in% c(12, 1, 2, 3) ~ "Winter (Dec-Mar)")) 

as.factor(bb_data$Season)
as.factor(bb_data$Year)
as.factor(bb_data$Month)

# Convert effort measurement to log-scale for offset
bb_data <- bb_data %>%
  mutate(
    monthly_effort_log = log(Monthly_cam_effort),
    monthly_effort_log_sc = as.numeric(scale(monthly_effort_log))
  )

# Re-level Season factor if relevant
bb_data$Season <- relevel(factor(bb_data$Season), ref = "Summer (Apr-Aug)")

## Grizzly Bear ## -------------------------------------------------------------
gb_data <- read.csv("C:/Users/persi/OneDrive - University of Victoria/MANUSCRIPT 1/Data/Detection Data/gb_matrix-2026-04-22.csv")

# Filter black bear data to remove hibernation period 
gb_data <- gb_data %>%
  filter(!is.na(Month)) %>%  # Remove NA values
  filter(!(month(Month) %in% c(12, 1, 2, 3)))

# Create new presence/absence column 
gb_data$Detections_binom <- ifelse(gb_data$Detections > 0, 1, 0)

# Check distribution 
table(gb_data$Detections_binom)

# Then join species matrix to covariate data
gb_data <- gb_data %>%
  left_join(all_covs, by = "SiteID") %>%
  dplyr::select(-ends_with("X"))

# Using timestamps, assign seasons based on the following: 
# - Summer = April through August 
# - Spawning Season = September through November 
# - Winter = December through March
gb_data <- gb_data %>%
  mutate(
    Season = case_when(
      lubridate::month(Month) %in% c(4, 5, 6, 7, 8) ~ "Summer (Apr-Aug)",
      lubridate::month(Month) %in% c(9, 10, 11) ~ "Spawning Season (Sept-Nov)",
      lubridate::month(Month) %in% c(12, 1, 2, 3) ~ "Winter (Dec-Mar)")) 

as.factor(gb_data$Season)
as.factor(gb_data$Year)
as.factor(gb_data$Month)

# Convert effort measurement to log-scale for offset
gb_data <- gb_data %>%
  mutate(
    monthly_effort_log = log(Monthly_cam_effort),
    monthly_effort_log_sc = as.numeric(scale(monthly_effort_log))
  )

# Re-level Season factor if relevant
gb_data$Season <- relevel(factor(gb_data$Season), ref = "Summer (Apr-Aug)")

## MODEL SETS BY SPECIES ## ----------------------------------------------------

## Wolf Candidate Model Set ## -------------------------------------------------
# null model (1)
W1 <- glmmTMB(Detections_binom ~ 1 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = wolf_data,
              ziformula = ~ 0,
              family = binomial)

# environmental model (2)
W2 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = wolf_data,
              ziformula = ~ 0,
              family = binomial)

# forest (stand) + environmental model (3)
W3 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m +
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_buff500 +
                second_growth_buff500 +
                old_growth_buff500 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = wolf_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) + environmental model (4)
W4 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = wolf_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) w/ interactions + environmental model (5)
W5 <- glmmTMB(Detections_binom ~
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                recent_harvest_wshed:Season +
                second_growth_wshed:Season +
                old_growth_wshed:Season +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = wolf_data,
              ziformula = ~ 0,
              family = binomial)


# Build AIC table
wolf.cand.models <- list("null model (1)" = W1,
                         "environmental model (2)" = W2, 
                         "forest (stand) + environmental model (3)" = W3,
                         "forest (watershed) + environmental model (4)" = W4,
                         "forest (watershed) w/ interactions + environmental model (5)" = W5)
wolf_selectionTable <- aictab(cand.set = wolf.cand.models, second.ord = FALSE)


## Deer Candidate Model Set ## -------------------------------------------------
# null model (1)
D1 <- glmmTMB(Detections_binom ~ 1 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

# environmental model (2)
D2 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

# forest (stand) + environmental model (3)
D3 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m +
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_buff500 +
                second_growth_buff500 +
                old_growth_buff500 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) + environmental model (4)
D4 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) w/ interactions + environmental model (5)
D5 <- glmmTMB(Detections_binom ~
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                recent_harvest_wshed:Season +
                second_growth_wshed:Season +
                old_growth_wshed:Season +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

#forest (stand) w/ interactions + environmental model (6)
D6 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_buff500 +
                second_growth_buff500 +
                old_growth_buff500 +
                recent_harvest_buff500:Season +
                second_growth_buff500:Season +
                old_growth_buff500:Season +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = deer_data,
              ziformula = ~ 0,
              family = binomial)

# Build AIC table
deer.cand.models <- list("null model (1)" = D1,
                         "environmental model (2)" = D2, 
                         "forest (stand) + environmental model (3)" = D3,
                         "forest (watershed) + environmental model (4)" = D4,
                         "forest (watershed) w/ interactions + environmental model (5)" = D5,
                         "forest (stand) w/ interactions + environmental model (6)" = D6)
deer_selectionTable <- aictab(cand.set = deer.cand.models, second.ord = FALSE)


## Black Bear Candidate Model Set ## -------------------------------------------
# null model (1)
B1 <- glmmTMB(Detections_binom ~ 1 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = bb_data,
              ziformula = ~ 0,
              family = binomial)

# environmental model (2)
B2 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = bb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (stand) + environmental model (3)
B3 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m +
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_buff500 +
                second_growth_buff500 +
                old_growth_buff500 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = bb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) + environmental model (4)
B4 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = bb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) w/ interactions + environmental model (5)
B5 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                recent_harvest_wshed:Season +
                second_growth_wshed:Season +
                old_growth_wshed:Season +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = bb_data,
              ziformula = ~ 0,
              family = binomial)


# Build AIC table
bb.cand.models <- list("null model (1)" = B1,
                       "environmental model (2)" = B2, 
                       "forest (stand) + environmental model (3)" = B3,
                       "forest (watershed) + environmental model (4)" = B4,
                       "forest (watershed) w/ interactions + environmental model (5)" = B5)
bb_selectionTable <- aictab(cand.set = bb.cand.models, second.ord = FALSE)

## Grizzly Bear Candidate Model Set ## -----------------------------------------
# null model (1)
G1 <- glmmTMB(Detections_binom ~ 1 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = gb_data,
              ziformula = ~ 0,
              family = binomial)

# environmental model (2)
G2 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = gb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (stand) + environmental model (3)
G3 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m +
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_buff500 +
                second_growth_buff500 +
                old_growth_buff500 +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = gb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) + environmental model (4)
G4 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = gb_data,
              ziformula = ~ 0,
              family = binomial)

# forest (watershed) w/ interactions + environmental model (5)
G5 <- glmmTMB(Detections_binom ~ 
                Season +
                coast_dist_m + 
                elevation_m +
                avg_sal_biomass + 
                avg_sal_biomass:Season + 
                recent_harvest_wshed +
                second_growth_wshed +
                old_growth_wshed +
                recent_harvest_wshed:Season +
                second_growth_wshed:Season +
                old_growth_wshed:Season +
                (1 | WTRSHD_FID/SiteID) +
                offset(monthly_effort_log),
              data = gb_data,
              ziformula = ~ 0,
              family = binomial)


# Build AIC table
gb.cand.models <- list("null model (1)" = G1,
                       "environmental model (2)" = G2, 
                       "forest (stand) + environmental model (3)" = G3,
                       "forest (watershed) + environmental model (4)" = G4,
                       "forest (watershed) w/ interactions + environmental model (5)" = G5)
gb_selectionTable <- aictab(cand.set = gb.cand.models, second.ord = FALSE)

## MODEL SELECTION ## ----------------------------------------------------------
wolf_selectionTable #(5)

deer_selectionTable #(5)

bb_selectionTable #(5) 

gb_selectionTable #(5) 

## MODEL SUMMARIES FOR TOP MODELS ## -------------------------------------------
# Wolf 
summary(W5)

# Deer
summary(D5)

# Black Bear 
summary(B5)

# Grizzly Bear 
summary(G5)

## DIAGNOSTICS ## --------------------------------------------------------------

## Calculate R2 ##
r2_values_wolf5 <- r2_nakagawa(W5)
print(r2_values_wolf5) # CR2 = 0.681, MR2 = 0.450

r2_values_deer5 <- r2_nakagawa(D5)
print(r2_values_deer5) # CR2 = 0.669, MR2 = 0.335

r2_values_bb5 <- r2_nakagawa(B5)
print(r2_values_bb5) # CR2 = 0.388, MR2 = 0.192

r2_values_gb5 <- r2_nakagawa(G5)
print(r2_values_gb5) # CR2 = 0.731, MR2 = 0.416

## Determine RE importance ## 
vc_wolf <- VarCorr(W5)
vc_wolf
vc_wolf$cond
sapply(vc_wolf$cond, function(x) attr(x, "stddev")^2)
vars_wolf <- sapply(vc_wolf$cond, function(x) attr(x, "stddev")^2)
total_var_wolf <- sum(vars_wolf)
vars_wolf / total_var_wolf # Wshed/SiteID = 0.8811349, Wshed = 0.1188651

vc_deer <- VarCorr(D5)
vc_deer
vc_deer$cond
sapply(vc_deer$cond, function(x) attr(x, "stddev")^2)
vars_deer <- sapply(vc_deer$cond, function(x) attr(x, "stddev")^2)
total_var_deer <- sum(vars_deer)
vars_deer / total_var_deer # Wshed/SiteID = 0.6068752, Wshed = 0.3931248 

vc_bb <- VarCorr(B5)
vc_bb
vc_bb$cond
sapply(vc_bb$cond, function(x) attr(x, "stddev")^2)
vars_bb <- sapply(vc_bb$cond, function(x) attr(x, "stddev")^2)
total_var_bb <- sum(vars_bb)
vars_bb / total_var_bb # Wshed/SiteID = 0.6222747, Wshed = 0.3777253 

vc_gb <- VarCorr(G5)
vc_gb
vc_gb$cond
sapply(vc_gb$cond, function(x) attr(x, "stddev")^2)
vars_gb <- sapply(vc_gb$cond, function(x) attr(x, "stddev")^2)
total_var_gb <- sum(vars_gb)
vars_gb / total_var_gb # Wshed/SiteID = 0.4996558, Wshed = 0.5003442 

## Residual Plots ## 
simres_wolf5 <-  simulateResiduals(fittedModel = W5, plot = TRUE) 
out <- testOutliers(simres_wolf5)
testDispersion(simres_wolf5)
plotResiduals(simres_wolf5, form = fitted(W5))

simres_deer5 <-  simulateResiduals(fittedModel = D5, plot = TRUE)
plotResiduals(simres_deer5, form = fitted(D5))

simres_bb5 <-  simulateResiduals(fittedModel = B5, plot = TRUE)
testOutliers(simres_bb5)
testUniformity(simres_bb5)
testZeroInflation(simres_bb5)
testDispersion(simres_bb5)
plotResiduals(simres_bb5, form = fitted(B5))

simres_gb5 <-  simulateResiduals(fittedModel = G5, plot = TRUE)
out <- testOutliers(simres_gb5)
out$observations
testDispersion(simres_gb5)
plotResiduals(simres_gb5, form = fitted(G5))

## Compute AUC ## 
wolf_data$predicted_probs <- predict(W5, type = "response")
auc_result_wolf <- roc(wolf_data$Detections_binom, wolf_data$predicted_probs)
print(auc_result_wolf) # AUC = 0.9121
plot(auc_result_wolf, col = "blue", lwd = 2, main = "Wolf ROC Curve")

deer_data$predicted_probs <- predict(D5, type = "response")
auc_result_deer <- roc(deer_data$Detections_binom, deer_data$predicted_probs)
print(auc_result_deer) # AUC = 0.8985
plot(auc_result_deer, col = "blue", lwd = 2, main = "Deer ROC Curve")

bb_data$predicted_probs <- predict(B5, type = "response")
auc_result_bb <- roc(bb_data$Detections_binom, bb_data$predicted_probs)
print(auc_result_bb) # AUC = 0.8316
plot(auc_result_bb, col = "blue", lwd = 2, main = "Black Bear ROC Curve")

gb_data$predicted_probs <- predict(G5, type = "response")
auc_result_gb <- roc(gb_data$Detections_binom, gb_data$predicted_probs)
print(auc_result_gb) # AUC = 0.9361
plot(auc_result_gb, col = "blue", lwd = 2, main = "Grizzly ROC Curve")

## VISUALIZE INTERACTIONS ## ---------------------------------------------------
plot(ggpredict(W5, terms = c("recent_harvest_wshed", "Season")))

plot(ggpredict(D5, terms = c("avg_sal_biomass", "Season")))
plot(ggpredict(D5, terms = c("old_growth_wshed", "Season")))

plot(ggpredict(B5, terms = c("avg_sal_biomass", "Season")))
plot(ggpredict(B5, terms = c("old_growth_wshed", "Season")))

plot(ggpredict(G5, terms = c("avg_sal_biomass", "Season")))
plot(ggpredict(G5, terms = c("recent_harvest_wshed", "Season")))
plot(ggpredict(G5, terms = c("second_growth_wshed", "Season")))


