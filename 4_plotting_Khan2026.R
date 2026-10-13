# Author: Persia Khan (persia.khan@gmail.com)
# Date Last Modified: 2026-01-12

################################################################################
# 4) PLOTTING - The following code was used to produce the figures in Khan et al. 
# (2026). This script builds on the model outputs generated in 
# '3_modelling_Khan2026'. Figures are listed in the order the appear within the 
# manuscript. 
################################################################################

library(ggplot2)
library(patchwork)

# Figure 2 - Coefficient Plot --------------------------------------------------

# Setup dataframes by species and bind for plotting
df1 <- tidy(D5, effects = "fixed") %>%
  mutate(Species = "Deer")

df2 <- tidy(W5, effects = "fixed") %>%
  mutate(Species = "Wolf")

df3 <- tidy(B5, effects = "fixed") %>%
  mutate(Species = "Black Bear")

df4 <- tidy(G5, effects = "fixed") %>%
  mutate(Species = "Grizzly Bear")

df_results <- bind_rows(df1, df2, df3, df4)

df_results$Species <- factor(df_results$Species, levels = c("Deer", "Wolf", "Black Bear", "Grizzly Bear"))

df_results <- df_results %>%
  dplyr::rename(Covariate = term, Estimate = estimate, SE = std.error, p_value = p.value) %>%
  dplyr::select(Species, Covariate, Estimate, SE, p_value)

df_results <- df_results %>%
  mutate(Significance = case_when(
    p_value < 0.05 ~ "p < 0.05",
    p_value >= 0.05 & p_value < 0.10 ~ "p = 0.05-0.10",
    p_value >= 0.10 ~ "p > 0.10"
  ))

df_results <- df_results %>%
  mutate(Significance = factor(Significance, levels = c("p < 0.05", "p = 0.05-0.10", "p > 0.10")))

df_results_filtered <- df_results %>%
  filter(Covariate %in% c("recent_harvest_wshed", 
                          "second_growth_wshed", 
                          "old_growth_wshed",
                          "coast_dist_m", 
                          "elevation_m", 
                          "SeasonSpawning Season (Sept-Nov):avg_sal_biomass",
                          "avg_sal_biomass:SeasonSpawning Season (Sept-Nov)"))

df_results_filtered <- df_results_filtered %>%
  mutate(Lower = Estimate - (1.96 * SE),  
         Upper = Estimate + (1.96 * SE))

# Build ggplot 
ggplot(df_results_filtered, aes(x = factor(Covariate, levels = c("recent_harvest_wshed", 
                                                                 "second_growth_wshed", 
                                                                 "old_growth_wshed", 
                                                                 "coast_dist_m", 
                                                                 "elevation_m",
                                                                 "SeasonSpawning Season (Sept-Nov):avg_sal_biomass",
                                                                 "avg_sal_biomass:SeasonSpawning Season (Sept-Nov)")),  
                                y = Estimate, shape = Species, color = Significance)) +  
  geom_point(position = position_dodge(width = 0.6), size = 9) + 
  geom_linerange(aes(ymin = Lower, ymax = Upper),  
                 position = position_dodge(width = 0.6), size = 1) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.2) +  
  scale_color_manual(values = c("p < 0.05" = "black", "p = 0.05-0.10" = "#808084", "p > 0.10" = "#C4C4C4")) + 
  scale_shape_manual(values = c(17, 18, 16, 15)) +  
  scale_x_discrete(labels = c("recent_harvest_wshed" = "Prop. Recent Harvest\n in Watershed", 
                              "second_growth_wshed" = "Prop. Second Growth\n in Watershed", 
                              "old_growth_wshed" = "Prop. Old Growth\n in Watershed",
                              "coast_dist_m" = "Distance to\n Coastline", 
                              "elevation_m" = "Elevation", 
                              "SeasonSpawning Season (Sept-Nov):avg_sal_biomass" = "Log Salmon Biomass\n x Spawning Season",
                              "avg_sal_biomass:SeasonSpawning Season (Sept-Nov)" = "Log Salmon Biomass\n x Spawning Season")) +
  theme_classic() +
  labs(x = "Fixed Effects", y = "Regression Coefficient") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18),  
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 16, face = "bold"),  
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 20),  
        legend.title = element_text(size = 20, face = "bold"),
  )

# Figure 3 - Interactions Plot -------------------------------------------------

# Create list of relevant predictors
x_labels <- list(
  "avg_sal_biomass" = "Average Salmon Biomass",
  "recent_harvest_wshed" = "Proportion of Recent Harvest (Watershed)",
  "second_growth_wshed" = "Proportion of Second Growth (Watershed)",
  "old_growth_wshed" = "Proportion of Old Growth (Watershed)"
)

# Plot function
create_interaction_plot <- function(model, predictor, title, show_legend = FALSE) {

  t <- ggpredict(model, terms = c(paste0(predictor, " [all]"), "Season"))
  
  if (!predictor %in% scaling_info$variable) {
    stop(paste("Scaling info missing for", predictor))
  }

  mean_val <- scaling_info$mean[scaling_info$variable == predictor]
  sd_val <- scaling_info$sd[scaling_info$variable == predictor]
  t$x <- (t$x * sd_val) + mean_val

  x_label <- ifelse(predictor %in% names(x_labels), x_labels[[predictor]], predictor)

  legend_pos <- if (show_legend) c(0.65, 0.88) else "none"
  
  p <- ggplot(t, aes(x = x, y = predicted, color = group)) +
    geom_line(size = 1.15) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = legend_pos,
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = alpha("white", 0.1)),
      legend.title = element_text(face = "bold", size = 8.5),
      legend.text = element_text(size = 8.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),  
      axis.text.y = element_text(size = 11),
      plot.margin = margin(5, 5, 5, 5, "pt")
    ) +
    ylab("Probability of Monthly Occurrence") +
    xlab(x_label) +
    scale_color_manual(
      values = c(
        "Summer (Apr-Aug)" = "#b3ba9e",
        "Winter (Dec-Mar)" = "#9786b3",
        "Spawning Season (Sept-Nov)" = "#cc9360"
      ),
      name = "Season"
    ) +
    scale_fill_manual(
      values = c(
        "Summer (Apr-Aug)" = "#b3ba9e",
        "Winter (Dec-Mar)" = "#9786b3",
        "Spawning Season (Sept-Nov)" = "#cc9360"
      ),
      name = "Season"
    )
  
  return(p)
}

# Assign values from species-specific top models
species_models <- list(
  "Deer" = D5,
  "Wolf" = W5,
  "Black Bear" = B5,
  "Grizzly Bear" = G5
)

# Define species and their corresponding predictor
species_predictors <- list(
  "Deer" = "old_growth_wshed",
  "Wolf" = "recent_harvest_wshed",
  "Black Bear" = "old_growth_wshed",
  "Grizzly Bear" = "second_growth_wshed"
)

# Store plots as individual objects to use 'patchwork'
first_plot <- TRUE  # Flag to track the first plot for showing legend

for (species in names(species_predictors)) {
  predictor <- species_predictors[[species]]
  
  if (!species %in% names(species_models)) {
    warning(paste("❌ No model found for", species, "- Check name matching"))
    next
  }
  
  model <- species_models[[species]]
  print(paste("✅ Processing:", species))
  print(paste("Model class for", species, ":", class(model)))
  print(paste("Using predictor:", predictor))
  
  plot_obj <- create_interaction_plot(
    model = model,
    predictor = predictor,
    title = paste(species),
    show_legend = first_plot
  )
  
  first_plot <- FALSE  # Set to FALSE after the first plot
  
  if (!is.null(plot_obj)) {
    plot_name <- paste0(gsub(" ", "_", tolower(species)), "_plot")
    assign(plot_name, plot_obj, envir = .GlobalEnv)
    print(paste("✅ Stored plot as:", plot_name))
  } else {
    warning(paste("❌ Skipping plot for", species, "due to errors."))
  }
}

plot_list <- c("wolf_plot", "black_bear_plot", "grizzly_bear_plot", "deer_plot")

existing_plots <- lapply(plot_list, function(p) {
  if (exists(p, envir = .GlobalEnv)) {
    print(paste("✅ Adding to final plot:", p))  # Debugging
    get(p)
  } else {
    print(paste("❌ Missing plot:", p))
    NULL
  }
})

existing_plots <- Filter(Negate(is.null), existing_plots)

if (length(existing_plots) > 0) {
  final_plot <- wrap_plots(existing_plots, ncol = 2)
  print("✅ Final plot created successfully!")
  print(final_plot)
} else {
  warning("❌ No plots were created successfully.")
}

# Build stacked plots with desired spacing using 'patchwork'
stacked_plot_A <- deer_plot / plot_spacer() / black_bear_plot +
  plot_layout(heights = c(1, 0.02, 1)) 
stacked_plot_B <- wolf_plot / plot_spacer() / grizzly_bear_plot +
  plot_layout(heights = c(1, 0.02, 1)) 

final_plot <- stacked_plot_A | plot_spacer() | stacked_plot_B 
final_plot +
  plot_layout(widths = c(1, 0.02, 1)) 

# Figure 3 - Species x Salmon Plot ---------------------------------------------

# Assign values from species-specific top models (for three species considered)
species_models <- list(
  "Deer" = D5,
  "Black Bear" = B5,
  "Grizzly Bear" = G5
)

# Define species and their corresponding predictor
species_predictors <- list(
  "Deer" = "avg_sal_biomass",
  "Black Bear" = "avg_sal_biomass",
  "Grizzly Bear" = "avg_sal_biomass"
)

# Store plots as individual objects
first_plot <- TRUE  # Flag to track the first plot for showing legend

for (species in names(species_predictors)) {
  predictor <- species_predictors[[species]]
  
  if (!species %in% names(species_models)) {
    warning(paste("❌ No model found for", species, "- Check name matching"))
    next
  }
  
  model <- species_models[[species]]
  print(paste("✅ Processing:", species))
  print(paste("Model class for", species, ":", class(model)))
  print(paste("Using predictor:", predictor))
  
  plot_obj <- create_interaction_plot(
    model = model,
    predictor = predictor,
    title = paste(species),
    show_legend = first_plot
  )
  
  first_plot <- FALSE 
  
  if (!is.null(plot_obj)) {
    plot_name <- paste0(gsub(" ", "_", tolower(species)), "_plot")
    assign(plot_name, plot_obj, envir = .GlobalEnv)
    print(paste("✅ Stored plot as:", plot_name))
  } else {
    warning(paste("❌ Skipping plot for", species, "due to errors."))
  }
}

plot_list <- c("wolf_plot", "black_bear_plot", "grizzly_bear_plot", "deer_plot")

# Retrieve only existing plots
existing_plots <- lapply(plot_list, function(p) {
  if (exists(p, envir = .GlobalEnv)) {
    print(paste("✅ Adding to final plot:", p))  # Debugging
    get(p)
  } else {
    print(paste("❌ Missing plot:", p))
    NULL
  }
})

existing_plots <- Filter(Negate(is.null), existing_plots)

if (length(existing_plots) > 0) {
  final_plot <- wrap_plots(existing_plots, ncol = 2)
  print("✅ Final plot created successfully!")
  print(final_plot)
} else {
  warning("❌ No plots were created successfully.")
}

# Assemble plot panels
sal_plots <- grizzly_bear_plot | black_bear_plot | deer_plot
sal_plots


# TABLE 3 - Master Table -------------------------------------------------------

extract_estimates_with_ci <- function(model, model_name = "Model") {

  is_glmmTMB <- "glmmTMB" %in% class(model)
  
  if (is_glmmTMB) {
    coefs <- as.data.frame(summary(model)$coefficients$cond)
    ci <- confint(model, parm = "beta_", level = 0.95)
  } else {
    coefs <- as.data.frame(summary(model)$coefficients)
    ci <- confint(model)
  }

  summary_df <- data.frame(
    Parameter = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    CI_lower = ci[, 1],
    CI_upper = ci[, 2],
    row.names = NULL
  )
  
  summary_df$Estimate_CI <- sprintf("%.3f (%.3f, %.3f)", summary_df$Estimate, summary_df$CI_lower, summary_df$CI_upper)
  
  summary_df$Model <- model_name

  summary_df[, c("Model", "Parameter", "Estimate_CI")]
}

results_list <- list(
  extract_estimates_with_ci(W5, "W5"),
  extract_estimates_with_ci(D5, "D5"),
  extract_estimates_with_ci(B5, "B5"),
  extract_estimates_with_ci(G5, "G5")
)

all_results <- do.call(rbind, results_list)

print(all_results)

wide_results <- all_results %>%
  pivot_wider(
    names_from = Model,
    values_from = Estimate_CI
  )

print(wide_results)

