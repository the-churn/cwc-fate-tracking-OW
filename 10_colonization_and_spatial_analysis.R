# ==============================================================================
# Script Name:  10_colonization_and_spatial_analysis.R
# Description:  Generates colonization year demographic figures (Figure 4), 
#               performs robust non-parametric statistical testing (Kruskal-
#               Wallis, Dunn's post-hoc, quantile regressions), and conducts 
#               spatial autocorrelation analyses (Global Moran's I) with 
#               Y-coordinate flipping and Monte Carlo 95% CI exports for GIS.
# Dependencies: ggplot2, ggridges, dplyr, readr, quantreg, dunn.test, spdep, sf
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES & SET DIRECTORIES
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggridges)
library(dplyr)
library(readr)
library(quantreg)
library(dunn.test)
library(spdep)
library(sf)

select <- dplyr::select  # Avoid masking by other packages

# Base Directory Setup
base_dir    <- "D:/PhD_Data(Large)/Submission_Dataset"
results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
tables_dir  <- file.path(base_dir, "tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# File Paths
final_results_csv <- file.path(results_dir, "Master_Dataset_Final_MonteCarlo_Results.csv")
fig_out_path      <- file.path(figures_dir, "Figure4_Colonization_Year_Ridgeline.png")
flipped_csv_out   <- file.path(results_dir, "Master_Dataset_Hotspot_Flipped.csv")

# Palette & Species Setup
main_cwc <- c("Madrepora oculata", "Primnoa msp.", "Desmophyllum pertusum")

coral_colors <- c(
  "Desmophyllum pertusum" = "#7A1F5C",  # Dark pinkish purple
  "Madrepora oculata"     = "#3B0F4B",  # Deep purple
  "Primnoa msp."          = "#F69336"   # Warm orange
)

# ------------------------------------------------------------------------------
# 2. LOAD & CLEAN DATA
# ------------------------------------------------------------------------------
if (!file.exists(final_results_csv)) {
  stop("Missing master Monte Carlo dataset. Please run Script 07 first.", call. = FALSE)
}

cat("Loading master dataset...\n")
colonization_year <- read_csv(final_results_csv, show_col_types = FALSE)

# Clean dataset for colonization year analysis
year_col <- "Colonization_Year_Est"
df_clean <- colonization_year %>%
  filter(Species %in% main_cwc, !is.na(.data[[year_col]])) %>%
  mutate(
    Species = factor(Species, levels = rev(main_cwc)),
    Colonization_Year_Est = as.numeric(.data[[year_col]])
  )

# ------------------------------------------------------------------------------
# 3. FIGURE 4: COLONIZATION YEAR RIDGELINE PLOT
# ------------------------------------------------------------------------------
cat("Processing Ridgeline geometry and percentile lines...\n")

ridge_scale <- 1.4

# Calculate dynamic X-axis limits and breaks based on true data range
min_data_year <- min(df_clean[[year_col]], na.rm = TRUE)
x_min_break   <- floor(min_data_year / 10) * 10  # Rounds down to nearest decade
x_max_break   <- 2020

# Use 20-year intervals if spanning >= 80 years to prevent label overlap, else 10-year
x_by_step     <- if ((x_max_break - x_min_break) >= 80) 20 else 10
x_breaks      <- seq(x_min_break, x_max_break, by = x_by_step)

# Calculate percentiles and medians per species
stats <- df_clean %>%
  group_by(Species) %>%
  summarise(
    n   = n(),
    p05 = quantile(.data[[year_col]], 0.05, na.rm = TRUE),
    p10 = quantile(.data[[year_col]], 0.10, na.rm = TRUE),
    p50 = median(.data[[year_col]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(Species, " (n = ", n, ")")) %>%
  arrange(Species)

stats$label <- factor(stats$label, levels = stats$label)
stats$y0    <- as.numeric(stats$label)

label_map <- setNames(as.character(stats$label), as.character(stats$Species))
df_clean$label <- factor(label_map[as.character(df_clean$Species)], levels = levels(stats$label))
fill_colors    <- setNames(coral_colors[as.character(stats$Species)], stats$label)

# Render base ridge to extract density curve boundary dynamically
base_ridge <- ggplot(df_clean, aes(x = .data[[year_col]], y = label, fill = label)) +
  geom_density_ridges(scale = ridge_scale, rel_min_height = 0.005)

ridge_geom <- ggplot_build(base_ridge)$data[[1]]
grp_to_y0  <- ridge_geom %>% group_by(group) %>% summarise(y0 = round(min(ymin), 3), .groups = "drop")
stats      <- stats %>% left_join(grp_to_y0, by = "y0")

get_curve_top <- function(grp, xval) {
  sub_df <- ridge_geom[ridge_geom$group == grp, ]
  sub_df <- sub_df[order(sub_df$x), ]
  approx(sub_df$x, sub_df$ymax, xout = xval, rule = 2)$y
}

stats <- stats %>%
  rowwise() %>%
  mutate(
    p05_top = get_curve_top(group, p05),
    p10_top = get_curve_top(group, p10),
    p50_top = get_curve_top(group, p50)
  ) %>%
  ungroup()

# Plot Ridgeline
fig4_plot <- ggplot(df_clean, aes(x = .data[[year_col]], y = label, fill = label)) +
  geom_density_ridges(scale = ridge_scale, rel_min_height = 0.005, colour = "black", linewidth = 0.5) +
  # 5th Percentile: Dotted Crimson Line
  geom_segment(data = stats, aes(x = p05, xend = p05, y = y0, yend = p05_top),
               inherit.aes = FALSE, colour = "#B3123E", linetype = "dotted", linewidth = 0.8) +
  # 10th Percentile: Dashed Crimson Line
  geom_segment(data = stats, aes(x = p10, xend = p10, y = y0, yend = p10_top),
               inherit.aes = FALSE, colour = "#B3123E", linetype = "dashed", linewidth = 0.8) +
  # Median: Solid White Line
  geom_segment(data = stats, aes(x = p50, xend = p50, y = y0, yend = p50_top),
               inherit.aes = FALSE, colour = "white", linewidth = 1.0) +
  # Median Label Box
  geom_label(data = stats, aes(x = p50, y = y0 + (p50_top - y0) * 0.75, label = round(p50), colour = label),
             inherit.aes = FALSE, fill = "white", fontface = "bold", size = 3.5,
             label.padding = unit(0.2, "lines"), label.size = NA) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  scale_colour_manual(values = fill_colors, guide = "none") +
  scale_x_continuous(
    limits = c(x_min_break - 2, 2025),
    breaks = x_breaks,
    expand = c(0.01, 0)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.1))) +
  labs(x = "Estimated Colonization Year", y = NULL) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(
    axis.text.y        = element_text(face = "italic", size = 10, colour = "black"),
    axis.title.x       = element_text(size = 10, hjust = 0.5, margin = margin(t = 8)),
    axis.text.x        = element_text(size = 9, colour = "black"),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(10, 15, 10, 10)
  )

ggsave(fig_out_path, plot = fig4_plot, width = 178, height = 120, units = "mm", dpi = 600, bg = "white")
cat("Saved Figure 4 to:", fig_out_path, "\n")

# ------------------------------------------------------------------------------
# 4. STATISTICAL DEMOGRAPHIC ANALYSES
# ------------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("STATISTICAL DEMOGRAPHIC TESTS (COLONIZATION YEAR)\n")
cat("======================================================================\n")

# 4A. Kruskal-Wallis & Dunn's Post-Hoc Test
kw_test <- kruskal.test(Colonization_Year_Est ~ Species, data = df_clean)
print(kw_test)

cat("\nDunn's Post-Hoc Test (Benjamini-Hochberg Adjusted):\n")
dunn_res <- dunn.test(x = df_clean$Colonization_Year_Est, g = df_clean$Species, method = "bh", altp = TRUE)

# 4B. Quantile Regressions (Baseline: Desmophyllum pertusum)
cat("\nQuantile Regression (5th Percentile, tau = 0.05) [Ref: D. pertusum]:\n")
q05_fit <- rq(Colonization_Year_Est ~ Species, data = df_clean, tau = 0.05, method = "fn")
print(summary(q05_fit, se = "boot", R = 1000))

cat("\nQuantile Regression (10th Percentile, tau = 0.10) [Ref: D. pertusum]:\n")
q10_fit <- rq(Colonization_Year_Est ~ Species, data = df_clean, tau = 0.10, method = "fn")
print(summary(q10_fit, se = "boot", R = 1000))

# 4C. Quantile Regressions (Direct Contrast Ref: Primnoa msp.)
cat("\nRe-leveling factor reference to Primnoa msp. for direct contrast...\n")
df_relevel <- df_clean %>%
  mutate(Species = relevel(factor(Species), ref = "Primnoa msp."))

cat("\nQuantile Regression (5th Percentile, tau = 0.05) [Ref: Primnoa msp.] :\n")
q05_primnoa_ref <- rq(Colonization_Year_Est ~ Species, data = df_relevel, tau = 0.05, method = "fn")
print(summary(q05_primnoa_ref, se = "boot", R = 1000))

cat("\nQuantile Regression (10th Percentile, tau = 0.10) [Ref: Primnoa msp.] :\n")
q10_primnoa_ref <- rq(Colonization_Year_Est ~ Species, data = df_relevel, tau = 0.10, method = "fn")
print(summary(q10_primnoa_ref, se = "boot", R = 1000))

# ------------------------------------------------------------------------------
# 5. SPATIAL AUTOCORRELATION & GIS EXPORT
# ------------------------------------------------------------------------------
if (all(c("Centroid_x", "Centroid_y") %in% colnames(colonization_year))) {
  cat("\n======================================================================\n")
  cat("SPATIAL ANALYSIS & GIS EXPORT\n")
  cat("======================================================================\n")
  
  # Y-Axis inversion for photogrammetry coordinates
  max_y <- max(colonization_year$Centroid_y, na.rm = TRUE)
  min_y <- min(colonization_year$Centroid_y, na.rm = TRUE)
  
  colonization_flipped <- colonization_year %>%
    filter(!is.na(Colonization_Year_Est), !is.na(Centroid_x), !is.na(Centroid_y)) %>%
    mutate(
      Centroid_y_flipped        = max_y - Centroid_y + min_y,
      Colonization_Year_Est    = as.numeric(Colonization_Year_Est),
      Age_2022_Median          = as.numeric(Age_2022_Median),
      
      # Inverted Monte Carlo 95% Confidence Bounds for Colonization Year
      Colonization_Year_Lower95  = 2022 - as.numeric(Age_2022_Upper95),
      Colonization_Year_Upper95  = 2022 - as.numeric(Age_2022_Lower95),
      Colonization_Year_CI_Range = Colonization_Year_Upper95 - Colonization_Year_Lower95
    )
  
  # Write clean dataset ready for ArcGIS Pro / QGIS
  write_csv(colonization_flipped, flipped_csv_out)
  cat("Saved Y-flipped dataset with 95% CIs for GIS to:", flipped_csv_out, "\n")
  
  # Global Moran's I per species (k = 5 nearest neighbors)
  for (sp in main_cwc) {
    sub_sp <- colonization_flipped %>% filter(Species == sp)
    if (nrow(sub_sp) > 6) {
      coords    <- cbind(sub_sp$Centroid_x, sub_sp$Centroid_y_flipped)
      knn       <- knearneigh(coords, k = 5)
      weights   <- nb2listw(knn2nb(knn), style = "W")
      moran_res <- moran.test(sub_sp$Colonization_Year_Est, weights)
      
      cat("\nGlobal Moran's I -", sp, ":\n")
      print(moran_res)
    }
  }
} else {
  cat("\nSpatial centroid columns not found. Skipping Moran's I and GIS flips.\n")
}

cat("\n======================================================================\n")
cat("Script 10 Execution Complete!\n")
