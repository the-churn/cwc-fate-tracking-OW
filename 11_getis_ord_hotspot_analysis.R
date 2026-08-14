# ==============================================================================
# Script Name:  11_getis_ord_hotspot_analysis.R
# Description:  Performs per-species Getis-Ord Gi* spatial hotspot analysis on
#               colonization dates across vertical canyon wall coordinates.
#               Exports GIS spatial outputs and generates a 3-panel stacked, 
#               large-format publication figure with a single universal legend.
# Dependencies: sf, spdep, ggplot2, dplyr, readr, patchwork, grid
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(sf)
library(spdep)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(grid)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

input_csv <- file.path(base_dir, "results", "Master_Dataset_Hotspot_Flipped.csv")
if (!file.exists(input_csv)) {
  input_csv <- file.path(base_dir, "results", "Master_Dataset_Final_MonteCarlo_Results.csv")
}

results_csv <- file.path(base_dir, "results", "Master_Dataset_GetisOrd_Results.csv")
output_fig  <- file.path(base_dir, "figures", "Figure5_GetisOrd_Hotspots_2Col.png")

dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

K_NEIGHBORS   <- 8
valid_species <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

# ------------------------------------------------------------------------------
# 3. LOAD & PREPARE SPATIAL DATA (METER CONVERSION)
# ------------------------------------------------------------------------------
cat("Loading dataset for spatial analysis from:\n  ", input_csv, "\n")
raw_data <- read_csv(input_csv, show_col_types = FALSE)

x_col <- "Centroid_x"
y_col <- "Centroid_y_flipped"

cat("Using spatial coordinates: X =", x_col, "| Y =", y_col, "\n")

spatial_df <- raw_data %>%
  filter(!is.na(Colonization_Year_Est), !is.na(.data[[x_col]]), !is.na(.data[[y_col]])) %>%
  filter(Species %in% valid_species) %>%
  mutate(
    # Convert mm to meters if raw units are > 500
    x_m = if (max(.data[[x_col]], na.rm = TRUE) > 500) .data[[x_col]] / 1000 else .data[[x_col]],
    y_m = if (max(.data[[y_col]], na.rm = TRUE) > 500) .data[[y_col]] / 1000 else .data[[y_col]]
  )

# ------------------------------------------------------------------------------
# 4. GETIS-ORD Gi* CALCULATION LOOP (PER SPECIES)
# ------------------------------------------------------------------------------
compute_species_gi <- function(sp_name, df, k_nn = 8) {
  sub_df <- df %>% filter(Species == sp_name)
  n_pts  <- nrow(sub_df)
  
  if (n_pts < (k_nn + 1)) {
    warning("Not enough points for ", sp_name, " to construct ", k_nn, " neighbors. Skipping.")
    return(NULL)
  }
  
  coords <- as.matrix(sub_df[, c("x_m", "y_m")])
  
  knn_nb  <- knearneigh(coords, k = k_nn)
  nb      <- knn2nb(knn_nb)
  nb_self <- include.self(nb)
  lw_self <- nb2listw(nb_self, style = "B", zero.policy = TRUE)
  
  gi_vec <- localG(sub_df$Colonization_Year_Est, lw_self, zero.policy = TRUE)
  gi_z   <- as.numeric(gi_vec)
  
  sub_df %>%
    mutate(
      Gi_Zscore = gi_z,
      Gi_Class  = case_when(
        gi_z >=  2.58 ~ "Hotspot (99% CI)",
        gi_z >=  1.96 & gi_z <  2.58 ~ "Hotspot (95% CI)",
        gi_z >=  1.65 & gi_z <  1.96 ~ "Hotspot (90% CI)",
        gi_z <= -2.58 ~ "Coldspot (99% CI)",
        gi_z <= -1.96 & gi_z > -2.58 ~ "Coldspot (95% CI)",
        gi_z <= -1.65 & gi_z > -1.96 ~ "Coldspot (90% CI)",
        TRUE ~ "Not Significant"
      ),
      Gi_Class = factor(Gi_Class, levels = c(
        "Hotspot (99% CI)", "Hotspot (95% CI)", "Hotspot (90% CI)",
        "Not Significant",
        "Coldspot (90% CI)", "Coldspot (95% CI)", "Coldspot (99% CI)"
      ))
    )
}

cat("\nComputing Getis-Ord Gi* across species...\n")
gi_results_list <- lapply(valid_species, function(sp) {
  cat("  -> Processing:", sp, "\n")
  compute_species_gi(sp, spatial_df, k_nn = K_NEIGHBORS)
})

final_gi_df <- bind_rows(gi_results_list)
write_csv(final_gi_df, results_csv)
cat("Exported GIS spatial results to:\n  ", results_csv, "\n")

# Determine shared global bounding box for perfect vertical alignment
x_min <- 0
x_max <- max(final_gi_df$x_m, na.rm = TRUE) * 1.02  # ~15.2 m
y_min <- min(final_gi_df$y_m, na.rm = TRUE) * 0.90  # ~0.4 m
y_max <- max(final_gi_df$y_m, na.rm = TRUE) * 1.08  # ~2.8 m

# ------------------------------------------------------------------------------
# 5. BUILD LARGE PUBLICATION FIGURE
# ------------------------------------------------------------------------------
hotspot_colors <- c(
  "Hotspot (99% CI)"   = "#D73027",  # Deep Red
  "Hotspot (95% CI)"   = "#FC8D59",  # Medium Orange
  "Hotspot (90% CI)"   = "#FEE090",  # Light Yellow
  "Not Significant"    = "#E5E7EB",  # Muted Light Grey
  "Coldspot (90% CI)"  = "#E0F3F8",  # Soft Cyan
  "Coldspot (95% CI)"  = "#91BFDB",  # Medium Blue
  "Coldspot (99% CI)"  = "#4575B4"   # Deep Navy
)

species_labels <- c(
  "Madrepora oculata"     = "Madrepora oculata",
  "Desmophyllum pertusum" = "Desmophyllum pertusum",
  "Primnoa msp."          = "Primnoa msp."
)

make_species_panel <- function(sp_name, show_x = FALSE) {
  p_df   <- final_gi_df %>% filter(Species == sp_name)
  bg_df  <- p_df %>% filter(Gi_Class == "Not Significant")
  sig_df <- p_df %>% filter(Gi_Class != "Not Significant")
  
  p <- ggplot() +
    # Faint non-significant background points
    geom_point(
      data = bg_df,
      aes(x = x_m, y = y_m),
      color = "#D1D5DB", alpha = 0.4, size = 1.2, stroke = 0
    ) +
    # Highlighted Hotspots/Coldspots (larger points with crisp dark outlines)
    geom_point(
      data = sig_df,
      aes(x = x_m, y = y_m, fill = Gi_Class),
      shape = 21, color = "#111827", stroke = 0.35, size = 2.2, alpha = 0.95
    ) +
    scale_fill_manual(
      values = hotspot_colors,
      drop = FALSE,
      name = "Getis-Ord Gi* Spatial Confidence"
    ) +
    # Lock coordinates to identical bounding box WITHOUT aspect ratio squishing
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    # Inset species title directly inside top-left of plot area
    annotate(
      "text", x = x_min + 0.3, y = y_max - 0.25,
      label = species_labels[[sp_name]],
      fontface = "bold.italic", size = 3.3, hjust = 0, color = "black"
    ) +
    labs(
      x = if (show_x) "Horizontal Position (m)" else NULL,
      y = "Elevation (m)"
    ) +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(
      axis.title        = element_text(size = 8.5, face = "bold", color = "black"),
      axis.text         = element_text(size = 7.5, color = "black"),
      axis.line         = element_line(linewidth = 0.4, color = "black"),
      axis.ticks        = element_line(linewidth = 0.4, color = "black"),
      panel.background  = element_rect(fill = "#FAFAFA", color = "#D1D5DB"),
      plot.margin       = margin(t = 3, r = 8, b = 3, l = 5, unit = "pt")
    )
  
  return(p)
}

cat("\nAssembling high-visibility stacked figure...\n")

panel_madr <- make_species_panel("Madrepora oculata", show_x = FALSE)
panel_desm <- make_species_panel("Desmophyllum pertusum", show_x = FALSE)
panel_prim <- make_species_panel("Primnoa msp.", show_x = TRUE)

# Combine 3 panels vertically with ONE universal horizontal legend at the bottom
combined_figure <- (panel_madr / panel_desm / panel_prim) +
  plot_layout(guides = "collect") &
  theme(
    legend.position  = "bottom",
    legend.direction = "horizontal",
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    legend.margin    = margin(t = 5, b = 2, unit = "pt"),
    plot.tag         = element_text(size = 11, face = "bold", family = "sans")
  ) &
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(size = 3)
    )
  )

# Add Panel A, B, C tags cleanly
combined_figure <- combined_figure +
  plot_annotation(tag_levels = 'A')

# ------------------------------------------------------------------------------
# 6. EXPORT HIGH-RESOLUTION PRINT FIGURE
# ------------------------------------------------------------------------------
ggsave(
  filename = output_fig,
  plot     = combined_figure,
  width    = 178,     # Full 2-column width for Cell Press / Current Biology (mm)
  height   = 160,     # Increased height so panels are tall, wide, and clear
  units    = "mm",
  dpi      = 600
)

cat("Successfully generated and exported Getis-Ord Gi* figure to:\n  ", output_fig, "\n")
cat("======================================================================\n")
cat("Script 11 Execution Complete!\n")
