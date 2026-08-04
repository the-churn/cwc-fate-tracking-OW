# ==============================================================================
# Script Name:  08_figure3_growth_and_demographics.R
# Description:  Assembles the two-panel Figure 3 for Current Biology: (A) size-
#               at-age growth trajectories with 95% Monte Carlo CI, and (B) a
#               ridgeline of population age demographics for 2022. Combines
#               both panels with patchwork and exports at the journal's
#               2-column (178 mm) print specification.
# Dependencies: ggplot2, ggridges, dplyr, readr, patchwork, grid
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggridges)
library(dplyr)
library(readr)
library(patchwork)
library(grid)   # unit() used in theme() below

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
curves_csv          <- file.path("outputs", "results", "size_age_curves.csv")
observed_csv         <- file.path("outputs", "results", "observed_points.csv")
final_results_csv    <- file.path("outputs", "results", "Master_Dataset_Final_MonteCarlo_Results.csv")
figures_dir           <- file.path("outputs", "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. LOAD UPSTREAM OUTPUTS (SCRIPTS 06-07)
# ------------------------------------------------------------------------------
required_files <- c(curves_csv, observed_csv, final_results_csv)
if (any(!file.exists(required_files))) {
  stop("Missing upstream outputs. Run 06_monte_carlo_age_simulation.R and ",
       "07_size_age_curves.R first.", call. = FALSE)
}

size_age_curves      <- read_csv(curves_csv, show_col_types = FALSE)
observed_points        <- read_csv(observed_csv, show_col_types = FALSE)
final_master_dataset <- read_csv(final_results_csv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 4. UNIFIED COLOR PALETTE (SHARED ACROSS BOTH PANELS)
# ------------------------------------------------------------------------------
coral_colors <- c(
  "Desmophyllum pertusum" = "#7A1F5C",  # dark pinkish purple
  "Madrepora oculata"     = "#3B0F4B",  # deep purple
  "Primnoa msp."          = "#F69336"   # warm orange
)

coral_labels <- c(
  "Desmophyllum pertusum" = expression(italic("D. pertusum")),
  "Madrepora oculata"     = expression(italic("M. oculata")),
  "Primnoa msp."          = expression(italic("Primnoa") ~ "msp.")
)

# ------------------------------------------------------------------------------
# 5. PANEL A -- SIZE-AT-AGE GROWTH TRAJECTORIES
# ------------------------------------------------------------------------------
# Ontogenetic collapse-threshold ages (D. pertusum, M. oculata, Primnoa msp.):
# the estimated age at which each species reaches its 35.9 / 59.9 / 187.8 cm^2
# collapse-size threshold. Carried over as fixed values from the source
# notebook -- these were originally interpolated from size_age_curves against
# those thresholds. Recompute directly if the underlying data/model changes
# (happy to wire that up as its own step if useful).
collapse_threshold_ages <- c(11.91, 21.54, 20.18)

panel_a <- ggplot(size_age_curves, aes(x = Age_Med, y = Size, color = Species, fill = Species)) +
  geom_point(
    data = observed_points,
    aes(x = Age, y = Size, color = Species),
    alpha = 0.25, size = 0.8, inherit.aes = FALSE
  ) +
  geom_vline(
    xintercept = collapse_threshold_ages,
    linetype = "dashed",
    color = c("#7A1F5C", "#3B0F4B", "#F69336"),
    alpha = 0.6, linewidth = 0.4
  ) +
  geom_ribbon(aes(xmin = Age_Lo95, xmax = Age_Hi95), alpha = 0.22, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1,000"),
    expand = c(0.02, 0)
  ) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_color_manual(values = coral_colors, labels = coral_labels) +
  scale_fill_manual(values = coral_colors, labels = coral_labels) +
  labs(
    x = "Estimated Age (years)",
    y = expression(paste("Colony Area (cm"^2*", log scale)"))
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.title        = element_text(size = 8, color = "black"),
    axis.text         = element_text(size = 7, color = "black"),
    axis.line         = element_line(linewidth = 0.4, color = "black"),
    axis.ticks        = element_line(linewidth = 0.4, color = "black"),
    legend.position   = c(0.82, 0.25),
    legend.title      = element_blank(),
    legend.text       = element_text(size = 7),
    legend.key.size   = unit(0.4, "cm"),
    legend.background = element_blank(),
    plot.margin       = margin(t = 5, r = 8, b = 5, l = 5, unit = "pt")
  )

# ------------------------------------------------------------------------------
# 6. PANEL B -- POPULATION AGE DEMOGRAPHICS (2022 RIDGELINE)
# ------------------------------------------------------------------------------
ridge_data <- final_master_dataset %>%
  filter(!is.na(Age_2022_Median), Age_2022_Median > 0, !is.na(Species)) %>%
  mutate(
    Species_Formatted = factor(Species, levels = c(
      "Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp."
    ))
  )

medians_df <- ridge_data %>%
  group_by(Species_Formatted) %>%
  summarise(med_x = round(median(Age_2022_Median, na.rm = TRUE), 1), .groups = "drop") %>%
  mutate(
    y_start = as.numeric(Species_Formatted),
    y_end   = as.numeric(Species_Formatted) + 0.85
  )

panel_b <- ggplot(ridge_data, aes(x = Age_2022_Median, y = Species_Formatted, fill = Species_Formatted)) +
  geom_density_ridges(
    scale = 1.2, rel_min_height = 0.005, alpha = 0.65,
    linewidth = 0.4, color = "black"
  ) +
  geom_segment(
    data = medians_df,
    aes(x = med_x, xend = med_x, y = y_start, yend = y_end),
    color = "black", linewidth = 0.6, inherit.aes = FALSE
  ) +
  geom_density_ridges(
    aes(point_color = Species_Formatted),
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.2, height = 0),
    point_shape = 21, point_size = 0.6, point_fill = "white",
    point_alpha = 0.25, alpha = 0, scale = 1.2, stroke = 0.3
  ) +
  scale_fill_manual(values = coral_colors) +
  scale_discrete_manual(aesthetics = "point_color", values = coral_colors) +
  scale_x_continuous(breaks = seq(0, 60, by = 20), limits = c(0, 65), expand = c(0.02, 0)) +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.25))) +
  labs(x = "Estimated Colony Age in 2022 (years)", y = NULL) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 8, color = "black"),
    axis.text.x  = element_text(size = 7, color = "black"),
    axis.text.y  = element_text(size = 7, face = "italic", color = "black"),
    axis.line    = element_line(linewidth = 0.4, color = "black"),
    axis.ticks   = element_line(linewidth = 0.4, color = "black"),
    legend.position = "none",
    plot.margin  = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")
  )

# ------------------------------------------------------------------------------
# 7. COMBINE PANELS & EXPORT (CURRENT BIOLOGY 2-COLUMN: 178 MM)
# ------------------------------------------------------------------------------
cat("Assembling final two-panel figure...\n")

combined_figure <- (panel_a | panel_b) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(plot.tag = element_text(size = 10, face = "bold", family = "sans"))
  )

print(combined_figure)

output_fig <- file.path(figures_dir, "Figure1_Growth_and_Demographics_2Col.PNG")

ggsave(
  filename = output_fig,
  plot     = combined_figure,
  width    = 178,
  height   = 65,
  units    = "mm",
  dpi      = 600,
  type     = "cairo"
)

cat("Successfully generated and exported:\n  ", output_fig, "\n")
