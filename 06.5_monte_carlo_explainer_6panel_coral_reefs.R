# ==============================================================================
# Script Name:  06.5_monte_carlo_explainer_6panel_coral_reefs.R
# Description:  Builds a 2x3 (6-panel) publication-ready figure formatted for
#               Coral Reefs detailing the Monte Carlo age integration model:
#                 - Top Row (Panels a-c): Growth trajectory propagation
#                 - Bottom Row (Panels d-f): Resulting age & colonization year
#                   probability distributions at target sizes
# Dependencies: dplyr, ggplot2, patchwork, ragg, MASS, nlme
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(patchwork)
library(ragg)
library(MASS)   # mvrnorm()
library(nlme)   # coef()/vcov() S3 methods for the saved gnls model

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & MANUAL INPUT KNOBS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

model_path      <- file.path(base_dir, "models", "m_growth_gnls_SSasymp.rds")  # from Script 04
mc_inputs_path  <- file.path(base_dir, "models", "monte_carlo_inputs.rds")     # optional, from Script 06
output_fig_path <- file.path(base_dir, "Figures", "Fig_monte_carlo_explainer_CoralReefs.png")

FONT_FAMILY  <- "sans"
group_levels <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")

species_palette <- c(
  "Madrepora oculata"     = "#722082",
  "Desmophyllum pertusum" = "#B63679",
  "Primnoa msp."          = "#fb9f3a"
)

# ==============================================================================
# *** MANUAL SIZE ESTIMATION KNOBS ***
# Set target colony size(s) in cm^2 for age evaluation.
# ==============================================================================
TARGET_SIZES <- c(
  "Madrepora oculata"     = 57,  # cm^2
  "Desmophyllum pertusum" = 33,  # cm^2
  "Primnoa msp."          = 155  # cm^2
)

# --- Simulation Parameters ---------------------------------------------------
N_SIMS_SPAGHETTI <- 200     # Spaghetti trajectories plotted in Row 1
N_SIMS_FULL      <- 5000    # Draws behind Row 2's histograms
SURVEY_YEAR_LATE <- 2022    # Feeds secondary colonization-year axis
SEED             <- 42

set.seed(SEED)

if (length(TARGET_SIZES) == 1) {
  TARGET_SIZES <- setNames(rep(TARGET_SIZES, 3), group_levels)
}

# ------------------------------------------------------------------------------
# 3. LOAD FITTED MODEL & DETECTION FLOORS
# ------------------------------------------------------------------------------
if (!file.exists(model_path)) {
  stop("Could not find saved GNLS model object at '", model_path, "'.\nExecute Script 04 first.", call. = FALSE)
}
m_gnls_final <- readRDS(model_path)

param_means <- coef(m_gnls_final)
param_vcov  <- vcov(m_gnls_final)

A0_list       <- list()
max_size_list <- list()

if (file.exists(mc_inputs_path)) {
  cat("Loading Script 06 outputs for audited detection floors (A0)...\n")
  mc_inputs <- readRDS(mc_inputs_path)
  for (sp in group_levels) {
    A0_list[[sp]] <- unname(mc_inputs$A0_species[sp])
    obs_sizes     <- c(mc_inputs$stand_age_audited$A1[mc_inputs$stand_age_audited$Species == sp],
                       mc_inputs$stand_age_audited$A2[mc_inputs$stand_age_audited$Species == sp])
    max_size_list[[sp]] <- max(obs_sizes[!is.na(obs_sizes) & obs_sizes > 0], na.rm = TRUE)
  }
} else {
  cat("Script 06 outputs not found -- using illustrative defaults...\n")
  for (sp in group_levels) {
    A0_list[[sp]]       <- 5    # cm^2 illustrative floor
    max_size_list[[sp]] <- 250  # cm^2 illustrative upper limit
  }
}

# ------------------------------------------------------------------------------
# 4. MODEL INTEGRATION FUNCTIONS
# ------------------------------------------------------------------------------
get_draw_params <- function(draw_row, species_name) {
  if (!species_name %in% group_levels) return(NULL)

  asym <- draw_row["Asym.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = 0,
    "Desmophyllum pertusum" = draw_row["Asym.fSpeciesDesmophyllum pertusum"],
    "Primnoa msp."          = draw_row["Asym.fSpeciesPrimnoa msp."]
  )

  r0 <- draw_row["R0.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = 0,
    "Desmophyllum pertusum" = draw_row["R0.fSpeciesDesmophyllum pertusum"],
    "Primnoa msp."          = draw_row["R0.fSpeciesPrimnoa msp."]
  )

  lrc <- draw_row["lrc.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = 0,
    "Desmophyllum pertusum" = draw_row["lrc.fSpeciesDesmophyllum pertusum"],
    "Primnoa msp."          = draw_row["lrc.fSpeciesPrimnoa msp."]
  )

  list(asym = asym, r0 = r0, lrc = lrc)
}

calc_age_draw <- function(target_size, species_name, draw_row, A0) {
  if (is.na(target_size) || is.na(A0) || A0 <= 0) return(NA_real_)
  if (target_size <= A0) return(0)

  p <- get_draw_params(draw_row, species_name)
  if (is.null(p)) return(NA_real_)

  f_integrand <- function(s) {
    rgr <- p$asym + (p$r0 - p$asym) * exp(-exp(p$lrc) * s)
    if (any(rgr <= 0, na.rm = TRUE)) return(rep(NA_real_, length(s)))
    1 / (s * rgr)
  }

  tryCatch({
    val <- integrate(f_integrand, lower = A0, upper = target_size)$value
    if (is.nan(val) || val < 0) NA_real_ else val
  }, error = function(e) NA_real_)
}

# ------------------------------------------------------------------------------
# 5. GENERATE MONTE CARLO DRAWS
# ------------------------------------------------------------------------------
cat("Drawing parameter sets from multivariate normal distribution...\n")
draws_spaghetti <- mvrnorm(n = N_SIMS_SPAGHETTI, mu = param_means, Sigma = param_vcov)
draws_full      <- mvrnorm(n = N_SIMS_FULL,      mu = param_means, Sigma = param_vcov)

plots_row1 <- list()
plots_row2 <- list()

tag_chars_row1 <- c("a", "b", "c")
tag_chars_row2 <- c("d", "e", "f")

# ------------------------------------------------------------------------------
# 6. BUILD PANELS FOR EACH SPECIES
# ------------------------------------------------------------------------------
for (i in seq_along(group_levels)) {
  sp          <- group_levels[i]
  target_size <- TARGET_SIZES[[sp]]
  A0          <- A0_list[[sp]]
  max_sz      <- max_size_list[[sp]]
  sp_color    <- species_palette[[sp]]
  
  cat(sprintf("\nProcessing %s (Target Size = %.1f cm^2, A0 = %.1f cm^2)...\n", sp, target_size, A0))
  
  # --- Row 1 Data: Spaghetti Trajectories -------------------------------------
  size_grid  <- exp(seq(log(A0), log(max_sz), length.out = 60))
  
  curve_mean <- data.frame(
    Size = size_grid,
    Age  = sapply(size_grid, calc_age_draw, species_name = sp, draw_row = param_means, A0 = A0)
  )
  
  spaghetti_df <- bind_rows(lapply(seq_len(N_SIMS_SPAGHETTI), function(k) {
    data.frame(
      Draw = k,
      Size = size_grid,
      Age  = sapply(size_grid, calc_age_draw, species_name = sp, draw_row = draws_spaghetti[k, ], A0 = A0)
    )
  }))
  
  # --- Row 1 Plot: Trajectory & Uncertainty Propagation -----------------------
  p_top <- ggplot(spaghetti_df, aes(Size, Age)) +
    geom_line(aes(group = Draw), color = sp_color, alpha = 0.08, linewidth = 0.35) +
    geom_line(data = curve_mean, color = "black", linewidth = 0.7) +
    geom_vline(xintercept = target_size, linetype = "dashed", color = "grey20", linewidth = 0.4) +
    labs(
      tag      = tag_chars_row1[i],
      title    = bquote(bolditalic(.(sp))),
      subtitle = "Growth trajectory propagation",
      x        = expression(bold(paste("Colony size (cm"^2, ")"))),
      y        = expression(bold("Implied age (yr)"))
    ) +
    theme_classic(base_size = 8, base_family = FONT_FAMILY) +
    theme(
      plot.tag          = element_text(face = "bold", size = 9, vjust = 1),
      plot.title        = element_text(size = 8.5, hjust = 0),
      plot.subtitle     = element_text(color = "grey40", size = 7),
      axis.title        = element_text(size = 7.5, face = "bold", color = "black"),
      axis.text         = element_text(size = 7, color = "black"),
      axis.line         = element_line(linewidth = 0.3, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
    )
  
  # --- Row 2 Data: Full Age Distribution --------------------------------------
  ages_at_target <- sapply(seq_len(N_SIMS_FULL), function(k) {
    calc_age_draw(target_size, sp, draws_full[k, ], A0)
  })
  ages_at_target <- ages_at_target[!is.na(ages_at_target)]
  
  age_med <- median(ages_at_target)
  age_lo  <- unname(quantile(ages_at_target, 0.025))
  age_hi  <- unname(quantile(ages_at_target, 0.975))
  
  # --- Row 2 Plot: Age Distribution Histogram --------------------------------
  p_bottom <- ggplot(data.frame(Age = ages_at_target), aes(Age)) +
    geom_histogram(fill = sp_color, alpha = 0.75, bins = 35, color = "white", linewidth = 0.1) +
    geom_vline(xintercept = age_med, color = "black", linewidth = 0.6) +
    geom_vline(xintercept = c(age_lo, age_hi), color = "black", linetype = "dotted", linewidth = 0.5) +
    scale_x_continuous(
      sec.axis = sec_axis(~ SURVEY_YEAR_LATE - ., name = "Colonization year")
    ) +
    labs(
      tag      = tag_chars_row2[i],
      title    = sprintf("Age estimate at %.0f cm\u00b2", target_size),
      subtitle = sprintf("%.1f yr [%.1f\u2013%.1f]", age_med, age_lo, age_hi),
      x        = expression(bold("Implied age (yr)")),
      y        = expression(bold("Monte Carlo draws"))
    ) +
    theme_classic(base_size = 8, base_family = FONT_FAMILY) +
    theme(
      plot.tag          = element_text(face = "bold", size = 9, vjust = 1),
      plot.title        = element_text(face = "bold", size = 8),
      plot.subtitle     = element_text(color = "grey30", size = 7.5, face = "bold"),
      axis.title        = element_text(size = 7.5, face = "bold", color = "black"),
      axis.title.x.top  = element_text(size = 7, face = "bold", color = "black"),
      axis.text         = element_text(size = 7, color = "black"),
      axis.line         = element_line(linewidth = 0.3, color = "black"),
      axis.ticks        = element_line(linewidth = 0.3, color = "black"),
      plot.margin       = margin(t = 4, r = 4, b = 4, l = 4)
    )
  
  plots_row1[[sp]] <- p_top
  plots_row2[[sp]] <- p_bottom
}

# ------------------------------------------------------------------------------
# 7. ASSEMBLE COMPOSITE 2x3 LAYOUT & EXPORT FOR CORAL REEFS
# ------------------------------------------------------------------------------
cat("\nAssembling 6-panel composite layout...\n")

row1_layout <- plots_row1[[1]] | plots_row1[[2]] | plots_row1[[3]]
row2_layout <- plots_row2[[1]] | plots_row2[[2]] | plots_row2[[3]]

# Combined composite with plot layout spacing but NO overall main title
master_composite <- (row1_layout / row2_layout) & 
  theme(plot.margin = margin(t = 3, r = 3, b = 3, l = 3))

dir.create(dirname(output_fig_path), recursive = TRUE, showWarnings = FALSE)

cat("Exporting Coral Reefs submission graphic (174 mm width, 600 DPI PNG)...\n")
ggsave(
  filename = output_fig_path,
  plot     = master_composite,
  width    = 174,   # Standard Springer double-column width (mm)
  height   = 125,   # Proportional height for 2x3 panel layout (mm)
  units    = "mm",
  dpi      = 600,
  device   = ragg::agg_png,
  bg       = "white"
)

cat("\n======================================================================\n")
cat("Execution Complete!\n")
cat("Coral Reefs figure exported to:\n  ", output_fig_path, "\n")
cat("======================================================================\n")
