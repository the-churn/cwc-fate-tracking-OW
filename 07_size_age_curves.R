# ==============================================================================
# Script Name:  07_size_age_curves.R
# Description:  Builds per-species size-at-age curves (median + 95% Monte Carlo
#               CI) by inverting the GNLS growth model across a log-spaced size
#               grid, using the parameter draws and audited dataset produced by
#               script 06. Also assembles the observed colony-level age-size
#               points used as an overlay in the script 08 figure.
# Dependencies: dplyr, readr
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
mc_inputs_path       <- file.path("outputs", "models", "monte_carlo_inputs.rds")
final_results_csv    <- file.path("outputs", "results", "Master_Dataset_Final_MonteCarlo_Results.csv")
curves_out_csv        <- file.path("outputs", "results", "size_age_curves.csv")
observed_out_csv       <- file.path("outputs", "results", "observed_points.csv")

dir.create(file.path("outputs", "results"), recursive = TRUE, showWarnings = FALSE)

curve_n_sims <- 5000   # draws used per grid point; keep in sync with N_SIMS in script 06
curve_grid_n <- 40     # points per species along the log-spaced size grid

# ------------------------------------------------------------------------------
# 3. LOAD UPSTREAM OUTPUTS (SCRIPT 06)
# ------------------------------------------------------------------------------
if (!file.exists(mc_inputs_path) || !file.exists(final_results_csv)) {
  stop("Missing outputs from script 06. Run 06_monte_carlo_age_simulation.R first.", call. = FALSE)
}

mc_inputs         <- readRDS(mc_inputs_path)
stand_age_audited <- mc_inputs$stand_age_audited
param_draws       <- mc_inputs$param_draws
A0_species        <- mc_inputs$A0_species
valid_species      <- mc_inputs$valid_species

final_master_dataset <- read_csv(final_results_csv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 4. AGE-INTEGRATION HELPERS
# ------------------------------------------------------------------------------
# Duplicated from script 06 so this script can run standalone from saved
# outputs (no need to re-fit or re-source anything). Keep in sync with script
# 06 if the growth-model form or parameter names ever change.
get_draw_params <- function(draw_row, species_name) {
  if (!species_name %in% valid_species) return(NULL)

  asym <- draw_row["Asym.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = draw_row["Asym.fSpeciesMadrepora oculata"],
    "Primnoa msp."          = draw_row["Asym.fSpeciesPrimnoa msp."],
    "Desmophyllum pertusum" = 0)

  r0 <- draw_row["R0.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = draw_row["R0.fSpeciesMadrepora oculata"],
    "Primnoa msp."          = draw_row["R0.fSpeciesPrimnoa msp."],
    "Desmophyllum pertusum" = 0)

  lrc <- draw_row["lrc.(Intercept)"] + switch(species_name,
    "Madrepora oculata"     = draw_row["lrc.fSpeciesMadrepora oculata"],
    "Primnoa msp."          = draw_row["lrc.fSpeciesPrimnoa msp."],
    "Desmophyllum pertusum" = 0)

  list(asym = asym, r0 = r0, lrc = lrc)
}

calc_age_draw <- function(target_size, species_name, draw_row, A0_vec) {
  if (is.na(target_size) || is.na(species_name) || !species_name %in% valid_species) {
    return(NA_real_)
  }
  A0 <- A0_vec[species_name]
  if (is.na(A0) || A0 <= 0) return(NA_real_)
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
# 5. BUILD SIZE-AGE CURVES (MEDIAN + 95% CI PER SPECIES)
# ------------------------------------------------------------------------------
build_size_age_curve <- function(species_name, df, param_draws, A0_vec, n_grid = 40) {
  n_sims_curve <- nrow(param_draws)
  A0 <- A0_vec[species_name]

  obs_sizes <- c(df$A1[df$Species == species_name], df$A2[df$Species == species_name])
  obs_sizes <- obs_sizes[!is.na(obs_sizes) & obs_sizes > 0]
  max_size  <- max(obs_sizes, na.rm = TRUE)

  size_grid <- exp(seq(log(A0), log(max_size), length.out = n_grid))

  grid_results <- lapply(size_grid, function(sz) {
    ages <- sapply(seq_len(n_sims_curve), function(k) {
      calc_age_draw(sz, species_name, param_draws[k, ], A0_vec)
    })
    data.frame(
      Species  = species_name,
      Size     = sz,
      Age_Med  = median(ages, na.rm = TRUE),
      Age_Lo95 = quantile(ages, 0.025, na.rm = TRUE),
      Age_Hi95 = quantile(ages, 0.975, na.rm = TRUE)
    )
  })
  bind_rows(grid_results)
}

cat("Building size-age curves (", curve_n_sims, "draws x", curve_grid_n, "grid points per species)...\n")

size_age_curves <- bind_rows(lapply(valid_species, function(sp) {
  build_size_age_curve(sp, stand_age_audited, param_draws[1:curve_n_sims, ], A0_species, n_grid = curve_grid_n)
}))

write_csv(size_age_curves, curves_out_csv)
cat("Saved size-age curves to:\n  ", curves_out_csv, "\n")

# ------------------------------------------------------------------------------
# 6. OBSERVED (COLONY-LEVEL) AGE-SIZE POINTS
# ------------------------------------------------------------------------------
observed_points <- bind_rows(
  final_master_dataset %>%
    filter(!is.na(Age_2015_Median)) %>%
    transmute(Species, Size = A1, Age = Age_2015_Median),
  final_master_dataset %>%
    filter(!is.na(Age_2022_Median)) %>%
    transmute(Species, Size = A2, Age = Age_2022_Median)
) %>%
  filter(Species %in% valid_species, !is.na(Size), Size > 0)

write_csv(observed_points, observed_out_csv)
cat("Saved observed age-size points to:\n  ", observed_out_csv, "\n")
