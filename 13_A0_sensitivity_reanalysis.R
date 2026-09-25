# ==============================================================================
# Script Name:  13_A0_sensitivity_reanalysis.R
# Description:  Re-runs the detection-floor (A0) sensitivity analysis described
#               in the manuscript's Sensitivity analysis section, at an explicit,
#               documented shift magnitude (resolves AR34: the original shift
#               magnitude did not survive the pasted manuscript text, and the
#               reviewer asked for a re-run at double it). A0_SHIFT_PCT below is
#               the one number that magnitude depends on -- change it here and
#               every downstream statistic updates. Currently set to a 200%
#               increase (A0_shifted = A0_original * 3) per instruction.
#
#               Two checks are reproduced, both against the audited/tracked
#               dataset and the SAME param_draws used in Scripts 06/07/09 -- so
#               the only thing that changes between "baseline" and "shifted"
#               below is A0 itself; parameter uncertainty draws are held fixed
#               for a clean, isolated sensitivity comparison:
#
#               Check 1 (elapsed-time consistency): for colonies tracked in both
#               2015 and 2022, does the model-implied elapsed time between the
#               two integrated ages (Delta_Age_Median = Age_2022 - Age_2015)
#               still track the TRUE 7.0-year survey interval once A0 is
#               shifted? Reproduces the "predicted elapsed time... closely
#               matched the true 7.0-year survey interval" statistic, at both
#               A0 values, per species.
#
#               Check 2 (adult vs. sub-threshold age shift): for every 2022
#               colony, split at each species' Athresh (Script 09) using that
#               colony's own observed area, how much does its point-estimate age
#               move between baseline and shifted A0? Reproduces the "altered
#               median adult ages (A > Athresh) by only X%... sub-threshold
#               individuals (A < Athresh)... median relative change of Y%,
#               corresponding to an absolute shift of Z years" statistics.
#
#               NOTE ON RUNTIME: Section 4 below recomputes a full N_SIMS-draw
#               Monte Carlo age for every colony at the shifted A0 -- the same
#               order of computation as Script 06's own run. Use N_SIMS_CHECK to
#               subsample param_draws for a faster pass while iterating; leave
#               it NULL for the final, publication-grade run.
#
# Dependencies: dplyr, readr, purrr, ggplot2 (same packages as Scripts 06/07/09)
# ==============================================================================

library(dplyr)
library(readr)
library(purrr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION & FILE PATHS  (same base_dir convention as Scripts 01-12)
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

mc_inputs_path     <- file.path(base_dir, "models", "monte_carlo_inputs.rds")                        # Script 06
final_results_csv  <- file.path(base_dir, "results", "Master_Dataset_Final_MonteCarlo_Results.csv")  # Script 06
athresh_path       <- file.path(base_dir, "models", "athresh_sizes.rds")                             # Script 09 (optional)

results_dir <- file.path(base_dir, "results")
figures_dir <- file.path(base_dir, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

elapsed_out_csv <- file.path(results_dir, "A0_sensitivity_elapsed_time.csv")
shift_out_csv   <- file.path(results_dir, "A0_sensitivity_age_shift.csv")
shift_out_fig   <- file.path(figures_dir, "FigureS_A0_Sensitivity.png")

# THE shift magnitude -- the one number the manuscript text was missing. 2.00 =
# a 200% increase, i.e. A0_shifted = A0_original * (1 + 2.00) = A0_original * 3.
# Change this single constant to test a different magnitude.
A0_SHIFT_PCT         <- 1.00
SURVEY_INTERVAL_YRS  <- 7.0    # true 2015-to-2022 gap; matches Script 06

# Matches Script 06's Monte Carlo draw count by default (NULL = use every draw
# already in param_draws). Set an integer (e.g. 500) to subsample while testing
# the pipeline -- see the runtime note above.
N_SIMS_CHECK <- NULL

# ------------------------------------------------------------------------------
# 2. LOAD SCRIPT 06 / 09 OUTPUTS
# ------------------------------------------------------------------------------
if (!file.exists(mc_inputs_path) || !file.exists(final_results_csv)) {
  stop("Missing outputs from Script 06. Run 06_monte_carlo_age_simulation.R first.", call. = FALSE)
}

mc_inputs         <- readRDS(mc_inputs_path)
stand_age_audited <- mc_inputs$stand_age_audited
param_draws       <- mc_inputs$param_draws
A0_species        <- mc_inputs$A0_species
valid_species     <- mc_inputs$valid_species

# TagLab.Genet.Id round-trips through Script 06's CSV export, and readr can
# silently re-infer it as numeric on read if every id in the file happens to be
# plain numeric (the same gotcha flagged in Script 12) -- force character on
# both sides before the join in Section 4 or it can fail silently.
final_master_dataset <- read_csv(final_results_csv, show_col_types = FALSE) %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id))

stand_age_audited <- stand_age_audited %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id))

athresh_sizes <- if (file.exists(athresh_path)) {
  readRDS(athresh_path)
} else {
  # Same fallback Script 09 uses if athresh_sizes.rds isn't present
  c("Desmophyllum pertusum" = 33.4, "Madrepora oculata" = 56.7, "Primnoa msp." = 155.4)
}

if (!is.null(N_SIMS_CHECK)) {
  param_draws <- param_draws[seq_len(min(N_SIMS_CHECK, nrow(param_draws))), , drop = FALSE]
}

A0_species_shifted <- A0_species * (1 + A0_SHIFT_PCT)

cat("======================================================================\n")
cat("A0 SENSITIVITY RE-ANALYSIS  (shift = +", A0_SHIFT_PCT * 100, "%)\n", sep = "")
cat("======================================================================\n")
cat("Baseline A0 (cm^2):\n"); print(round(A0_species, 3))
cat("Shifted A0 (cm^2):\n");  print(round(A0_species_shifted, 3))
cat("Draws per colony:", nrow(param_draws), "\n\n")

# ------------------------------------------------------------------------------
# 3. AGE-INTEGRATION HELPERS
#    Identical to Scripts 06/07/09, so results here are directly comparable to
#    the already-published baseline numbers.
# ------------------------------------------------------------------------------
get_draw_params <- function(draw_row, species_name) {
  if (!species_name %in% valid_species) return(NULL)
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

calc_age_draw <- function(target_size, species_name, draw_row, A0_vec) {
  if (is.na(target_size) || is.na(species_name) || !species_name %in% valid_species) return(NA_real_)
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

# Monte Carlo median age (across param_draws) for one colony at one A0 vector
calc_age_median <- function(target_size, species_name, draws, A0_vec) {
  if (is.na(target_size) || is.na(species_name) || !species_name %in% valid_species) return(NA_real_)
  ages <- sapply(seq_len(nrow(draws)), function(k) calc_age_draw(target_size, species_name, draws[k, ], A0_vec))
  median(ages, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# 4. RE-RUN 2015 AND 2022 AGES UNDER THE SHIFTED A0
#    Baseline numbers are read straight from Script 06's own output
#    (final_master_dataset) rather than recomputed, so Check 1/2 baselines are
#    guaranteed identical to the already-reported manuscript figures. Mirrors
#    Script 06's own for-loop + progress-bar structure.
# ------------------------------------------------------------------------------
n_rows <- nrow(stand_age_audited)
age_2015_shift <- age_2022_shift <- rep(NA_real_, n_rows)

cat("Recomputing 2015 and 2022 ages for", n_rows, "colonies under the shifted A0...\n")
pb <- txtProgressBar(min = 0, max = n_rows, style = 3)

for (i in seq_len(n_rows)) {
  sp <- stand_age_audited$Species[i]
  if (!is.na(sp) && sp %in% valid_species) {
    sz_2015 <- stand_age_audited$A1[i]
    sz_2022 <- stand_age_audited$A2[i]
    if (!is.na(sz_2015) && sz_2015 > 0) {
      age_2015_shift[i] <- calc_age_median(sz_2015, sp, param_draws, A0_species_shifted)
    }
    if (!is.na(sz_2022) && sz_2022 > 0) {
      age_2022_shift[i] <- calc_age_median(sz_2022, sp, param_draws, A0_species_shifted)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

shifted <- stand_age_audited %>%
  mutate(
    Age_2015_Median_shifted = age_2015_shift,
    Age_2022_Median_shifted = age_2022_shift
  )

comparison <- final_master_dataset %>%
  dplyr::select(TagLab.Genet.Id, Species, A1, A2, tracked_both,
                Age_2015_Median, Age_2022_Median, Delta_Age_Median) %>%
  left_join(
    shifted %>% dplyr::select(TagLab.Genet.Id, Age_2015_Median_shifted, Age_2022_Median_shifted),
    by = "TagLab.Genet.Id"
  ) %>%
  mutate(
    Delta_Age_Median_shifted = ifelse(
      tracked_both & !is.na(Age_2015_Median_shifted) & !is.na(Age_2022_Median_shifted) & A2 > A1,
      Age_2022_Median_shifted - Age_2015_Median_shifted, NA_real_
    )
  )

# ------------------------------------------------------------------------------
# 5. CHECK 1 -- MODEL-IMPLIED ELAPSED TIME vs. TRUE 7.0-YEAR SURVEY INTERVAL
# ------------------------------------------------------------------------------
elapsed_check <- comparison %>%
  filter(tracked_both, !is.na(Delta_Age_Median)) %>%
  group_by(Species) %>%
  summarise(
    n                         = n(),
    True_Interval_Yrs         = SURVEY_INTERVAL_YRS,
    Median_Delta_Age_Baseline = round(median(Delta_Age_Median, na.rm = TRUE), 1),
    Median_Delta_Age_Shifted  = round(median(Delta_Age_Median_shifted, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  mutate(
    Baseline_minus_True = round(Median_Delta_Age_Baseline - True_Interval_Yrs, 1),
    Shifted_minus_True  = round(Median_Delta_Age_Shifted  - True_Interval_Yrs, 1)
  )

cat("\n======================================================================\n")
cat("CHECK 1: MODEL-IMPLIED ELAPSED TIME vs. TRUE", SURVEY_INTERVAL_YRS, "-YEAR INTERVAL\n")
cat("======================================================================\n")
print(elapsed_check)

write_csv(elapsed_check, elapsed_out_csv)
cat("Saved to:\n  ", elapsed_out_csv, "\n")

# ------------------------------------------------------------------------------
# 6. CHECK 2 -- ADULT (A > Athresh) vs. SUB-THRESHOLD (A < Athresh) AGE SHIFT
#    Uses each colony's own 2022 area against ITS species' Athresh (Script 09).
# ------------------------------------------------------------------------------
shift_check_colonies <- comparison %>%
  filter(Species %in% valid_species, !is.na(A2), A2 > 0,
         !is.na(Age_2022_Median), !is.na(Age_2022_Median_shifted)) %>%
  mutate(
    Athresh_cm2   = athresh_sizes[Species],
    Size_Class    = ifelse(A2 > Athresh_cm2, "Adult (A > Athresh)", "Sub-threshold (A < Athresh)"),
    Abs_Shift_Yrs = Age_2022_Median_shifted - Age_2022_Median,
    # Guard against division by zero for colonies whose baseline age is
    # exactly 0 (i.e. A2 <= baseline A0 -- possible near the 3rd-percentile floor)
    Pct_Shift     = ifelse(Age_2022_Median > 0, 100 * Abs_Shift_Yrs / Age_2022_Median, NA_real_)
  )

shift_summary <- shift_check_colonies %>%
  group_by(Size_Class) %>%
  summarise(
    n                    = n(),
    Median_Pct_Shift     = round(median(Pct_Shift, na.rm = TRUE), 2),
    Mean_Pct_Shift       = round(mean(Pct_Shift, na.rm = TRUE), 2),
    Median_Abs_Shift_Yrs = round(median(Abs_Shift_Yrs, na.rm = TRUE), 2),
    .groups = "drop"
  )

shift_summary_by_species <- shift_check_colonies %>%
  group_by(Species, Size_Class) %>%
  summarise(
    n                    = n(),
    Median_Pct_Shift     = round(median(Pct_Shift, na.rm = TRUE), 2),
    Median_Abs_Shift_Yrs = round(median(Abs_Shift_Yrs, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("\n======================================================================\n")
cat("CHECK 2: AGE SHIFT UNDER A0 +", A0_SHIFT_PCT * 100, "% -- ADULT vs. SUB-THRESHOLD\n", sep = "")
cat("======================================================================\n")
cat("Pooled across species:\n")
print(shift_summary)
cat("\nBy species:\n")
print(shift_summary_by_species, n = Inf)

write_csv(
  shift_check_colonies %>%
    dplyr::select(TagLab.Genet.Id, Species, A2, Athresh_cm2, Size_Class,
                  Age_2022_Median, Age_2022_Median_shifted, Abs_Shift_Yrs, Pct_Shift),
  shift_out_csv
)
cat("\nSaved per-colony results to:\n  ", shift_out_csv, "\n")

# ------------------------------------------------------------------------------
# 7. DIAGNOSTIC FIGURE -- BASELINE vs. SHIFTED AGE, PER COLONY
# ------------------------------------------------------------------------------
species_palette <- c(
  "Madrepora oculata"     = "#722082",
  "Desmophyllum pertusum" = "#B63679",
  "Primnoa msp."          = "#fb9f3a"
)

fig_df <- shift_check_colonies %>%
  mutate(Species = factor(Species, levels = names(species_palette)))

p_shift <- ggplot(fig_df, aes(x = Age_2022_Median, y = Age_2022_Median_shifted, color = Species)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.5, size = 1.4) +
  facet_wrap(~ Size_Class) +
  scale_color_manual(values = species_palette, name = "Species") +
  labs(
    title = paste0("Colony age under baseline vs. A0 +", A0_SHIFT_PCT * 100, "% detection floor"),
    x = "Baseline A0 -- 2022 median age (yr)",
    y = "Shifted A0 -- 2022 median age (yr)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text        = element_text(face = "bold"),
    legend.position   = "bottom"
  )

ggsave(shift_out_fig, p_shift, width = 178, height = 100, units = "mm", dpi = 600)
cat("\nSaved diagnostic figure to:\n  ", shift_out_fig, "\n")

cat("======================================================================\n")
cat("Script 13 execution complete!\n")
cat("======================================================================\n")
