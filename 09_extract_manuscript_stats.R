# ==============================================================================
# Script Name:  09_extract_manuscript_stats.R
# Description:  Extracts exact medians, 95% CIs, sample sizes, and demographic 
#               proportions from upstream outputs to populate manuscript text.
# ==============================================================================

library(dplyr)
library(readr)
library(purrr)

base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

# ------------------------------------------------------------------------------
# 1. LOAD UPSTREAM OUTPUTS
# ------------------------------------------------------------------------------
final_master_path <- file.path(base_dir, "results", "Master_Dataset_Final_MonteCarlo_Results.csv")
mc_inputs_path    <- file.path(base_dir, "models", "monte_carlo_inputs.rds")
athresh_path      <- file.path(base_dir, "models", "athresh_sizes.rds")

if (!file.exists(final_master_path) || !file.exists(mc_inputs_path)) {
  stop("Missing required outputs from Script 06/07. Run previous scripts first.")
}

final_df   <- read_csv(final_master_path, show_col_types = FALSE)
mc_inputs  <- readRDS(mc_inputs_path)

param_draws   <- mc_inputs$param_draws
A0_species    <- mc_inputs$A0_species
valid_species <- mc_inputs$valid_species

# Load A_thresh sizes (or default to fitted script values if file not present)
athresh_sizes <- if (file.exists(athresh_path)) {
  readRDS(athresh_path)
} else {
  c("Desmophyllum pertusum" = 33.4, "Madrepora oculata" = 56.7, "Primnoa msp." = 155.4)
}

# ------------------------------------------------------------------------------
# 2. HELPER: MONTE CARLO AGE AT A_THRESH WITH 95% CI
# ------------------------------------------------------------------------------
get_draw_params <- function(draw_row, species_name) {
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

calc_athresh_age_draws <- function(species_name, target_size, param_draws, A0_vec) {
  A0 <- A0_vec[species_name]
  n_draws <- nrow(param_draws)
  
  ages <- sapply(seq_len(n_draws), function(k) {
    p <- get_draw_params(param_draws[k, ], species_name)
    f_integrand <- function(s) {
      rgr <- p$asym + (p$r0 - p$asym) * exp(-exp(p$lrc) * s)
      if (any(rgr <= 0, na.rm = TRUE)) return(rep(NA_real_, length(s)))
      1 / (s * rgr)
    }
    tryCatch(integrate(f_integrand, lower = A0, upper = target_size)$value, error = function(e) NA_real_)
  })
  
  data.frame(
    Species     = species_name,
    Athresh_cm2 = target_size,
    Median_Age  = round(median(ages, na.rm = TRUE), 1),
    Lower_95    = round(quantile(ages, 0.025, na.rm = TRUE), 2),
    Upper_95    = round(quantile(ages, 0.975, na.rm = TRUE), 2)
  )
}

# ------------------------------------------------------------------------------
# 3. PRINT STATISTICAL SUMMARY FOR MANUSCRIPT
# ------------------------------------------------------------------------------
cat("\n======================================================================\n")
cat("1. AGE AT A_THRESH TRANSITION (MEDIAN & 95% CI)\n")
cat("======================================================================\n")

athresh_summary <- map_dfr(names(athresh_sizes), function(sp) {
  calc_athresh_age_draws(sp, athresh_sizes[[sp]], param_draws, A0_species)
})
print(athresh_summary)

cat("\n======================================================================\n")
cat("2. 2022 POPULATION AGE DEMOGRAPHICS (FIGURE 4B / RIDGELINE)\n")
cat("======================================================================\n")

pop_2022_summary <- final_df %>%
  filter(!is.na(Age_2022_Median), Species %in% valid_species) %>%
  group_by(Species) %>%
  summarise(
    Total_Colonies_2022 = n(),
    Median_Age_2022     = round(median(Age_2022_Median, na.rm = TRUE), 1),
    IQR_Age_2022        = round(IQR(Age_2022_Median, na.rm = TRUE), 1),
    Min_Age             = round(min(Age_2022_Median, na.rm = TRUE), 1),
    Max_Age             = round(max(Age_2022_Median, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(pop_2022_summary)

cat("\nTotal Survey Population Size (2022):", sum(pop_2022_summary$Total_Colonies_2022), "\n")

cat("\n======================================================================\n")
cat("3. TRACKED GROWTH COLONY SAMPLE SIZES (TABLE 1 / SCRIPT 04)\n")
cat("======================================================================\n")

tracked_summary <- final_df %>%
  filter(tracked_both == TRUE, Species %in% valid_species) %>%
  count(Species, name = "Tracked_Growth_Pairs")
print(tracked_summary)

cat("\n======================================================================\n")
cat("4. MATURATION / THRESHOLD PROPORTIONS (% OF POPULATION > A_THRESH AGE)\n")
cat("======================================================================\n")

# ------------------------------------------------------------------------------
# 4. MATURATION / THRESHOLD PROPORTIONS (% OF POPULATION > A_THRESH AGE)
# ------------------------------------------------------------------------------

maturation_summary <- final_df %>%
  filter(!is.na(Age_2022_Median), Species %in% valid_species) %>%
  left_join(
    athresh_summary %>% dplyr::select(Species, Athresh_Age = Median_Age), 
    by = "Species"
  ) %>%
  group_by(Species) %>%
  summarise(
    Total_n             = n(),
    Athresh_Age_Yrs     = first(Athresh_Age),
    n_Above_Threshold   = sum(Age_2022_Median >= Athresh_Age),
    Pct_Above_Threshold = round(100 * mean(Age_2022_Median >= Athresh_Age), 1),
    .groups = "drop"
  )

print(maturation_summary)
cat("======================================================================\n")
