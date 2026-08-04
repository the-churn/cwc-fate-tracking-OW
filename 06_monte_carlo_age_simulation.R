# ==============================================================================
# Script Name:  06_monte_carlo_age_simulation.R
# Description:  Draws parameter sets from the fitted GNLS growth model (script
#               05) and Monte Carlo-integrates colony age from colony size, for
#               the 2015 and 2022 surveys. Includes the recruit-audit step
#               (candidates flagged automatically, reviewed manually in
#               TagLab, merged back in) and produces the final per-colony,
#               uncertainty-bounded age dataset used by scripts 07-08.
# Dependencies: openxlsx, dplyr, writexl, readxl, readr, MASS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(openxlsx)
library(dplyr)
library(writexl)
library(readxl)
library(readr)
library(MASS)   # mvrnorm()

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
raw_data_path        <- file.path("data", "Error_Propagation", "Master_CWC_Tracked_MDC_Age.xlsx")
master_csv_path       <- file.path("data", "master_data_combined.csv")
audit_template_path   <- file.path("data", "audit", "TagLab_QuickJump_Audit.xlsx")
model_path             <- file.path("outputs", "models", "m_gnls_final.rds")
mc_inputs_path         <- file.path("outputs", "models", "monte_carlo_inputs.rds")
final_results_csv     <- file.path("outputs", "results", "Master_Dataset_Final_MonteCarlo_Results.csv")

dir.create(file.path("data", "audit"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("outputs", "models"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("outputs", "results"), recursive = TRUE, showWarnings = FALSE)

valid_species          <- c("Desmophyllum pertusum", "Madrepora oculata", "Primnoa msp.")
N_SIMS                 <- 5000   # Monte Carlo draws from the GNLS parameter covariance
P_FLOOR                <- 0.03   # detection-floor percentile used to set species-specific A0
AGE_OUTLIER_THRESHOLD   <- 10     # years; flags 2022-only "recruits" for manual TagLab audit
SURVEY_YEAR_EARLY       <- 2015
SURVEY_YEAR_LATE        <- 2022
SURVEY_INTERVAL_YRS     <- SURVEY_YEAR_LATE - SURVEY_YEAR_EARLY
SEED                    <- 42

set.seed(SEED)

# ------------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# Standardize species labels. Run once, immediately after loading data.
clean_species_names <- function(df) {
  df %>%
    mutate(
      Species = trimws(Species),
      Species = case_when(
        Species %in% c("Primnoa msp.1", "Primnoa msp.5", "Primnoa msp") ~ "Primnoa msp.",
        Species %in% c("D. pertusum", "Desmophyllum pertusum")         ~ "Desmophyllum pertusum",
        Species == "Madrepora oculata"                                  ~ "Madrepora oculata",
        TRUE ~ Species  # keeps "Coral Recruit" or other labels intact
      )
    )
}

# Flag colonies tracked across both survey years + growth direction.
# Call again any time A1/A2 change (e.g. after the audit merge).
flag_tracked <- function(df) {
  df %>%
    mutate(
      tracked_both = !is.na(A1) & !is.na(A2) & A1 > 0 & A2 > 0,
      delta_area   = ifelse(tracked_both, A2 - A1, NA_real_),
      is_growing   = ifelse(tracked_both, delta_area > 0, NA)
    )
}

# Species-specific detection floor (A0): the size below which a colony is
# assumed undetectable / age = 0.
compute_A0 <- function(df, valid_species, p_floor) {
  df_long <- data.frame(
    Species = c(df$Species, df$Species),
    Area    = c(df$A1, df$A2)
  )
  A0_summary <- df_long %>%
    filter(!is.na(Species), Species %in% valid_species, !is.na(Area), Area > 0) %>%
    group_by(Species) %>%
    summarise(min_size = quantile(Area, probs = p_floor, na.rm = TRUE), .groups = "drop")
  setNames(A0_summary$min_size, A0_summary$Species)
}

# Extract species-specific (Asym, R0, lrc) for one Monte Carlo draw.
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

# Core numerical integration: age implied by a given size, for one draw.
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

# Runs the full Monte Carlo age simulation over every row of df. n_sims is
# inferred from nrow(param_draws) -- pass a subset of rows for a cheaper pass.
simulate_ages <- function(df, param_draws, A0_vec) {
  n_sims <- nrow(param_draws)
  n_rows <- nrow(df)

  age_2015_med <- age_2015_lo <- age_2015_hi <- rep(NA_real_, n_rows)
  age_2022_med <- age_2022_lo <- age_2022_hi <- rep(NA_real_, n_rows)
  delta_t_med  <- delta_t_lo  <- delta_t_hi  <- rep(NA_real_, n_rows)

  cat("  Running Monte Carlo simulation across", n_rows, "rows (", n_sims, "draws)...\n")
  pb <- txtProgressBar(min = 0, max = n_rows, style = 3)

  for (i in seq_len(n_rows)) {
    sp <- df$Species[i]
    if (!is.na(sp) && sp %in% valid_species) {
      sz_2015 <- df$A1[i]
      sz_2022 <- df$A2[i]

      sim_2015 <- if (!is.na(sz_2015) && sz_2015 > 0) {
        sapply(seq_len(n_sims), function(k) calc_age_draw(sz_2015, sp, param_draws[k, ], A0_vec))
      } else NULL

      sim_2022 <- if (!is.na(sz_2022) && sz_2022 > 0) {
        sapply(seq_len(n_sims), function(k) calc_age_draw(sz_2022, sp, param_draws[k, ], A0_vec))
      } else NULL

      if (!is.null(sim_2015)) {
        age_2015_med[i] <- median(sim_2015, na.rm = TRUE)
        age_2015_lo[i]  <- quantile(sim_2015, probs = 0.025, na.rm = TRUE)
        age_2015_hi[i]  <- quantile(sim_2015, probs = 0.975, na.rm = TRUE)
      }
      if (!is.null(sim_2022)) {
        age_2022_med[i] <- median(sim_2022, na.rm = TRUE)
        age_2022_lo[i]  <- quantile(sim_2022, probs = 0.025, na.rm = TRUE)
        age_2022_hi[i]  <- quantile(sim_2022, probs = 0.975, na.rm = TRUE)
      }
      if (!is.null(sim_2015) && !is.null(sim_2022) && sz_2022 > sz_2015) {
        d <- sim_2022 - sim_2015
        delta_t_med[i] <- median(d, na.rm = TRUE)
        delta_t_lo[i]  <- quantile(d, probs = 0.025, na.rm = TRUE)
        delta_t_hi[i]  <- quantile(d, probs = 0.975, na.rm = TRUE)
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  df %>%
    mutate(
      Age_2015_Median       = age_2015_med,
      Age_2015_Lower95      = age_2015_lo,
      Age_2015_Upper95      = age_2015_hi,
      Age_2022_Median       = age_2022_med,
      Age_2022_Lower95      = age_2022_lo,
      Age_2022_Upper95      = age_2022_hi,
      Delta_Age_Median      = delta_t_med,
      Delta_Age_Lower95     = delta_t_lo,
      Delta_Age_Upper95     = delta_t_hi,
      Colonization_Year_Est = ifelse(!is.na(age_2022_med), SURVEY_YEAR_LATE - age_2022_med, NA_real_)
    )
}

# ------------------------------------------------------------------------------
# 4. LOAD FITTED MODEL (SCRIPT 05) & PREPARE RAW DATA
# ------------------------------------------------------------------------------
if (!file.exists(model_path)) {
  stop(
    "Could not find a saved model at '", model_path, "'.\n",
    "Run 05_growth_rate_publication_figure.R first -- it fits m_gnls_final and\n",
    "should save it there via saveRDS(). See the one-line addition for script 05.",
    call. = FALSE
  )
}
m_gnls_final <- readRDS(model_path)

cat("Loading and preparing raw tracked-colony dataset...\n")
stand_age_raw <- read.xlsx(raw_data_path) %>%
  clean_species_names() %>%
  flag_tracked()

cat("Total rows:", nrow(stand_age_raw), "\n")
cat("Colonies tracked", SURVEY_YEAR_EARLY, "&", SURVEY_YEAR_LATE, ":",
    sum(stand_age_raw$tracked_both, na.rm = TRUE), "\n")

A0_species <- compute_A0(stand_age_raw, valid_species, P_FLOOR)
cat("Species-specific detection floor (A0 at", P_FLOOR * 100, "th percentile):\n")
print(round(A0_species, 3))

# ------------------------------------------------------------------------------
# 5. MONTE CARLO PARAMETER DRAWS (FROM GNLS MODEL)
# ------------------------------------------------------------------------------
param_means <- coef(m_gnls_final)
param_vcov  <- vcov(m_gnls_final)

# All N_SIMS draws generated once and reused everywhere downstream, so the
# exploratory pass and the final pass come from the same distribution.
param_draws <- mvrnorm(n = N_SIMS, mu = param_means, Sigma = param_vcov)

# ------------------------------------------------------------------------------
# 6. RECRUIT AUDIT: BUILD TEMPLATE IF NEEDED, OTHERWISE SKIP
# ------------------------------------------------------------------------------
# Colonies present only in 2022 with an implausibly old estimated age are
# candidates for having been missed/occluded in the 2015 survey rather than
# true recruits. This step only needs to run once -- if an audit file already
# exists (i.e. you've already reviewed it in TagLab), it's skipped entirely.
if (!file.exists(audit_template_path)) {

  cat("\nNo audit file found -- running exploratory pass to flag recruit outliers...\n")
  n_sims_explore <- min(500, N_SIMS)
  stand_age_explore <- simulate_ages(stand_age_raw, param_draws[1:n_sims_explore, ], A0_species)

  outlier_recruits <- stand_age_explore %>%
    filter((is.na(A1) | A1 == 0) & !is.na(A2) & A2 > 0) %>%
    filter(!is.na(Age_2022_Median), Age_2022_Median >= AGE_OUTLIER_THRESHOLD) %>%
    arrange(desc(Age_2022_Median))

  cat("Flagged", nrow(outlier_recruits), "recruit outliers (age >=", AGE_OUTLIER_THRESHOLD, "yrs):\n")
  print(table(outlier_recruits$Species))

  master_df <- read_csv(master_csv_path, show_col_types = FALSE)
  id_map_2022 <- master_df %>%
    filter(Year == SURVEY_YEAR_LATE) %>%
    mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
    distinct(TagLab.Genet.Id, .keep_all = TRUE) %>%
    dplyr::select(TagLab.Genet.Id, TagLab.Id_2022 = TagLab.Id)

  audit_template <- outlier_recruits %>%
    mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
    left_join(id_map_2022, by = "TagLab.Genet.Id") %>%
    transmute(
      TagLab.Genet.Id,
      TagLab.Id_2022,                # search this ID directly in TagLab
      Species,
      A2_Area_2022 = A2,
      Est_Age_2022 = round(Age_2022_Median, 1),
      Audit_Status = NA_character_,  # fill: "confirmed_recruit" / "occluded" / "artifact" / "growth"
      Match_2015   = NA_character_,  # if "growth": the matching 2015 TagLab.Id
      A1_2015      = NA_real_        # if "growth": the matched 2015 area
    )

  write_xlsx(audit_template, audit_template_path)
  cat("\nAudit template written to:\n  ", audit_template_path, "\n")
  cat("Open it, look up each TagLab.Id_2022 in TagLab, and fill in Audit_Status\n")
  cat("(and Match_2015 / A1_2015 where relevant). Save the file, then re-run this script.\n")
  stop("Stopping here for the manual audit step -- see instructions above.", call. = FALSE)
}

cat("\nExisting audit file found at:\n  ", audit_template_path, "\n  Merging results and continuing.\n")

# ------------------------------------------------------------------------------
# 7. MERGE AUDIT RESULTS -> BUILD AUDITED DATASET
# ------------------------------------------------------------------------------
audit_results <- read_excel(audit_template_path) %>%
  mutate(
    TagLab.Genet.Id = as.character(TagLab.Genet.Id),
    A1_2015         = as.numeric(A1_2015)
  )

stand_age_audited <- stand_age_raw %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
  left_join(
    audit_results %>% dplyr::select(TagLab.Genet.Id, Audit_Status, Match_2015, A1_2015),
    by = "TagLab.Genet.Id"
  ) %>%
  mutate(
    # "growth" colonies get their true 2015 area restored -- only if a
    # confirmed match was actually recorded (an unconfirmed guess must NOT
    # be treated as evidence).
    A1 = case_when(
      !is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015) ~ A1_2015,
      TRUE ~ A1
    ),
    Known_Non_Recruit = (!is.na(Audit_Status) & Audit_Status == "occluded") |
                        (!is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015))
  ) %>%
  filter(is.na(Audit_Status) | Audit_Status != "artifact") %>%
  flag_tracked()  # recompute now that A1 may have changed

cat("Audited dataset ready:", nrow(stand_age_audited), "rows (",
    sum(!is.na(stand_age_audited$Audit_Status)), "reviewed).\n")

# ------------------------------------------------------------------------------
# 8. FINAL MONTE CARLO SIMULATION & EXPORT
# ------------------------------------------------------------------------------
cat("\nRunning final simulation (N_SIMS =", N_SIMS, ")...\n")
final_master_dataset <- simulate_ages(stand_age_audited, param_draws, A0_species) %>%
  mutate(
    Annual_Growth_cm2_yr = ifelse(!is.na(A1) & !is.na(A2) & A1 > 0,
                                   round((A2 - A1) / SURVEY_INTERVAL_YRS, 2), NA_real_)
  )

cat("\nSimulation complete -- species summary:\n")
print(
  final_master_dataset %>%
    group_by(Species) %>%
    summarise(
      Total_Colonies     = n(),
      Tracked_Colonies   = sum(tracked_both, na.rm = TRUE),
      Median_Age_2022    = round(median(Age_2022_Median, na.rm = TRUE), 1),
      Mean_Growth_cm2_yr = round(mean(Annual_Growth_cm2_yr, na.rm = TRUE), 2),
      .groups = "drop"
    )
)

write_csv(final_master_dataset, final_results_csv)
cat("\nSaved final Monte Carlo age dataset to:\n  ", final_results_csv, "\n")

# ------------------------------------------------------------------------------
# 9. SAVE DOWNSTREAM INPUTS FOR SCRIPTS 07-08
# ------------------------------------------------------------------------------
saveRDS(
  list(
    stand_age_audited = stand_age_audited,
    param_draws       = param_draws,
    A0_species        = A0_species,
    valid_species      = valid_species
  ),
  mc_inputs_path
)
cat("Saved Monte Carlo inputs (for scripts 07-08) to:\n  ", mc_inputs_path, "\n")
