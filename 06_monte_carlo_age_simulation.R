# ==============================================================================
# Script Name:  06_monte_carlo_age_simulation.R
# Description:  Draws parameter sets from the fitted SSasymp GNLS growth model
#               (Script 04/05) and Monte Carlo-integrates colony age from colony 
#               planar size for 2015 and 2022 surveys. Handles recruit-audit 
#               workflow (TagLab integration) and exports uncertainty-bounded 
#               age dataset for scripts 07-08.
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
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

# Ensure subdirectories exist under base_dir
dir.create(file.path(base_dir, "models"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "results"), recursive = TRUE, showWarnings = FALSE)

# Input Paths
raw_data_path       <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")
master_csv_path     <- file.path(base_dir, "master_data_combined.csv")
model_path          <- file.path(base_dir, "models", "m_growth_gnls_SSasymp.rds")

# Audit Path (Saved directly in base_dir so a fresh file is created)
audit_template_path <- file.path(base_dir, "TagLab_QuickJump_Audit.xlsx")

# Output Paths
mc_inputs_path      <- file.path(base_dir, "models", "monte_carlo_inputs.rds")
final_results_csv   <- file.path(base_dir, "results", "Master_Dataset_Final_MonteCarlo_Results.csv")

# Hyperparameters
valid_species          <- c("Madrepora oculata", "Desmophyllum pertusum", "Primnoa msp.")
N_SIMS                 <- 5000   # Monte Carlo draws from parameter covariance matrix
P_FLOOR                <- 0.03   # Detection-floor percentile (sets species-specific A0)
AGE_OUTLIER_THRESHOLD  <- 10     # Years; flags 2022-only "recruits" for manual TagLab audit
SURVEY_YEAR_EARLY      <- 2015
SURVEY_YEAR_LATE       <- 2022
SURVEY_INTERVAL_YRS    <- SURVEY_YEAR_LATE - SURVEY_YEAR_EARLY
SEED                   <- 42

set.seed(SEED)

# ------------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# Standardize species labels across input datasets
clean_species_names <- function(df) {
  df %>%
    mutate(
      Species = trimws(Species),
      Species = case_when(
        Species %in% c("Primnoa msp.1", "Primnoa msp.5", "Primnoa msp.") ~ "Primnoa msp.",
        Species %in% c("D. pertusum", "Desmophyllum pertusum")         ~ "Desmophyllum pertusum",
        Species == "Madrepora oculata"                                 ~ "Madrepora oculata",
        TRUE ~ Species
      )
    )
}

# Flag tracking status and area changes
flag_tracked <- function(df) {
  df %>%
    mutate(
      tracked_both = !is.na(A1) & !is.na(A2) & A1 > 0 & A2 > 0,
      delta_area   = ifelse(tracked_both, A2 - A1, NA_real_),
      is_growing   = ifelse(tracked_both, delta_area > 0, NA)
    )
}

# Species-specific detection floor (A0)
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

# Extract parameter set for one Monte Carlo draw
# Matches fitted GNLS coefficients (Madrepora oculata = Intercept Baseline)
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

# Core numerical integration: age implied by size for one draw
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

# Monte Carlo integration runner
simulate_ages <- function(df, param_draws, A0_vec) {
  n_sims <- nrow(param_draws)
  n_rows <- nrow(df)

  age_2015_med <- age_2015_lo <- age_2015_hi <- rep(NA_real_, n_rows)
  age_2022_med <- age_2022_lo <- age_2022_hi <- rep(NA_real_, n_rows)
  delta_t_med  <- delta_t_lo  <- delta_t_hi  <- rep(NA_real_, n_rows)

  cat("Running Monte Carlo simulation across", n_rows, "rows (", n_sims, "draws per row)...\n")
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
# 4. LOAD FITTED MODEL & PREPARE RAW DATA
# ------------------------------------------------------------------------------
if (!file.exists(model_path)) {
  stop("Could not find saved GNLS model object at '", model_path, "'.\nExecute Script 04/05 first.", call. = FALSE)
}
m_gnls_final <- readRDS(model_path)

cat("Loading raw tracked colony dataset from:", raw_data_path, "\n")
stand_age_raw <- read.xlsx(raw_data_path) %>%
  clean_species_names() %>%
  flag_tracked()

cat("Total rows:", nrow(stand_age_raw), "\n")
cat("Colonies tracked in both years:", sum(stand_age_raw$tracked_both, na.rm = TRUE), "\n")

A0_species <- compute_A0(stand_age_raw, valid_species, P_FLOOR)
cat("Species-specific detection floor (A0 at", P_FLOOR * 100, "th percentile):\n")
print(round(A0_species, 3))

# ------------------------------------------------------------------------------
# 5. MONTE CARLO PARAMETER DRAWS
# ------------------------------------------------------------------------------
param_means <- coef(m_gnls_final)
param_vcov  <- vcov(m_gnls_final)

param_draws <- mvrnorm(n = N_SIMS, mu = param_means, Sigma = param_vcov)

# ------------------------------------------------------------------------------
# 6. RECRUIT AUDIT WORKFLOW
# ------------------------------------------------------------------------------
if (!file.exists(audit_template_path)) {

  cat("\nNo existing audit file found. Running exploratory pass to flag recruits...\n")
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
      TagLab.Id_2022,
      Species,
      A2_Area_2022 = A2,
      Est_Age_2022 = round(Age_2022_Median, 1),
      Audit_Status = NA_character_,
      Match_2015   = NA_character_,
      A1_2015      = NA_real_
    )

  write_xlsx(audit_template, audit_template_path)
  cat("\nAudit template saved to:\n  ", audit_template_path, "\n")
  stop("Execution paused for manual TagLab review.", call. = FALSE)
}

cat("\nFound completed audit template at:\n  ", audit_template_path, "\nMerging audit decisions...\n")

# ------------------------------------------------------------------------------
# 7. MERGE AUDIT DECISIONS
# ------------------------------------------------------------------------------
audit_results <- read_excel(audit_template_path) %>%
  mutate(
    TagLab.Genet.Id = as.character(TagLab.Genet.Id),
    A1_2015         = as.numeric(A1_2015)
  )

n_before_audit_merge <- nrow(stand_age_raw)

stand_age_audited <- stand_age_raw %>%
  mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
  left_join(
    audit_results %>% dplyr::select(TagLab.Genet.Id, Audit_Status, Match_2015, A1_2015),
    by = "TagLab.Genet.Id"
  )

if (nrow(stand_age_audited) != n_before_audit_merge) {
  dup_ids <- stand_age_raw %>%
    mutate(TagLab.Genet.Id = as.character(TagLab.Genet.Id)) %>%
    count(TagLab.Genet.Id) %>%
    filter(n > 1) %>%
    pull(TagLab.Genet.Id)
  stop(
    "Audit merge duplicate error: TagLab.Genet.Id non-unique on join.\n",
    "Duplicated ID(s): ", paste(dup_ids, collapse = ", "),
    call. = FALSE
  )
}

stand_age_audited <- stand_age_audited %>%
  mutate(
    A1 = case_when(
      !is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015) ~ A1_2015,
      TRUE ~ A1
    ),
    Known_Non_Recruit = (!is.na(Audit_Status) & Audit_Status == "occluded") |
                        (!is.na(Audit_Status) & Audit_Status == "growth" & !is.na(A1_2015))
  ) %>%
  filter(is.na(Audit_Status) | Audit_Status != "artifact") %>%
  flag_tracked()

# ------------------------------------------------------------------------------
# 8. FINAL MONTE CARLO SIMULATION & EXPORT
# ------------------------------------------------------------------------------
cat("\nRunning final Monte Carlo simulation (N =", N_SIMS, "draws)...\n")
final_master_dataset <- simulate_ages(stand_age_audited, param_draws, A0_species) %>%
  mutate(
    Annual_Growth_cm2_yr = ifelse(!is.na(A1) & !is.na(A2) & A1 > 0,
                                  round((A2 - A1) / SURVEY_INTERVAL_YRS, 2), NA_real_)
  )

cat("\nSpecies-level summary:\n")
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
cat("\nSaved final Monte Carlo results to:\n  ", final_results_csv, "\n")

# ------------------------------------------------------------------------------
# 9. SAVE DOWNSTREAM INPUTS FOR SCRIPTS 07-08
# ------------------------------------------------------------------------------
saveRDS(
  list(
    stand_age_audited = stand_age_audited,
    param_draws        = param_draws,
    A0_species         = A0_species,
    valid_species      = valid_species
  ),
  mc_inputs_path
)
cat("Saved Monte Carlo workspace inputs for Scripts 07-08 to:\n  ", mc_inputs_path, "\n")
cat("======================================================================\n")
cat("Script 06 Execution Complete!\n")
