# ==============================================================================
# Script Name:  02_error_propagation_and_MDC.R
# Description:  Models empirical TagLab measurement error using repeat CNN 
#                annotations (CNN Error.xlsx), classifies ecological state 
#                transitions (persistence, mortality, recruitment, occluded), 
#                propagates dynamic measurement error to calculate 95% Minimum 
#                Detectable Change (MDC95), and exports 
#                Master_CWC_Tracked_MDC.xlsx for downstream plotting and 
#                statistical analyses (Script 03).
# Dependencies: readxl, readr, dplyr, writexl
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LOAD LIBRARIES
# ------------------------------------------------------------------------------
library(readxl)
library(readr)
library(dplyr)
library(writexl)

# ------------------------------------------------------------------------------
# 2. CONFIGURATION & FILE PATHS
# ------------------------------------------------------------------------------
base_dir <- "D:/PhD_Data(Large)/Submission_Dataset"

cnn_error_path    <- file.path(base_dir, "CNN Error.xlsx")
input_csv_path    <- file.path(base_dir, "master_CWC_tracked.csv")
output_xlsx_path  <- file.path(base_dir, "Master_CWC_Tracked_MDC.xlsx")
summary_xlsx_path <- file.path(base_dir, "Master_CWC_Species_Summary.xlsx")

# Sampling duration (2015 to 2022)
time_interval_years <- 7

# ------------------------------------------------------------------------------
# 3. MODEL EMPIRICAL MEASUREMENT ERROR (TAGLAB REPEAT SEGMENTATION)
# ------------------------------------------------------------------------------
cat("Loading CNN error calibration data...\n")
cnn_df <- read_excel(cnn_error_path)

# Standardize lower-case column names from raw repeat segmentation file
colony_stats <- cnn_df %>%
  group_by(colony_id, species, timepoint) %>%
  summarise(
    mean_area = mean(area, na.rm = TRUE),
    sd_area   = sd(area, na.rm = TRUE),
    cv_area   = sd_area / mean_area,
    .groups   = "drop"
  )

# Fit linear empirical error model: SD(area) ~ Mean(area) + Species
model_final <- lm(sd_area ~ mean_area + species, data = colony_stats)

cat("Empirical Error Model Summary:\n")
print(summary(model_final))

# Extract training species factor levels
valid_species_levels <- levels(factor(colony_stats$species))

# ------------------------------------------------------------------------------
# 4. LOAD TRACKED DATASET & CLASSIFY STATE TRANSITIONS
# ------------------------------------------------------------------------------
cat("\nProcessing tracked master dataset from Script 01...\n")
df <- read_csv(input_csv_path, show_col_types = FALSE)

# Guard: state_transition below is derived FROM Script 01's `fate` column,
# not re-derived from raw area. Re-deriving from area alone
# (present_2022 = Area_2022 > 0) is exactly what silently reclassified
# occluded genets as "mortality" before -- it can't distinguish "not
# measured" from "measured and gone." If `fate` is missing, this is almost
# certainly a stale CSV from before the occlusion fix.
if (!"fate" %in% names(df)) {
  stop("`fate` column not found in ", input_csv_path,
       " -- re-run the updated Script 01 before this script.")
}

# Map initial/final area and establish presence/absence rules
df <- df %>%
  mutate(
    A1 = Area_2015,
    A2 = Area_2022,
    
    present_2015 = !is.na(A1) & A1 > 0,
    present_2022 = !is.na(A2) & A2 > 0,
    
    # Classify state transitions FROM Script 01's `fate`, with a new
    # "occluded" category alongside the naming Script 03 already expects.
    # CHECK SCRIPT 03 before trusting its output for these genets: a
    # hardcoded factor(levels = c("persistence","mortality","recruitment"))
    # or a state_transition %in% c(...) filter that doesn't include
    # "occluded" will silently drop or NA them out rather than error.
    state_transition = case_when(
      fate == "Occluded"                                      ~ "occluded",
      fate == "Total Mortality"                                ~ "mortality",
      fate == "Recruit"                                        ~ "recruitment",
      fate %in% c("Survivor", "Fission", "Fusion", "Complex")  ~ "persistence",
      TRUE                                                     ~ "unknown"
    )
  )

# ------------------------------------------------------------------------------
# 5. DYNAMIC ERROR PROPAGATION & DETECTABLE CHANGE (MDC95)
# ------------------------------------------------------------------------------
# Harmonize species mapping for model prediction:
# 1. Map 'Primnoa msp.1' -> 'Primnoa msp.5' (or matching Primnoa training level)
# 2. Set 'Coral Recruit' -> NA (recruits are excluded from MDC predictions)
target_primnoa <- if ("Primnoa msp.5" %in% valid_species_levels) "Primnoa msp.5" else grep("Primnoa", valid_species_levels, value = TRUE)[1]

df_pred_mapped <- df %>%
  mutate(
    species_for_pred = case_when(
      Species == "Primnoa msp.1" ~ target_primnoa,
      Species == "Coral Recruit" ~ NA_character_,
      TRUE                       ~ Species
    ),
    species_for_pred = factor(species_for_pred, levels = valid_species_levels)
  )

# Build prediction dataframes using df_pred_mapped$species_for_pred
pred_df_2015 <- data.frame(
  mean_area = df_pred_mapped$A1, 
  species   = df_pred_mapped$species_for_pred
)

pred_df_2022 <- data.frame(
  mean_area = df_pred_mapped$A2, 
  species   = df_pred_mapped$species_for_pred
)

df <- df %>%
  mutate(
    # Predict dynamic measurement uncertainty (sigma) per colony based on
    # area & species. For occluded genets, A2 is NA (given Script 01's
    # occlusion patch), so sigma2 predicts to NA automatically. The explicit
    # `state_transition == "occluded"` guards below are a second line of
    # defense in case this ever runs against an older tracked CSV where A2
    # was still coded as 0 instead of NA.
    sigma1 = predict(model_final, newdata = pred_df_2015),
    sigma2 = predict(model_final, newdata = pred_df_2022),
    
    # Total change in planar area over 7-year monitoring window (cm²)
    change_area = A2 - A1,
    
    # Root Sum of Squares (quadrature sum) for error propagation.
    # Occluded genets have no valid 2022 measurement, so there's no change to
    # propagate error for.
    sigma_change = ifelse(Species == "Coral Recruit" | state_transition == "occluded",
                           NA, sqrt(sigma1^2 + sigma2^2)),
    
    # Absolute 95% Minimum Detectable Change threshold (cm²)
    MDC95 = ifelse(Species == "Coral Recruit" | state_transition == "occluded",
                    NA, 1.96 * sigma_change),
    
    # Binary detectability boolean (evaluated for persistent colonies, excluding recruits)
    detectable_change = case_when(
      Species == "Coral Recruit"        ~ NA,
      state_transition == "occluded"    ~ NA,
      state_transition == "persistence" ~ abs(change_area) > MDC95,
      TRUE                               ~ NA
    ),
    
    # ANNUALIZED METRIC ALIGNMENT FOR SCRIPT 03 ---------------------------------
    Annual_Area_Change = ifelse(state_transition == "occluded", NA_real_, change_area / time_interval_years),
    MDC95_Annual        = ifelse(state_transition == "occluded", NA_real_, MDC95 / time_interval_years)
  )

# ------------------------------------------------------------------------------
# 6. SPECIES-LEVEL ECOLOGICAL SUMMARY
# ------------------------------------------------------------------------------
species_summary <- df %>%
  group_by(Species) %>%
  summarise(
    n                  = n(),
    mean_change        = mean(change_area[state_transition == "persistence"], na.rm = TRUE),
    mean_MDC           = mean(MDC95[state_transition == "persistence"], na.rm = TRUE),
    percent_detectable = mean(detectable_change == TRUE, na.rm = TRUE) * 100,
    # NOTE: denominator is total tracked n (including occluded genets), so
    # mortality_rate/recruitment_rate read as "share of the full 2015
    # cohort," not "share of genets with a resolved 2022 fate."
    # occlusion_rate is reported separately so occluded genets stay visible
    # rather than being folded into either mortality or persistence. If you
    # want mortality relative to only the resolved cohort, divide by
    # (n - occluded_n) instead.
    mortality_rate     = sum(state_transition == "mortality") / n,
    recruitment_rate   = sum(state_transition == "recruitment") / n,
    occluded_n         = sum(state_transition == "occluded"),
    occlusion_rate     = sum(state_transition == "occluded") / n,
    .groups            = "drop"
  )

cat("\nSpecies Demographic & Error Summary:\n")
print(species_summary)

# ------------------------------------------------------------------------------
# 7. EXPORT DATASETS
# ------------------------------------------------------------------------------
dir.create(dirname(output_xlsx_path), recursive = TRUE, showWarnings = FALSE)

write_xlsx(df, output_xlsx_path)
write_xlsx(species_summary, summary_xlsx_path)

cat("\n======================================================================\n")
cat("Script 02 execution complete!\n")
cat("Master MDC file saved to:\n ", output_xlsx_path, "\n")
cat("Species Summary saved to:\n ", summary_xlsx_path, "\n")
cat("======================================================================\n")
