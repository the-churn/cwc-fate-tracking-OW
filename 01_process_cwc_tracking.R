# ==============================================================================
# Script Name:  01_process_cwc_tracking.R
# Description:  Processes TagLab annotations (2015 vs 2022), standardizes 
#               taxa/geometries, tracks cold-water coral (CWC) genets, calculates
#               demographic fates, and performs perimeter-based MDC error propagation.
# Dependencies: tidyverse, stringr
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION & INPUT PATHS
# ------------------------------------------------------------------------------
# Set base path to your local data directory (Replace with your local path)
base_dir <- "D:/PhD_Data(Large)"

# Define subdirectories relative to base path
data_dir   <- file.path(base_dir, "Submission_Dataset")
output_dir <- file.path(base_dir, "Submission_Dataset")

# Define explicit input/output file paths
path_2015    <- file.path(data_dir, "OW15_1fps_full.csv")
path_2022    <- file.path(data_dir, "OW22_1fps_coreg_full.csv")
path_master  <- file.path(output_dir, "master_data_combined.csv")
path_points  <- file.path(output_dir, "Master_Point_Data.csv")
path_tracked <- file.path(output_dir, "master_CWC_tracked.csv")

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Exclusion list for genet.IDs that were not within overlapping area (2022 vs 2015) after visual auditing
exclude_genet_ids <- c(
  941, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1008, 1005, 1004, 
  1006, 1007, 984, 971, 975, 973, 974, 966, 965, 1086, 989, 993, 
  994, 995, 1115
)

# Genet.IDs confirmed present in 2015 that were occluded (not annotatable) in
# the 2022 survey imagery/mesh -- e.g. overgrowth by another colony/sponge/algae,
# sediment drape, structural collapse/rubble, or a camera-geometry shadow --
# after visual auditing. These are within the overlapping survey area (unlike
# exclude_genet_ids above) but must NOT be scored as mortality: they get a
# distinct "Occluded" fate and are treated as censored, not dead.
occluded_2022_genet_ids <- c(
  46, 67, 292, 371, 447, 469, 470, 497, 498, 499, 500, 501, 502, 507, 599, 828
)

# Spatial uncertainty parameters (in cm)
sigma_scale <- 0.100  # Agisoft RMS error (0.001 m)
sigma_gsd   <- 0.216  # Ground Sampling Distance resolution (2.16 mm) of orthomosaic
sigma_edge  <- sqrt(sigma_scale^2 + sigma_gsd^2) # ~0.238 cm combined edge uncertainty


# ------------------------------------------------------------------------------
# 2. INGESTION & DATASET STANDARDIZATION
# ------------------------------------------------------------------------------
data_2015 <- read_csv(path_2015, show_col_types = FALSE)
data_2022 <- read_csv(path_2022, show_col_types = FALSE)

# Combine datasets and classify geometry types
master_dataset <- bind_rows(
  data_2015 %>% mutate(Year = 2015),
  data_2022 %>% mutate(Year = 2022)
) %>%
  mutate(
    TagLab.Type = case_when(
      TagLab.Class.name %in% c("Madrepora oculata", "D. pertusum", "Coral Recruit") ~ "Region",
      str_detect(TagLab.Class.name, "^Primnoa") ~ "Region",
      TRUE ~ "Points"
    )
  ) %>%
  relocate(Year, .before = everything())

# Export master raw dataset
write_csv(master_dataset, path_master)

# Export point dataset (dropping region-specific spatial metrics)
master_point_data <- master_dataset %>% 
  filter(TagLab.Type == "Points") %>% 
  select(-TagLab.Perimeter, -TagLab.Area, -TagLab.Surf.area)

write_csv(master_point_data, path_points)


# ------------------------------------------------------------------------------
# 3. REGION CLEANING & ID HANDLING
# ------------------------------------------------------------------------------
regions_cleaned <- master_dataset %>%
  filter(TagLab.Type == "Region") %>%
  filter(!is.na(Image.name) & Image.name != "") %>%
  filter(!TagLab.Genet.Id %in% exclude_genet_ids) %>%
  mutate(
    # Assign persistent Join_ID across sampling years
    Join_ID = ifelse(
      TagLab.Genet.Id == 0 | is.na(TagLab.Genet.Id), 
      paste0(Year, "_Single_", row_number()), 
      as.character(TagLab.Genet.Id)
    )
  )

# ------------------------------------------------------------------------------
# 4. ANNUAL COHORT AGGREGATION
# ------------------------------------------------------------------------------
agg_2015 <- regions_cleaned %>%
  filter(Year == 2015) %>%
  group_by(Join_ID) %>%
  summarise(
    Species         = first(TagLab.Class.name),
    Area_2015       = sum(TagLab.Area, na.rm = TRUE),
    Surf_Area_2015  = sum(TagLab.Surf.area, na.rm = TRUE),
    Perimeter_2015  = sum(TagLab.Perimeter, na.rm = TRUE), 
    Centroid_x_2015 = mean(TagLab.Centroid.x, na.rm = TRUE),
    Centroid_y_2015 = mean(TagLab.Centroid.y, na.rm = TRUE),
    Poly_Count_2015 = n(),
    .groups         = "drop"
  )

agg_2022 <- regions_cleaned %>%
  filter(Year == 2022) %>%
  group_by(Join_ID) %>%
  summarise(
    Species         = first(TagLab.Class.name),
    Area_2022       = sum(TagLab.Area, na.rm = TRUE),
    Surf_Area_2022  = sum(TagLab.Surf.area, na.rm = TRUE),
    Perimeter_2022  = sum(TagLab.Perimeter, na.rm = TRUE), 
    Centroid_x_2022 = mean(TagLab.Centroid.x, na.rm = TRUE),
    Centroid_y_2022 = mean(TagLab.Centroid.y, na.rm = TRUE),
    Poly_Count_2022 = n(),
    .groups         = "drop"
  )


# ------------------------------------------------------------------------------
# 5. DEMOGRAPHIC FATE ANALYSIS & MDC ERROR PROPAGATION
# ------------------------------------------------------------------------------
master_CWC_flagged <- full_join(agg_2015, agg_2022, by = "Join_ID", suffix = c("_2015", "_2022")) %>%
  mutate(
    Species = coalesce(Species_2015, Species_2022),
    c15     = replace_na(Poly_Count_2015, 0),
    c22     = replace_na(Poly_Count_2022, 0),
    
    # Demographic fate classification (7-year window)
    # NOTE: the occlusion check runs FIRST. case_when() returns the value of
    # the first TRUE condition, so this must precede "Total Mortality" or it
    # will never fire -- occluded genets would silently fall through to
    # "Total Mortality" exactly as before.
    fate = case_when(
      Join_ID %in% as.character(occluded_2022_genet_ids) & c15 > 0 & c22 == 0 ~ "Occluded",
      c15 == 0 & c22 > 0  ~ "Recruit",
      c15 > 0  & c22 == 0 ~ "Total Mortality",
      c15 == 1 & c22 > 1  ~ "Fission",
      c15 > 1  & c22 == 1 ~ "Fusion",
      c15 > 1  & c22 > 1  ~ "Complex",
      c15 == 1 & c22 == 1 ~ "Survivor",
      TRUE                ~ "Unknown"
    ),
    
    # Fill missing spatial values for non-overlapping years.
    # Exception: for "Occluded" genets, a 2022 value of 0 would misrepresent
    # "not measured" as "measured and found empty" -- keep these NA so any
    # arithmetic on them (change, ratio, growth rate, MDC) propagates NA
    # automatically instead of silently fabricating a number, whether that
    # arithmetic happens here or in a downstream script (e.g. Script 02) that
    # reads these columns directly.
    Area_2015      = replace_na(Area_2015, 0),
    Area_2022      = ifelse(fate == "Occluded", NA_real_, replace_na(Area_2022, 0)),
    Surf_Area_2015 = replace_na(Surf_Area_2015, 0),
    Surf_Area_2022 = ifelse(fate == "Occluded", NA_real_, replace_na(Surf_Area_2022, 0)),
    Perimeter_2015 = replace_na(Perimeter_2015, 0),
    Perimeter_2022 = ifelse(fate == "Occluded", NA_real_, replace_na(Perimeter_2022, 0)),
    Centroid_x     = coalesce(Centroid_x_2015, Centroid_x_2022),
    Centroid_y     = coalesce(Centroid_y_2015, Centroid_y_2022),
    
    # Growth metrics (Annualized over 7 years)
    # For "Occluded" genets, Area_2022 = 0 is a non-detection artifact, not a
    # real measurement -- computing change/RGR from it would fabricate a false
    # "shrank to nothing" signal. Force these to NA instead of a number.
    Annual_Area_Change      = ifelse(fate == "Occluded", NA_real_, (Area_2022 - Area_2015) / 7),
    Annual_Surf_Area_Change = ifelse(fate == "Occluded", NA_real_, (Surf_Area_2022 - Surf_Area_2015) / 7),
    
    RGR_Planar  = ifelse(fate != "Occluded" & Area_2015 > 0 & Area_2022 > 0,
                          (log(Area_2022) - log(Area_2015)) / 7, NA),
    RGR_Surface = ifelse(fate != "Occluded" & Surf_Area_2015 > 0 & Surf_Area_2022 > 0,
                          (log(Surf_Area_2022) - log(Surf_Area_2015)) / 7, NA),
    
    # Minimum Detectable Change (MDC) perimeter propagation
    Area_MDC_Total  = 1.96 * sigma_edge * sqrt(Perimeter_2015^2 + Perimeter_2022^2),
    Area_MDC_Annual = Area_MDC_Total / 7,
    
    # Detectability classification
    # Occluded genets get their own label rather than falling into
    # "Detectable Change" (via the old Total Mortality branch) or silently
    # resolving to "Below Detection Limit" once Annual_Area_Change is NA.
    Status_Planar = case_when(
      fate == "Occluded"                        ~ "Not Assessed (Occluded)",
      fate %in% c("Recruit", "Total Mortality")  ~ "Detectable Change",
      abs(Annual_Area_Change) > Area_MDC_Annual  ~ "Detectable Change",
      TRUE                                        ~ "Below Detection Limit"
    )
  )

# ------------------------------------------------------------------------------
# 5b. VALIDATION -- confirm the occlusion override matched what you expect
# ------------------------------------------------------------------------------
occlusion_audit <- master_CWC_flagged %>%
  filter(Join_ID %in% as.character(occluded_2022_genet_ids)) %>%
  select(Join_ID, c15, c22, fate)

unmatched_pattern <- occlusion_audit %>% filter(fate != "Occluded")
if (nrow(unmatched_pattern) > 0) {
  warning(
    "These occluded_2022_genet_ids did NOT receive the 'Occluded' fate -- ",
    "check for a typo, or a genet that actually had a 2022 detection ",
    "(c22 > 0) and shouldn't be on this list: ",
    paste(unmatched_pattern$Join_ID, collapse = ", ")
  )
}

missing_ids <- setdiff(as.character(occluded_2022_genet_ids), occlusion_audit$Join_ID)
if (length(missing_ids) > 0) {
  warning(
    "These occluded_2022_genet_ids were not found in master_CWC_flagged at all -- ",
    "check exclude_genet_ids, ID typos, or whether they were filtered upstream: ",
    paste(missing_ids, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 5c. FINAL COLUMN SELECTION
# ------------------------------------------------------------------------------
master_CWC_tracked <- master_CWC_flagged %>%
  select(
    TagLab.Genet.Id = Join_ID, Species, fate, Status_Planar, 
    Area_MDC_Annual, Area_MDC_Total, Centroid_x, Centroid_y,
    Area_2015, Area_2022, Annual_Area_Change, RGR_Planar,
    Surf_Area_2015, Surf_Area_2022, Annual_Surf_Area_Change, RGR_Surface,
    Perimeter_2015, Perimeter_2022
  )

# Export processed CWC tracking dataset
write_csv(master_CWC_tracked, path_tracked)

cat("Pipeline complete. Processed dataset saved to:\n", path_tracked, "\n")
cat(nrow(occlusion_audit), "occluded genets reclassified; see occlusion_audit for details.\n")
